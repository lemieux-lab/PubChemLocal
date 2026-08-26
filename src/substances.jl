# Build the `substances` tables from PubChem's bulk Substance SDF dump
# (ftp.ncbi.nlm.nih.gov/pubchem/Substance/CURRENT-Full/SDF/*.sdf.gz).
#
# Substance records are the as-deposited-by-source records: they don't carry
# a structure worth recomputing inchi/inchikey/mkey from (that's what the
# `compounds` table is for), but they do carry the heterogeneous external
# identifiers (vendor catalog numbers, CAS numbers, source names, ...) that
# people actually search by, plus a link to the CID(s) PubChem resolved the
# substance to. This is the local stand-in for what
# https://pubchem.ncbi.nlm.nih.gov/rest/pug/.../xref/RegistryID/<id>/...
# does over the network (see the old prep_cigs.jl `pubchem_smiles`).
#
# Tag names below confirmed against a real Substance .sdf.gz on the cluster
# node (2026-08-26): SID_TAG, SOURCE_NAME_TAG and the REGID entry in
# DEFAULT_IDENTIFIER_TAGS all matched as guessed. CID_ASSOCIATIONS_TAG did
# not — there is no flat "PUBCHEM_CID" tag; see its docstring below.

using DataFrames, Arrow, SHA

const SID_TAG = "PUBCHEM_SUBSTANCE_ID"
const SOURCE_NAME_TAG = "PUBCHEM_EXT_DATASOURCE_NAME"

"""
The substance-to-compound link. Unlike a simple `PUBCHEM_CID` field, real
records hold one `<cid> <assoc_type>` pair per line, e.g.:

```
> <PUBCHEM_CID_ASSOCIATIONS>
15685509  1
```

— and can hold more than one line/CID per substance. `assoc_type` is a
PubChem-internal numeric code (same-connectivity vs. mixture-component vs.
same-stereo, etc.) whose exact mapping isn't decoded here; every association
is captured as-is in the `cid_links` table (see [`parse_substances`](@ref))
rather than guessing which type code means "the" canonical CID. If a
real build turns out to over-match (e.g. a mixture's mixture-component
associations pulling in unrelated CIDs), filter `cid_links` by
`assoc_type` before building the lookup index — once the type codes are
pinned down.
"""
const CID_ASSOCIATIONS_TAG = "PUBCHEM_CID_ASSOCIATIONS"

"""
Tags indexed into the `identifiers` table by default. `PUBCHEM_EXT_DATASOURCE_REGID`
is the source's own ID for the substance (catalog number, CAS number, etc.),
it's what PubChem's `xref/RegistryID` endpoint matches against, i.e. what
Carl was querying through the API. Pass a broader/narrower `identifier_tags`
to [`parse_substances`](@ref)/[`build_substances_tables`](@ref) to change
this (e.g. add `"PUBCHEM_SUBSTANCE_SYNONYM"` to also index synonyms).
"""
const DEFAULT_IDENTIFIER_TAGS = ["PUBCHEM_EXT_DATASOURCE_REGID"]

"""
    parse_substances(records; identifier_tags=DEFAULT_IDENTIFIER_TAGS)
        -> (substances=DataFrame, cid_links=DataFrame, identifiers=DataFrame)

`records` is a vector of raw SDF record strings. Returns three tables
sharing `sid` as a key:

- `substances`: one row per SID — `sid`, `source`.
- `cid_links`: long format, one row per (sid, cid, assoc_type) triple, from
  [`CID_ASSOCIATIONS_TAG`](@ref) — a substance can link to more than one CID.
- `identifiers`: long format, one row per (sid, tag, value) triple, only for
  tags in `identifier_tags` — a substance can have several values for the
  same tag (e.g. more than one synonym), each gets its own row. This is the
  table [`identify`](@ref) searches.

Both `cid_links` and `identifiers` are kept long/flexible rather than
one-wide-column-per-tag (as `compounds` does) because both are genuinely
multi-valued per substance and heterogeneous per source; a wide table would
be mostly missing values and would need a schema change every time a new
source's tag shows up.
"""
function parse_substances(records::AbstractVector{<:AbstractString};
                           identifier_tags::AbstractVector{<:AbstractString}=DEFAULT_IDENTIFIER_TAGS)
    sid = String[]
    source = String[]

    link_sid = String[]
    link_cid = String[]
    link_type = String[]

    id_sid = String[]
    id_tag = String[]
    id_value = String[]

    for record in records
        tags = parse_tags(record)
        this_sid = haskey(tags, SID_TAG) ? first(tags[SID_TAG]) : ""
        push!(sid, this_sid)
        push!(source, haskey(tags, SOURCE_NAME_TAG) ? first(tags[SOURCE_NAME_TAG]) : "")

        for assoc in get(tags, CID_ASSOCIATIONS_TAG, String[])
            fields = split(assoc)
            isempty(fields) && continue
            push!(link_sid, this_sid)
            push!(link_cid, fields[1])
            push!(link_type, length(fields) > 1 ? fields[2] : "")
        end

        for tag in identifier_tags
            haskey(tags, tag) || continue
            for value in tags[tag]
                push!(id_sid, this_sid)
                push!(id_tag, tag)
                push!(id_value, value)
            end
        end
    end

    substances = DataFrame(sid=sid, source=source)
    cid_links = DataFrame(sid=link_sid, cid=link_cid, assoc_type=link_type)
    identifiers = DataFrame(sid=id_sid, tag=id_tag, value=id_value)
    return (substances=substances, cid_links=cid_links, identifiers=identifiers)
end

const _SUBSTANCE_RESULT = NamedTuple{(:substances, :cid_links, :identifiers),Tuple{DataFrame,DataFrame,DataFrame}}

"""
    substances_from_sdf(fn; blk_size=20_000, identifier_tags=DEFAULT_IDENTIFIER_TAGS)
        -> (substances=DataFrame, cid_links=DataFrame, identifiers=DataFrame)

Substance-side equivalent of [`compounds_from_sdf`](@ref).
"""
function substances_from_sdf(fn::AbstractString; blk_size::Integer=20_000,
                              identifier_tags::AbstractVector{<:AbstractString}=DEFAULT_IDENTIFIER_TAGS)
    records = split_records(read_gz_text(fn))
    blocks = chunk(records, blk_size)
    results = Vector{_SUBSTANCE_RESULT}(undef, length(blocks))

    @sync for (i, blk) in enumerate(blocks)
        Threads.@spawn results[i] = parse_substances(blk; identifier_tags)
    end

    substances = reduce(vcat, (r.substances for r in results); cols=:union)
    cid_links = reduce(vcat, (r.cid_links for r in results); cols=:union)
    identifiers = reduce(vcat, (r.identifiers for r in results); cols=:union)
    return (substances=substances, cid_links=cid_links, identifiers=identifiers)
end

_sub_out_fn(out_dir, prefix, i, sub_fn) =
    joinpath(out_dir, "$(prefix)-$(lpad(i, 2, '0'))-$(bytes2hex(sha256(join(sub_fn, ':')))[1:7]).arrow")

"""
    build_substances_tables(data_dir; out_dir=data_dir, identifier_tags=DEFAULT_IDENTIFIER_TAGS,
                             file_blk_size=50, record_blk_size=20_000, block_tasks=2, file_tasks=4)

Substance-side equivalent of [`build_compounds_table`](@ref): parses every
`.sdf.gz` in `data_dir` and writes, per group of `file_blk_size` input
files, `substances-NN-<hash>.arrow`, `cid_links-NN-<hash>.arrow` and
`identifiers-NN-<hash>.arrow` under `out_dir`. Same resumability/checkpoint
scheme (skip if the substances output already exists, `.ongoing`
touch-files per input file).
"""
function build_substances_tables(data_dir::AbstractString; out_dir::AbstractString=data_dir,
                                  identifier_tags::AbstractVector{<:AbstractString}=DEFAULT_IDENTIFIER_TAGS,
                                  file_blk_size::Integer=50, record_blk_size::Integer=20_000,
                                  block_tasks::Integer=2, file_tasks::Integer=4)
    all_fn = filter(endswith(".gz"), readdir(data_dir))
    file_blocks = chunk(all_fn, file_blk_size)

    ch_blk = Channel{Int}(length(file_blocks))
    foreach(i -> put!(ch_blk, i), eachindex(file_blocks))
    close(ch_blk)

    Threads.foreach(ch_blk; ntasks=block_tasks) do blk_i
        sub_fn = file_blocks[blk_i]
        sub_out = _sub_out_fn(out_dir, "substances", blk_i, sub_fn)
        link_out = _sub_out_fn(out_dir, "cid_links", blk_i, sub_fn)
        id_out = _sub_out_fn(out_dir, "identifiers", blk_i, sub_fn)
        if isfile(sub_out)
            @info "Skipping block $blk_i of size $(length(sub_fn)), output already exists." sub_out
            return
        end
        @info "Starting block $blk_i of size $(length(sub_fn))..."

        results = Vector{_SUBSTANCE_RESULT}(undef, length(sub_fn))
        ch_file = Channel{Int}(length(sub_fn))
        foreach(i -> put!(ch_file, i), eachindex(sub_fn))
        close(ch_file)

        Threads.foreach(ch_file; ntasks=file_tasks) do i
            fn = joinpath(data_dir, sub_fn[i])
            @info "Starting file $(sub_fn[i]) [$i/$(length(sub_fn))]..."
            touch("$fn.ongoing")
            results[i] = substances_from_sdf(fn; blk_size=record_blk_size, identifier_tags)
            rm("$fn.ongoing")
            @info "Done with file $(sub_fn[i]) [$i/$(length(sub_fn))]."
        end

        substances = reduce(vcat, (r.substances for r in results); cols=:union)
        cid_links = reduce(vcat, (r.cid_links for r in results); cols=:union)
        identifiers = reduce(vcat, (r.identifiers for r in results); cols=:union)

        Arrow.write("$sub_out.tmp", substances)
        Arrow.write("$link_out.tmp", cid_links)
        Arrow.write("$id_out.tmp", identifiers)
        mv("$sub_out.tmp", sub_out)
        mv("$link_out.tmp", link_out)
        mv("$id_out.tmp", id_out)
        @info "Done with block $blk_i." sub_out link_out id_out
    end
    return nothing
end
