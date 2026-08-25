# Build the `substances` tables from PubChem's bulk Substance SDF dump
# (ftp.ncbi.nlm.nih.gov/pubchem/Substance/CURRENT-Full/SDF/*.sdf.gz).
#
# Substance records are the as-deposited-by-source records: they don't carry
# a structure worth recomputing inchi/inchikey/mkey from (that's what the
# `compounds` table is for), but they do carry the heterogeneous external
# identifiers (vendor catalog numbers, CAS numbers, source names, ...) that
# people actually search by, plus a link to the standardized CID(s) PubChem
# resolved the substance to. This is the local stand-in for what
# https://pubchem.ncbi.nlm.nih.gov/rest/pug/.../xref/RegistryID/<id>/...
# does over the network (see the old prep_cigs.jl `pubchem_smiles`).
#
# *** TAG NAMES BELOW ARE UNVERIFIED — confirm against a real
# *** Substance .sdf.gz on the cluster node that can see
# *** /scratch/lemieuxs/pubchem/ before trusting any ingested output.
# They're PubChem's documented/conventional field names, but this repo has
# no local sample to check them against.

using DataFrames, Arrow, SHA

const SID_TAG = "PUBCHEM_SUBSTANCE_ID"
const CID_TAG = "PUBCHEM_CID"                        # CID(s) this substance standardizes to
const SOURCE_NAME_TAG = "PUBCHEM_EXT_DATASOURCE_NAME"

"""
Tags indexed into the `identifiers` table by default. `PUBCHEM_EXT_DATASOURCE_REGID`
is the source's own ID for the substance (catalog number, CAS number, etc.)
— it's what PubChem's `xref/RegistryID` endpoint matches against, i.e. what
Carl was querying through the API. Pass a broader/narrower `identifier_tags`
to [`parse_substances`](@ref)/[`build_substances_tables`](@ref) to change
this (e.g. add `"PUBCHEM_SUBSTANCE_SYNONYM"` to also index synonyms).
"""
const DEFAULT_IDENTIFIER_TAGS = ["PUBCHEM_EXT_DATASOURCE_REGID"]

"""
    parse_substances(records; identifier_tags=DEFAULT_IDENTIFIER_TAGS)
        -> (substances=DataFrame, identifiers=DataFrame)

`records` is a vector of raw SDF record strings. Returns two tables sharing
`sid` as a key:

- `substances`: one row per SID — `sid`, `cid` (empty string if the
  substance didn't standardize to a compound), `source`.
- `identifiers`: long format, one row per (sid, tag, value) triple, only for
  tags in `identifier_tags` — a substance can have several values for the
  same tag (e.g. more than one synonym), each gets its own row. This is the
  table [`identify`](@ref) searches.

Kept long/flexible rather than one-wide-column-per-tag (as `compounds` does)
because Substance tags are heterogeneous per source and often multi-valued;
a wide table would be mostly missing values and would need a schema change
every time a new source's tag shows up.
"""
function parse_substances(records::AbstractVector{<:AbstractString};
                           identifier_tags::AbstractVector{<:AbstractString}=DEFAULT_IDENTIFIER_TAGS)
    sid = String[]
    cid = String[]
    source = String[]

    id_sid = String[]
    id_tag = String[]
    id_value = String[]

    tagset = Set(identifier_tags)

    for record in records
        tags = parse_tags(record)
        this_sid = haskey(tags, SID_TAG) ? first(tags[SID_TAG]) : ""
        push!(sid, this_sid)
        push!(cid, haskey(tags, CID_TAG) ? first(tags[CID_TAG]) : "")
        push!(source, haskey(tags, SOURCE_NAME_TAG) ? first(tags[SOURCE_NAME_TAG]) : "")

        for tag in identifier_tags
            haskey(tags, tag) || continue
            for value in tags[tag]
                push!(id_sid, this_sid)
                push!(id_tag, tag)
                push!(id_value, value)
            end
        end
    end

    substances = DataFrame(sid=sid, cid=cid, source=source)
    identifiers = DataFrame(sid=id_sid, tag=id_tag, value=id_value)
    return (substances=substances, identifiers=identifiers)
end

"""
    substances_from_sdf(fn; blk_size=20_000, identifier_tags=DEFAULT_IDENTIFIER_TAGS)
        -> (substances=DataFrame, identifiers=DataFrame)

Substance-side equivalent of [`compounds_from_sdf`](@ref).
"""
function substances_from_sdf(fn::AbstractString; blk_size::Integer=20_000,
                              identifier_tags::AbstractVector{<:AbstractString}=DEFAULT_IDENTIFIER_TAGS)
    records = split_records(read_gz_text(fn))
    blocks = chunk(records, blk_size)
    results = Vector{NamedTuple{(:substances, :identifiers),Tuple{DataFrame,DataFrame}}}(undef, length(blocks))

    @sync for (i, blk) in enumerate(blocks)
        Threads.@spawn results[i] = parse_substances(blk; identifier_tags)
    end

    substances = reduce(vcat, (r.substances for r in results); cols=:union)
    identifiers = reduce(vcat, (r.identifiers for r in results); cols=:union)
    return (substances=substances, identifiers=identifiers)
end

_sub_out_fn(out_dir, prefix, i, sub_fn) =
    joinpath(out_dir, "$(prefix)-$(lpad(i, 2, '0'))-$(bytes2hex(sha256(join(sub_fn, ':')))[1:7]).arrow")

"""
    build_substances_tables(data_dir; out_dir=data_dir, identifier_tags=DEFAULT_IDENTIFIER_TAGS,
                             file_blk_size=50, record_blk_size=20_000, block_tasks=2, file_tasks=4)

Substance-side equivalent of [`build_compounds_table`](@ref): parses every
`.sdf.gz` in `data_dir` and writes, per group of `file_blk_size` input
files, a `substances-NN-<hash>.arrow` and matching `identifiers-NN-<hash>.arrow`
under `out_dir`. Same resumability/checkpoint scheme (skip if the substances
output already exists, `.ongoing` touch-files per input file).
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
        id_out = _sub_out_fn(out_dir, "identifiers", blk_i, sub_fn)
        if isfile(sub_out)
            @info "Skipping block $blk_i of size $(length(sub_fn)), output already exists." sub_out
            return
        end
        @info "Starting block $blk_i of size $(length(sub_fn))..."

        results = Vector{NamedTuple{(:substances, :identifiers),Tuple{DataFrame,DataFrame}}}(undef, length(sub_fn))
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
        identifiers = reduce(vcat, (r.identifiers for r in results); cols=:union)

        Arrow.write("$sub_out.tmp", substances)
        Arrow.write("$id_out.tmp", identifiers)
        mv("$sub_out.tmp", sub_out)
        mv("$id_out.tmp", id_out)
        @info "Done with block $blk_i." sub_out id_out
    end
    return nothing
end
