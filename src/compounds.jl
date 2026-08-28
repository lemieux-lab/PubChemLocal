#=
Build the `compounds` table from PubChem's bulk Compound SDF dump
(ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/*.sdf.gz).

One row per CID: every SDF property tag found on the record, plus the
three computed columns (inchi, inchikey, mkey) from identifiers.jl. A few
tags are multi-line in real data (see parse_compounds's docstring), their
lines are newline-joined rather than truncated.

This is a straight refactor of the original PubChem.jl script, same
threaded block-processing / checkpoint-file / resumable-block design, just
reorganized into functions that take their paths as arguments instead of a
hardcoded, host-specific constant.
=#

using DataFrames, Arrow, SHA

## Parsing ##

"""
    parse_compounds(records) -> DataFrame

`records` is a vector of raw SDF record strings (as produced by
[`split_records`](@ref)). Returns one row per record: every tag found via
[`parse_tags`](@ref) plus `inchi`, `inchikey`, `mkey` computed from the
record's mol block.

Most Compound tags are single-valued, but not all: `PUBCHEM_BONDANNOTATIONS`,
`PUBCHEM_COORDINATE_TYPE` and `PUBCHEM_NONSTANDARDBOND` are confirmed
multi-line in real data (2026-08-26 sample), and there may be others. Since
this table is wide (one column per tag), a genuinely multi-valued tag's
lines are newline-joined into that one cell rather than truncated to the
first line. Full content is kept, just not exploded into its own row/table
the way `cid_links`/`identifiers` are for `substances` (none of these three
feed `inchi`/`inchikey`/`mkey`, which come from RDKit on the mol block, not
from these tags, and they're not identifiers anyone searches compounds by).
"""
function parse_compounds(records::AbstractVector{<:AbstractString})
    n = length(records)
    props = Dict{String,Vector{String}}()
    inchi = Vector{String}(undef, n)
    inchikey = Vector{String}(undef, n)
    mkey = Vector{String}(undef, n)

    for (i, record) in enumerate(records)
        for (tag, values) in parse_tags(record)
            col = get!(() -> fill("", n), props, tag)
            col[i] = join(values, '\n')
        end
        ids = identifiers_from_molblock(record)
        inchi[i] = ids.inchi
        inchikey[i] = ids.inchikey
        mkey[i] = ids.mkey
    end

    df = DataFrame(props)
    df.inchi = inchi
    df.inchikey = inchikey
    df.mkey = mkey
    return df
end

## Per-file processing ##

"""
    compounds_from_sdf(fn; blk_size=20_000) -> DataFrame

Decompress and parse a single `.sdf.gz` Compound bulk file, splitting its
records into chunks of `blk_size` processed on separate tasks
(`Threads.@spawn`) and concatenating the results (`cols=:union`, since a
chunk boundary can still land such that one chunk sees a tag another
doesn't).
"""
function compounds_from_sdf(fn::AbstractString; blk_size::Integer=20_000)
    records = split_records(read_gz_text(fn))
    blocks = chunk(records, blk_size)
    dfs = Vector{DataFrame}(undef, length(blocks))

    @sync for (i, blk) in enumerate(blocks)
        Threads.@spawn dfs[i] = parse_compounds(blk)
    end

    df = DataFrame()
    for tdf in dfs
        append!(df, tdf, cols=:union)
    end
    return df
end

## Full-corpus build ##

_checksum(sub_fn) = bytes2hex(sha256(join(sub_fn, ':')))[1:7]
_out_fn(out_dir, i, sub_fn) = joinpath(out_dir, "compounds-$(lpad(i, 2, '0'))-$(_checksum(sub_fn)).arrow")

"""
    build_compounds_table(data_dir; out_dir=data_dir, file_blk_size=50,
                           record_blk_size=20_000, block_tasks=2, file_tasks=4)

Parse every `.sdf.gz` file in `data_dir` into the `compounds` table and
write it out as one Arrow file per group of `file_blk_size` input files
(default 50) under `out_dir`.

Resumable: a group is skipped if its (checksum-named) output Arrow file
already exists, and each input file gets a `<file>.ongoing` touch-file
while it's being processed, so a killed/resumed run can tell "not started"
from "in progress" from "done", same scheme as the original script. This
does not clean up stale `.ongoing` files from a run that died mid-file.
Check for those by hand before resuming after a crash.

Needs to run on a host that can see `data_dir` (PubChem's raw Compound
dumps currently live under `/scratch/lemieuxs/pubchem/compounds`, which is
only mounted on one node of the cluster).
"""
function build_compounds_table(data_dir::AbstractString; out_dir::AbstractString=data_dir,
                                file_blk_size::Integer=50, record_blk_size::Integer=20_000,
                                block_tasks::Integer=2, file_tasks::Integer=4)
    all_fn = filter(endswith(".gz"), readdir(data_dir))
    file_blocks = chunk(all_fn, file_blk_size)

    ch_blk = Channel{Int}(length(file_blocks))
    foreach(i -> put!(ch_blk, i), eachindex(file_blocks))
    close(ch_blk)

    Threads.foreach(ch_blk; ntasks=block_tasks) do blk_i
        sub_fn = file_blocks[blk_i]
        blk_fn = _out_fn(out_dir, blk_i, sub_fn)
        if isfile(blk_fn)
            @info "Skipping block $blk_i of size $(length(sub_fn)), output already exists." blk_fn
            return
        end
        @info "Starting block $blk_i of size $(length(sub_fn))..."

        dfs = Vector{DataFrame}(undef, length(sub_fn))
        ch_file = Channel{Int}(length(sub_fn))
        foreach(i -> put!(ch_file, i), eachindex(sub_fn))
        close(ch_file)

        Threads.foreach(ch_file; ntasks=file_tasks) do i
            fn = joinpath(data_dir, sub_fn[i])
            @info "Starting file $(sub_fn[i]) [$i/$(length(sub_fn))]..."
            touch("$fn.ongoing")
            dfs[i] = compounds_from_sdf(fn; blk_size=record_blk_size)
            rm("$fn.ongoing")
            @info "Done with file $(sub_fn[i]) [$i/$(length(sub_fn))]." rss_mb=rss_mb()
        end

        df = DataFrame()
        for tdf in dfs
            append!(df, tdf, cols=:union)
        end

        Arrow.write("$blk_fn.tmp", df)
        mv("$blk_fn.tmp", blk_fn)
        @info "Done with block $blk_i." blk_fn rss_mb=rss_mb()
    end
    return nothing
end
