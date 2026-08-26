# The actual point of the package: given an identifier as you'd hand it to
# PubChem's PUG REST `xref` endpoint (a vendor catalog number, a CAS number,
# ...), find the matching compound(s) locally instead of over the network.
#
# COMPOUND_CID_TAG below is confirmed against a real Compound .sdf.gz record
# (2026-08-26).
#
# Scale note: at full PubChem size (~118M compounds, ~300M+ substances)
# repeatedly `filter`-ing the raw tables per query is not viable. The index
# built here (Dict value -> cids, plus grouping compounds by cid) turns
# repeated identify() calls into O(1) dict lookups; it still requires both
# full tables to be materialized in memory, which may itself become the
# next bottleneck once running against the real bulk tables — if so, the
# natural next step is pushing this down into DuckDB (reads the Arrow
# files directly, no separate load step) rather than rebuilding the index
# in-process each session.

using DataFrames

const COMPOUND_CID_TAG = "PUBCHEM_COMPOUND_CID"

"""
    IdentifierIndex

Precomputed `identifiers -> cid_links -> compounds` lookup structure, built
once with [`build_identifier_index`](@ref) and then queried repeatedly via
[`identify`](@ref) without re-scanning the source tables each time.
"""
struct IdentifierIndex
    value_to_cids::Dict{String,Vector{String}}   # identifier value -> CID(s)
    compounds_by_cid::GroupedDataFrame            # grouped `compounds`, by cid_tag
    cid_tag::String
end

"""
    build_identifier_index(identifiers, cid_links, compounds;
                            tags=DEFAULT_IDENTIFIER_TAGS, cid_tag=COMPOUND_CID_TAG)
        -> IdentifierIndex

- `identifiers`: long-format table from [`parse_substances`](@ref)/
  [`build_substances_tables`](@ref) (columns `sid`, `tag`, `value`).
- `cid_links`: that same function's other output (columns `sid`, `cid`,
  `assoc_type`), carries the SID -> CID link(s) `identifiers` doesn't. A
  substance can link to more than one CID (see [`CID_ASSOCIATIONS_TAG`](@ref));
  all of them are used here, `assoc_type` isn't filtered on.
- `compounds`: the wide table from [`parse_compounds`](@ref)/
  [`build_compounds_table`](@ref).
- `tags`: which identifier tags to index (default: just the registry-ID tag
  Carl's PubChem-API lookups used). Must be a subset of the tags the
  `identifiers` table was actually built with.
"""
function build_identifier_index(identifiers::AbstractDataFrame, cid_links::AbstractDataFrame,
                                 compounds::AbstractDataFrame;
                                 tags::AbstractVector{<:AbstractString}=DEFAULT_IDENTIFIER_TAGS,
                                 cid_tag::AbstractString=COMPOUND_CID_TAG)
    tagset = Set(tags)
    sid_to_cids = Dict{String,Vector{String}}()
    for row in eachrow(cid_links)
        isempty(row.cid) && continue
        push!(get!(() -> String[], sid_to_cids, row.sid), row.cid)
    end

    value_to_cids = Dict{String,Vector{String}}()
    for row in eachrow(identifiers)
        row.tag in tagset || continue
        cids = get(sid_to_cids, row.sid, String[])
        isempty(cids) && continue
        dest = get!(() -> String[], value_to_cids, row.value)
        for cid in cids
            push!(dest, cid)
        end
    end

    compounds_by_cid = groupby(compounds, cid_tag)
    return IdentifierIndex(value_to_cids, compounds_by_cid, cid_tag)
end

"""
    identify(id, index::IdentifierIndex) -> DataFrame

Resolve one external `id` (e.g. a vendor catalog number or CAS number) to
compound rows, entirely against the locally-built, precomputed `index` —
the local replacement for `xref/RegistryID/<id>/...` against the PubChem
REST API. Returns an empty frame if `id` isn't found; possibly more than
one row if it resolves ambiguously to several CIDs.
"""
function identify(id::AbstractString, index::IdentifierIndex)
    cids = get(index.value_to_cids, id, String[])
    isempty(cids) && return DataFrame()

    frames = DataFrame[]
    for cid in unique(cids)
        key = (cid,)
        haskey(index.compounds_by_cid, key) || continue
        push!(frames, DataFrame(index.compounds_by_cid[key]))
    end
    isempty(frames) && return DataFrame()
    return reduce(vcat, frames; cols=:union)
end

"""
    identify(ids, index::IdentifierIndex) -> DataFrame

Batch form: resolves each id in `ids`, returning the union of matches with
an extra `queried_id` column recording which input id produced each row
(the mapping isn't guaranteed one-to-one in either direction).
"""
function identify(ids::AbstractVector{<:AbstractString}, index::IdentifierIndex)
    out = DataFrame()
    for id in ids
        matched = identify(id, index)
        isempty(matched) && continue
        matched = copy(matched)
        matched.queried_id = fill(id, nrow(matched))
        append!(out, matched; cols=:union)
    end
    return out
end
