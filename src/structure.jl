#=
The other direction from lookup.jl: given a structure (SMILES or a mol
block) rather than an external identifier, find the matching row(s) in the
local `compounds` table directly, no substances/identifiers involved.
Mirrors the two-tier matching from the old intersect.jl script: strict
match on `inchikey` (the main merge key), falling back to the looser
`mkey` (same atoms/connectivity/charge, ignoring stereochemistry and
isotopes) when nothing matches exactly.
=#

using DataFrames

## Index ##

"""
    StructureIndex

Precomputed `compounds` groupings by `inchikey` and by `mkey`, built once
with [`build_structure_index`](@ref) and then queried repeatedly via
[`identify_by_structure`](@ref)/[`identify_by_smiles`](@ref)/
[`identify_by_molblock`](@ref) without regrouping each call.

# Fields
- `by_inchikey`: `compounds` grouped by `inchikey`, for the strict tier.
- `by_mkey`: `compounds` grouped by `mkey`, for the loose tier.
"""
struct StructureIndex
    by_inchikey::GroupedDataFrame
    by_mkey::GroupedDataFrame
end

"""
    build_structure_index(compounds) -> StructureIndex

`compounds` is the wide table from [`parse_compounds`](@ref)/
[`build_compounds_table`](@ref) (must have `inchikey` and `mkey` columns,
which it always does).
"""
build_structure_index(compounds::AbstractDataFrame) =
    StructureIndex(groupby(compounds, :inchikey), groupby(compounds, :mkey))

function _group_lookup(gdf::GroupedDataFrame, key::AbstractString)
    key == "" && return nothing
    haskey(gdf, (key,)) || return nothing
    return DataFrame(gdf[(key,)])
end

## Lookup ##

"""
    identify_by_structure(ids::MolIdentifiers, index::StructureIndex; tier=:both) -> DataFrame

Resolve computed identifiers to `compounds` rows. `tier`:

- `:inchikey`: exact match only (same compound, including stereochemistry).
- `:mkey`: loose match only (same "graph": atoms, connectivity, charge).
- `:both` (default): try `:inchikey` first, if nothing matches, fall back
  to `:mkey`.

Returns an empty frame if `ids` couldn't be computed in the first place
(e.g. RDKit failed to parse the input) or nothing matches. The result
carries an extra `match_tier` column (`"inchikey"` or `"mkey"`) recording
which layer actually matched, since a `:both` match is silent about that
otherwise.
"""
function identify_by_structure(ids::MolIdentifiers, index::StructureIndex; tier::Symbol=:both)
    tier in (:inchikey, :mkey, :both) || throw(ArgumentError("tier must be :inchikey, :mkey or :both"))

    if tier in (:inchikey, :both)
        hit = _group_lookup(index.by_inchikey, ids.inchikey)
        if hit !== nothing
            hit.match_tier = fill("inchikey", nrow(hit))
            return hit
        end
        tier == :inchikey && return DataFrame()
    end

    hit = _group_lookup(index.by_mkey, ids.mkey)
    if hit !== nothing
        hit.match_tier = fill("mkey", nrow(hit))
        return hit
    end
    return DataFrame()
end

"""
    identify_by_smiles(smiles, index::StructureIndex; tier=:both) -> DataFrame
    identify_by_molblock(molblock, index::StructureIndex; tier=:both) -> DataFrame

Convenience wrappers computing identifiers from a SMILES string or an SDF
mol block (via [`identifiers_from_smiles`](@ref)/
[`identifiers_from_molblock`](@ref)) and resolving them against `index`.
"""
identify_by_smiles(smiles::AbstractString, index::StructureIndex; tier::Symbol=:both) =
    identify_by_structure(identifiers_from_smiles(smiles), index; tier)

identify_by_molblock(molblock::AbstractString, index::StructureIndex; tier::Symbol=:both) =
    identify_by_structure(identifiers_from_molblock(molblock), index; tier)
