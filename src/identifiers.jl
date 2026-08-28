#=
The three computed columns, in one place.

These used to be reimplemented three separate times across the lab's
scripts (PubChem.jl on mol blocks, intersect.jl on SMILES) and had already
drifted apart slightly. This is the single source of truth going forward.
=#

using RDKitMinimalLib

"""
    MolIdentifiers

# Fields
- `inchi`: full InChI, computed by RDKit (rdkitminimallib).
- `inchikey`: the standard InChIKey, the main merge key across tables.
- `mkey`: the first 14 characters of `inchikey`, i.e. just the hash of the
  InChI *connectivity layer* (atoms, connectivity, charge). A much looser
  key than `inchikey`: two compounds that differ only in stereochemistry or
  isotopes share an `mkey` but not an `inchikey`.

Any field is `""` when it could not be computed (e.g. RDKit failed to parse
the input structure).
"""
struct MolIdentifiers
    inchi::String
    inchikey::String
    mkey::String
end

function Base.show(io::IO, ids::MolIdentifiers)
    print(io, "MolIdentifiers(inchikey=", repr(ids.inchikey), ", mkey=", repr(ids.mkey), ")")
end

const EMPTY_IDENTIFIERS = MolIdentifiers("", "", "")

mkey_from_inchikey(inchikey::AbstractString) =
    length(inchikey) < 14 ? "" : inchikey[1:14]

function identifiers_from_inchi(inchi::AbstractString)
    isempty(inchi) && return EMPTY_IDENTIFIERS
    inchikey = get_inchikey_for_inchi(inchi)
    return MolIdentifiers(inchi, inchikey, mkey_from_inchikey(inchikey))
end

"""
    identifiers_from_molblock(molblock) -> MolIdentifiers

As used for the `compounds` table: `molblock` is the CTAB block from an SDF
record (the record text itself is fine, RDKit reads up to the first
`M  END` and ignores the trailing `> <TAG>` property blocks).
"""
function identifiers_from_molblock(molblock::AbstractString)
    inchi = get_inchi_for_molblock(molblock)
    return identifiers_from_inchi(inchi)
end

"""
    identifiers_from_smiles(smiles) -> MolIdentifiers

As used e.g. for vendor catalog SMILES columns. Returns `EMPTY_IDENTIFIERS`
if RDKit can't parse `smiles`.
"""
function identifiers_from_smiles(smiles::AbstractString)
    mol = get_mol(smiles)
    mol === nothing && return EMPTY_IDENTIFIERS
    return identifiers_from_inchi(get_inchi(mol))
end
