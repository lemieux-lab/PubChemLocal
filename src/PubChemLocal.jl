module PubChemLocal

include("sdf.jl")
include("identifiers.jl")
include("compounds.jl")
include("substances.jl")
include("lookup.jl")
include("structure.jl")

# sdf.jl
export read_gz_text, split_records, parse_tags

# identifiers.jl
export MolIdentifiers, identifiers_from_molblock, identifiers_from_smiles, identifiers_from_inchi

# compounds.jl
export parse_compounds, compounds_from_sdf, build_compounds_table

# substances.jl
export parse_substances, substances_from_sdf, build_substances_tables, DEFAULT_IDENTIFIER_TAGS

# lookup.jl
export IdentifierIndex, build_identifier_index, identify

# structure.jl
export StructureIndex, build_structure_index, identify_by_structure,
       identify_by_smiles, identify_by_molblock

end # module PubChemLocal
