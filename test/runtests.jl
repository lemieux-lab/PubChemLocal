using PubChemLocal
using DataFrames
using Test

# Minimal V2000 mol blocks (real, RDKit-parseable) for two distinct
# compounds, formatted as PubChem-style SDF records so the ingestion code
# is exercised against real structures rather than mocked-out identifiers.

const METHANOL = """
methanol
     RDKit          2D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
> <PUBCHEM_COMPOUND_CID>
887

> <PUBCHEM_IUPAC_NAME>
methanol

\$\$\$\$
"""

const ETHANOL = """
ethanol
     RDKit          2D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END
> <PUBCHEM_COMPOUND_CID>
702

> <PUBCHEM_IUPAC_NAME>
ethanol

\$\$\$\$
"""

@testset "PubChemLocal" begin

    @testset "sdf: split_records / parse_tags" begin
        text = METHANOL * ETHANOL
        records = PubChemLocal.split_records(text)
        @test length(records) == 2

        tags = PubChemLocal.parse_tags(records[1])
        @test tags["PUBCHEM_COMPOUND_CID"] == ["887"]
        @test tags["PUBCHEM_IUPAC_NAME"] == ["methanol"]

        # multi-line / multi-valued tag (e.g. PUBCHEM_SUBSTANCE_SYNONYM) must
        # keep every value, not just the first line of the block
        multi = "mol block placeholder\n\n> <PUBCHEM_SUBSTANCE_SYNONYM>\nfoo\nbar\nbaz\n"
        tags2 = PubChemLocal.parse_tags(multi)
        @test tags2["PUBCHEM_SUBSTANCE_SYNONYM"] == ["foo", "bar", "baz"]
    end

    @testset "identifiers: computed columns" begin
        ids = identifiers_from_molblock(METHANOL)
        @test !isempty(ids.inchi)
        @test !isempty(ids.inchikey)
        @test ids.mkey == ids.inchikey[1:14]

        # methanol and ethanol must not collide
        ids2 = identifiers_from_molblock(ETHANOL)
        @test ids.inchikey != ids2.inchikey
        @test ids.mkey != ids2.mkey

        # smiles path should agree with the molblock path for the same compound
        ids3 = identifiers_from_smiles("CO")
        @test ids3.inchikey == ids.inchikey

        bad = identifiers_from_smiles("not a smiles")
        @test bad.inchi == "" && bad.inchikey == "" && bad.mkey == ""
    end

    @testset "compounds table" begin
        records = PubChemLocal.split_records(METHANOL * ETHANOL)
        df = parse_compounds(records)
        @test nrow(df) == 2
        @test Set(df[!, "PUBCHEM_COMPOUND_CID"]) == Set(["887", "702"])
        @test all(!isempty, df.inchikey)
        @test length(unique(df.mkey)) == 2
    end

    @testset "substances tables + identify" begin
        sub1 = """
        substance placeholder

        \$\$\$\$
        > <PUBCHEM_SUBSTANCE_ID>
        1001

        > <PUBCHEM_CID>
        887

        > <PUBCHEM_EXT_DATASOURCE_NAME>
        TestVendor

        > <PUBCHEM_EXT_DATASOURCE_REGID>
        CAT-887

        """
        sub2 = """
        substance placeholder

        > <PUBCHEM_SUBSTANCE_ID>
        1002

        > <PUBCHEM_CID>
        702

        > <PUBCHEM_EXT_DATASOURCE_NAME>
        TestVendor

        > <PUBCHEM_EXT_DATASOURCE_REGID>
        CAT-702

        """
        records = [sub1, sub2]
        subs = parse_substances(records)
        @test nrow(subs.substances) == 2
        @test Set(subs.substances.cid) == Set(["887", "702"])
        @test nrow(subs.identifiers) == 2
        @test Set(subs.identifiers.value) == Set(["CAT-887", "CAT-702"])

        compounds = parse_compounds(PubChemLocal.split_records(METHANOL * ETHANOL))

        index = build_identifier_index(subs.identifiers, subs.substances, compounds)
        hit = identify("CAT-887", index)
        @test nrow(hit) == 1
        @test hit[1, "PUBCHEM_COMPOUND_CID"] == "887"

        miss = identify("no-such-id", index)
        @test nrow(miss) == 0

        batch = identify(["CAT-887", "CAT-702", "no-such-id"], index)
        @test nrow(batch) == 2
        @test Set(batch.queried_id) == Set(["CAT-887", "CAT-702"])
    end

end
