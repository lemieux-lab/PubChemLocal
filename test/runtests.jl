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

        # a multi-line tag (e.g. PUBCHEM_BONDANNOTATIONS in real data) must be
        # kept in full, newline-joined, not truncated to its first line
        multi = """
        one atom placeholder
             RDKit          2D

          1  0  0  0  0  0  0  0  0  0999 V2000
            0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
        M  END
        > <PUBCHEM_COMPOUND_CID>
        1

        > <PUBCHEM_BONDANNOTATIONS>
        8 9 6
        8 10 6

        \$\$\$\$
        """
        mdf = parse_compounds(PubChemLocal.split_records(multi))
        @test mdf[1, "PUBCHEM_BONDANNOTATIONS"] == "8 9 6\n8 10 6"
    end

    @testset "substances tables + identify" begin
        sub1 = """
        substance placeholder

        > <PUBCHEM_SUBSTANCE_ID>
        1001

        > <PUBCHEM_CID_ASSOCIATIONS>
        887  1

        > <PUBCHEM_EXT_DATASOURCE_NAME>
        TestVendor

        > <PUBCHEM_EXT_DATASOURCE_REGID>
        CAT-887

        """
        sub2 = """
        substance placeholder

        > <PUBCHEM_SUBSTANCE_ID>
        1002

        > <PUBCHEM_CID_ASSOCIATIONS>
        702  1

        > <PUBCHEM_EXT_DATASOURCE_NAME>
        TestVendor

        > <PUBCHEM_EXT_DATASOURCE_REGID>
        CAT-702

        """
        sub3 = """
        substance placeholder

        > <PUBCHEM_SUBSTANCE_ID>
        1003

        > <PUBCHEM_CID_ASSOCIATIONS>
        887  1
        702  2

        > <PUBCHEM_EXT_DATASOURCE_NAME>
        TestVendor

        > <PUBCHEM_EXT_DATASOURCE_REGID>
        CAT-BOTH

        """
        records = [sub1, sub2, sub3]
        subs = parse_substances(records)
        @test nrow(subs.substances) == 3
        @test nrow(subs.cid_links) == 4  # sub1: 1, sub2: 1, sub3: 2
        @test Set(subs.cid_links[subs.cid_links.sid .== "1003", :cid]) == Set(["887", "702"])
        @test nrow(subs.identifiers) == 3
        @test Set(subs.identifiers.value) == Set(["CAT-887", "CAT-702", "CAT-BOTH"])

        compounds = parse_compounds(PubChemLocal.split_records(METHANOL * ETHANOL))

        index = build_identifier_index(subs.identifiers, subs.cid_links, compounds)
        hit = identify("CAT-887", index)
        @test nrow(hit) == 1
        @test hit[1, "PUBCHEM_COMPOUND_CID"] == "887"

        both = identify("CAT-BOTH", index)  # one substance linked to two CIDs
        @test nrow(both) == 2
        @test Set(both[!, "PUBCHEM_COMPOUND_CID"]) == Set(["887", "702"])

        miss = identify("no-such-id", index)
        @test nrow(miss) == 0

        batch = identify(["CAT-887", "CAT-702", "no-such-id"], index)
        @test nrow(batch) == 2
        @test Set(batch.queried_id) == Set(["CAT-887", "CAT-702"])
    end

    @testset "structure lookup" begin
        compounds = parse_compounds(PubChemLocal.split_records(METHANOL * ETHANOL))
        index = build_structure_index(compounds)

        hit = identify_by_smiles("CO", index)  # methanol
        @test nrow(hit) == 1
        @test hit[1, "PUBCHEM_COMPOUND_CID"] == "887"
        @test hit[1, :match_tier] == "inchikey"

        hit2 = identify_by_molblock(ETHANOL, index)
        @test nrow(hit2) == 1
        @test hit2[1, "PUBCHEM_COMPOUND_CID"] == "702"

        miss = identify_by_smiles("c1ccccc1", index)  # benzene: not in this tiny table
        @test nrow(miss) == 0

        # methanol/ethanol have no stereoisomer to exercise a real
        # inchikey-miss/mkey-hit fallback, so drive that path against a
        # hand-built table instead
        fake = DataFrame(cid=["1"], inchikey=["AAAAAAAAAAAAAA-BBBBBBBBBB-C"], mkey=["AAAAAAAAAAAAAA"])
        fake_index = StructureIndex(groupby(fake, :inchikey), groupby(fake, :mkey))
        query = MolIdentifiers("inchi", "AAAAAAAAAAAAAA-ZZZZZZZZZZ-D", "AAAAAAAAAAAAAA")

        @test nrow(identify_by_structure(query, fake_index; tier=:inchikey)) == 0
        loose = identify_by_structure(query, fake_index)  # :both -> falls back to mkey
        @test nrow(loose) == 1
        @test loose[1, :match_tier] == "mkey"

        # an unresolvable query must never match an unresolvable compound
        # just because both have inchikey/mkey == ""
        unresolved = DataFrame(cid=["2"], inchikey=[""], mkey=[""])
        unresolved_index = StructureIndex(groupby(unresolved, :inchikey), groupby(unresolved, :mkey))
        @test nrow(identify_by_structure(PubChemLocal.EMPTY_IDENTIFIERS, unresolved_index)) == 0
    end

end
