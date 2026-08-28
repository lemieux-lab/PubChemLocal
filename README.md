# PubChemLocal.jl

Local replacement for the PubChem PUG REST `xref` workflow: build `compounds`
and `substances` tables once from PubChem's bulk `.sdf.gz` dumps, then
resolve identifiers to compound structures entirely offline, in either
direction: from an external identifier (vendor catalog number, CAS number,
...) to a compound, or from a structure (SMILES/mol block) to the matching
compound(s) already in `compounds`.

Status: early scaffold. Tag names are confirmed against real data
(2026-08-26): one real Compound and one real Substance record by hand, plus
a ~5000-record multi-valued-tag sweep on each (using this package's own
`parse_tags`) to check which tags are genuinely multi-line rather than
assuming. That sweep is what caught `PUBCHEM_CID_ASSOCIATIONS` (fixed, see
`cid_links` below) and `PUBCHEM_BONDANNOTATIONS`/`PUBCHEM_COORDINATE_TYPE`/
`PUBCHEM_NONSTANDARDBOND` (fixed, `parse_compounds` now joins multi-line
values instead of truncating). `SID_TAG`/`SOURCE_NAME_TAG` swept clean
(single-valued, as assumed). The builders still haven't run against a real
bulk file end to end, only synthetic fixtures (see `test/runtests.jl`). See
"Known gaps" below before trusting output from a real run.

## Layout

- `src/sdf.jl`: shared low-level SDF parsing (gunzip, record splitting,
  `> <TAG>` property block extraction).
- `src/identifiers.jl`: the three computed columns (`inchi`, `inchikey`,
  `mkey`) from a mol block or a SMILES string. Single source of truth,
  replacing three independent, slightly-drifted implementations that used to
  live in `PubChem.jl` and `intersect.jl`.
- `src/compounds.jl`: `compounds` table, one row per CID, from PubChem's
  Compound bulk dump. Refactor of the original `PubChem.jl` script (same
  threaded block-processing / resumable-block-file design), now exposed as
  functions instead of a hardcoded top-level script. A handful of tags are
  multi-line in real data; their lines are newline-joined into that column
  rather than truncated to the first line.
- `src/substances.jl`: `substances`/`cid_links`/`identifiers` tables, from
  PubChem's Substance bulk dump. New, no equivalent existed before.
  `substances` is one row per SID. `cid_links` is long-format (`sid`, `cid`,
  `assoc_type`), the SID→CID link, since a substance can associate to more
  than one CID. `identifiers` is long-format (`sid`, `tag`, `value`) so it
  stays flexible across sources' heterogeneous identifier tags without a
  schema change per source.
- `src/lookup.jl`: `identify` walks `identifiers → cid_links → compounds`
  given an external id to return matching compound rows.
  Local stand-in for `xref/RegistryID/<id>/...` against the PubChem API.
- `src/structure.jl`: the other direction. `identify_by_smiles`/
  `identify_by_molblock` compute `inchikey`/`mkey` for a query structure and
  match it against `compounds` directly, strict (`inchikey`) first, falling
  back to loose (`mkey`) if nothing matches exactly, the same two-tier
  strategy `intersect.jl` used by hand.

## Usage

```julia
using PubChemLocal, Arrow, DataFrames

# one-time setup, run on the cluster node (rerun whenever the bulk dump changes)
# that can see the raw .sdf.gz files, see scripts/build_compounds.jl and
# scripts/build_substances.jl
build_compounds_table("/scratch/lemieuxs/pubchem/compounds")
build_substances_tables("/scratch/lemieuxs/pubchem/substances")

# per-session: load the Arrow outputs and build the lookup index once
arrow_files(dir, prefix) = filter(f -> startswith(f, prefix), readdir(dir; join=true))
load_table(dir, prefix) = vcat(DataFrame.(Arrow.Table.(arrow_files(dir, prefix)))...; cols=:union)

compounds = load_table(out_dir, "compounds-")
cid_links = load_table(out_dir, "cid_links-")
identifiers = load_table(out_dir, "identifiers-")
index = build_identifier_index(identifiers, cid_links, compounds)

identify("CAT-887", index)                       # -> compound row(s), or empty
identify(["CAT-887", "CAS-123-45-6"], index)      # batch form

# the other direction: structure -> compound, no substances table involved
sindex = build_structure_index(compounds)
identify_by_smiles("CCO", sindex)                 # strict inchikey match, falls back to mkey
identify_by_smiles("CCO", sindex; tier=:inchikey) # strict only
identify_by_molblock(some_molblock, sindex; tier=:mkey)  # loose only
```

To search by more than the default identifier tag (see
`DEFAULT_IDENTIFIER_TAGS`), pass `identifier_tags` when building the
substances tables (it controls what gets indexed into the long
`identifiers` table in the first place), then the matching `tags` kwarg to
`build_identifier_index`.

## Known gaps / next decisions

1. `assoc_type` in `cid_links` isn't decoded. `PUBCHEM_CID_ASSOCIATIONS`
   gives each substance a `<cid> <type>` pair per line (a real example:
   `15685509  1`), and `parse_substances` captures every pair as-is rather
   than assuming a single "the" CID per substance. What the numeric type
   codes mean (same-connectivity vs. mixture-component vs. same-stereo,
   etc.) hasn't been pinned down, so `identify()` currently treats every
   association as an equally valid link. If a real build over-matches (e.g.
   a mixture's component associations pulling in unrelated CIDs for one
   identifier), that's the place to add a filter, once the type codes are
   figured out from a broader sample or PubChem's own documentation.
2. Scale. `build_identifier_index`/`build_structure_index` load full
   `compounds`/`substances`/`identifiers` tables into memory and build
   in-process `Dict`/`groupby` indexes. Fine for a subset, but at full
   PubChem size (~118M compounds, ~300M+ substances) this may need to move
   to something that doesn't require materializing everything per session,
   e.g. DuckDB querying the Arrow files directly. Deferred until we know the
   actual working-set size after a real build.