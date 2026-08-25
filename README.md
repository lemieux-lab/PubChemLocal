# PubChemLocal.jl

Local replacement for the PubChem PUG REST `xref` workflow: build `compounds`
and `substances` tables once from PubChem's bulk `.sdf.gz` dumps, then
resolve external identifiers (vendor catalog numbers, CAS numbers, ...) to
compound structures/InChIKeys entirely offline.

Status: **early scaffold**, structurally tested against synthetic fixtures
(see `test/runtests.jl`) but not yet run against real PubChem bulk data.
see "Known gaps" below before trusting output from a real run.

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
  functions instead of a hardcoded top-level script.
- `src/substances.jl`: `substances`/`identifiers` tables, from PubChem's
  Substance bulk dump. New, no equivalent existed before. `substances` is
  one row per SID with the SID→CID link. `identifiers` is long-format
  (`sid`, `tag`, `value`) so it stays flexible across sources' heterogeneous
  identifier tags without a schema change per source.
- `src/lookup.jl`: `identify` walks `identifiers → substances → compounds`
  given an external id to return matching compound rows.
  Local stand-in for `xref/RegistryID/<id>/...` against the PubChem API.

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
substances = load_table(out_dir, "substances-")
identifiers = load_table(out_dir, "identifiers-")
index = build_identifier_index(identifiers, substances, compounds)

identify("CAT-887", index)                       # -> compound row(s), or empty
identify(["CAT-887", "CAS-123-45-6"], index)      # batch form
```

To search by more than the default identifier tag (see
`DEFAULT_IDENTIFIER_TAGS`), pass `identifier_tags` when *building* the
substances tables (it controls what gets indexed into the long
`identifiers` table in the first place), then the matching `tags` kwarg to
`build_identifier_index`.

## Known gaps / next decisions

1. **Substance/Compound SDF tag names are unverified.** `SID_TAG`,
   `CID_TAG`, `SOURCE_NAME_TAG` (`substances.jl`), `DEFAULT_IDENTIFIER_TAGS`,
   and `COMPOUND_CID_TAG` (`lookup.jl`) are PubChem's documented/conventional
   field names, but nothing here has been run against a real bulk file yet
   (no sample was reachable from this account, only from the cluster node
   that mounts `/scratch/lemieuxs/pubchem/`). First real-data task: parse
   one Substance `.sdf.gz` record by hand and confirm/correct these
   constants before trusting a full build.
2. **Where do the raw Substance dumps live?** `scripts/build_substances.jl`
   guesses a `substances` dir alongside the existing `compounds` one:
   confirm/adjust.
3. **Scale.** `build_identifier_index` loads full `compounds`/`substances`/
   `identifiers` tables into memory and builds an in-process `Dict` +
   `groupby` index. Fine for a subset, but at full PubChem size (~118M
   compounds, ~300M+ substances) this may need to move to something that
   doesn't require materializing everything per session, e.g. DuckDB
   querying the Arrow files directly. Deferred until we know the actual
   working-set size after a real build.
4. **`test.jld2`** in the parent directory (2.2 GB) reports zero JLD2 keys
   when opened, looks like a stale/incomplete artifact from earlier
   exploration, not consumed by anything here.
