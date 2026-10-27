# Pilot investigation: issue #32

Issue: [#32](https://github.com/CostaLab/scMEGA/issues/32). Read the full report and
its one maintainer comment. Related #29 uses obsolete argument names (confirmed from its screenshot), a
distinct documentation/signature mismatch; #27 concerns split-layer selection in a different function.

## Finding and reproduction

**Confirmed bug, P0 scientific impact.** With SeuratObject 5.4.0, `LayerData` on an
Assay5 selects the requested cells but returns them in storage order. Its installed
`LayerData.StdAssay` implementation sorts matched indices. `CreatePairedObject`
assumed request order and unconditionally renamed those columns to `cell_name`.
The resulting RNA/ATAC counts can therefore belong to the wrong declared pair,
while RNA metadata/embeddings follow df.pair ordering. No warning detects that
reassignment. The report is not merely an API signature incompatibility.

`Rscript scripts/maintenance/reproduce_32.R` uses a 2-by-12 synthetic sparse counts
matrix with distinguishable columns. Requesting c12,...,c1 returns c1,...,c12;
renaming produces a different matrix from named selection. The script asserts
both the original failure and the correction's count preservation with the real
accessor. It does not require scMEGA installation or biological downloads.

The regression suite was also run against pristine `e6ce4f9:R/pair_cells.R`, placed
in /tmp/scmega-original/R. This produced 14 failed expectations and three test
errors in the five then-present test cases (two expected dense-conversion warnings,
no skips after Signac installation). Real ChromatinAssay construction was exercised;
the original code failed count invariants in full and partial pairing scenarios.
See [before log](evidence/issue-32-before-full.log). The initial pre-Signac run had
one integration skip, recorded separately; it is not claimed as full validation.

## Scientific invariants and contract

1. For each row in df.pair, RNA counts come from its RNA identifier and ATAC counts
   from its ATAC identifier, regardless of either object's column order.
2. Renaming output pairs does not change their original count vectors.
3. Gene/peak names and their vectors stay associated; sparse dimensions survive.
4. Output assays, metadata and embeddings have consistent pair IDs; genomic ranges
   remain correct. Existing Signac `min.cells=10` filtering is preserved.
5. Only declared pairs are retained. Unlisted source cells are permitted; referenced
   missing cells fail explicitly. This function does not infer an intersection or
   require RNA and ATAC barcodes to be identical.
6. Repeated source IDs are supported intentionally: PairCells/GetPairList can produce
   many-to-one matches. Output cell_name IDs must be unique. Valid Seurat inputs
   enforce source cell-name integrity; no new assay coercion is performed.

## Minimal correction

After extraction, verify all referenced names exist and explicitly index both
matrices with their respective df.pair columns using `drop = FALSE`, then rename.
Validate required pairing columns, non-empty/non-missing named identifiers and
unique output names. Factor identifiers are converted to character; numeric
positional identifiers are rejected. These errors prevent silent mispairing.

The public signature, defaults, normalization, matching algorithm, reference
annotations, peak filtering and data representation are unchanged. Correctly
ordered inputs preserve the same results. Previously misassigned counts change;
analyses built from affected paired objects should be regenerated and downstream
scientific conclusions rechecked. This is a scientific correction, not numerical
cleanup. The relevant roxygen/Rd argument documentation and example were repaired;
README and other scientific documentation were preserved.

## Tests and validation scope

Final focused tests cover legacy Assay/Assay5, dense/sparse input, identical and
permuted modality order, partial overlap, unused cells, absent referenced cells,
repeated source IDs, singleton extraction, factor IDs, invalid numeric/missing/empty
identifiers, duplicate output IDs, and count/feature integrity. Real Signac output
is tested using both Assay5 and ChromatinAssay ATAC sources, full and 10-of-12 pairs,
RNA metadata, embeddings and genomic ranges. The 10-cell subset respects the existing
min.cells=10 threshold; a full one-cell ChromatinAssay is not advertised as supported.

**117 expectations pass**, six test cases, zero failures/errors/skips; two recorded
SeuratObject warnings convert dense fixture matrices to dgCMatrix. R 4.5.3,
SeuratObject 5.4.0, Signac 1.17.1, Matrix 1.7-5, testthat 3.3.2, GenomicRanges 1.62.1.
See [after log](evidence/issue-32-after.log) and [session](evidence/pilot-session.log).

These are focused source tests, not installed-scMEGA tests. A child environment
supplies the actual Signac::CreateChromatinAssay constructor because the package has
a pre-existing missing Signac import. Matrix/Seurat accessor dispatch, Signac object
construction and metadata/reduction operations are real. The extraction tests
replace only the output constructor to inspect both matrices before construction.
The package's namespace/dependency defects remain visible and were not patched
around in installation/check.

`R CMD build --no-build-vignettes --no-manual .` passes. Installation and built-
tarball check fail on 12 unavailable Imports; check also reports undeclared Seurat.
No claim of passing full package CI, v4 support, alternate dependency stacks or
real-dataset validation is made. Split counts-layer joining remains outside this
fix; missing requested cells in the selected layer fail instead of misaligning.

The actual pilot diff was separately validated/applied to a clean local checkout
of infrastructure commit a9bea1e, then the trusted regression suite passed all
117 expectations again. [Patch identity](evidence/pilot-patch-validation.json) and
[round-trip log](evidence/pilot-patch-roundtrip.log) record this check.

## Review artifact

Infrastructure is on local branch `codex/maintenance-infrastructure`; the pilot is
on `codex/issue-32-createpairedobject`, based on that infrastructure. A local draft
PR description is [pilot-pr.md](pilot-pr.md). No remote branch/PR was created because
remote publication was not authorized. Package-wide checks must be repaired and
rerun before approval; automated publication correctly remains blocked.
