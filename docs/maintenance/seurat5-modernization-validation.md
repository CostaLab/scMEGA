# Seurat v5 modernization validation

## Scope and dependency

Branch: `codex/seurat5-memory-modernization`, based on `3b9e904`. The repository
already declares `Seurat (>= 5.0.0)` and `SeuratObject (>= 5.0.0)`; these correct
minimums are retained. This PR is independent of Draft PR #46 and does not duplicate
its delayed-scaling optimization. Both touch `CoembedData()`, so whichever merges
second must update from `main` and rerun the exact fixture comparison.

## APIs modernized

- one validated helper resolves a requested Seurat v5 layer and rejects ambiguous
  split layers;
- `GetTrajectory()`, `CreatePairedObject()`, and `CoembedData()` use the resolved
  layer rather than allowing first-match behavior;
- direct metadata/reduction slots are replaced with `object[[]]`, `Embeddings()`,
  and `Reductions()`;
- direct Signac motif slots are replaced with `GetMotifData()` plus explicit
  ChromatinAssay/motif validation;
- the stale `GetAssayData` namespace import is removed;
- stable Seurat/Signac release installation replaces old development-branch
  instructions.

No Seurat v4 version branch existed to remove. Legacy `Assay` support is retained
because it is valid in SeuratObject v5 and required for existing data. Signac
`ChromatinAssay` remains intact; it is never coerced to Assay5 or blindly joined.

## Split-layer contract

When more than one matching `counts.*` or `data.*` layer exists, scMEGA now lists
the layers and asks the caller to join only the intended disjoint layers. It does
not normalize again, choose one batch, or guess whether measurements are compatible.
Single layers, explicitly joined Assay5 layers, legacy Assay, and ChromatinAssay
are supported.

## Scientific equivalence

The deterministic Cardiomyocyte CoembedData fixture is exactly identical to
baseline for cells, metadata, assays, variable features, GeneActivity scale data,
RNA data, PCA, and UMAP. Maximum numerical difference is zero. Layer tests verify
feature/cell order and matrix values across all three supported assay classes.
Public motif access returns the same identifiers as the former slot access.

No scientific algorithm, model, threshold, default, cell, feature, peak, pairing,
coordinate, or regulatory interaction changes.

Validation results:

- 198 testthat expectations passed; 0 failed, skipped, or warned;
- 22 maintenance security tests passed;
- package installation and `R CMD build` passed;
- ordinary `R CMD check`: 0 errors, 0 warnings, 2 pre-existing notes;
- deterministic real-data CoembedData output is exact at every saved boundary.

## BPCells and sketching

BPCells was genuinely installed and evaluated, but not adopted because normalized
and trajectory outputs were not exactly equivalent and Signac/legacy constructors
materialized disk-backed counts. It does not solve dense CCA. Sketching was reviewed
but not executed because it subsamples cells and violates the complete-analysis
contract. See `seurat5-compatibility.md` and `evaluate_bpcells.R`.

## Full tutorial status

The current exact Seurat CCA is unchanged, so repeating the known memory-pressure
failure would be unsafe and provide no new evidence. The latest controlled full
attempt is the #46 run: stopped during `FindTransferAnchors()` at 33.80 GiB RSS and
20.23 GiB swap, before anchors returned. The full Cardiomyocyte tutorial remains
blocked on a 64 GB host. Subsampling is not introduced or recommended as an exact
solution.

## Review recommendation

The API/layer modernization is safe for focused review after the complete package
suite and hosted checks pass. BPCells and implicit/blockwise CCA require separate
research PRs only if exact output equivalence can be demonstrated. PR #46 is the
only current memory-reduction proposal and is green but modest; this PR does not
claim to solve the tutorial.
