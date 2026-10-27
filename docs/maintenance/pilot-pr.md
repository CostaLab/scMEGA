# fix: resolve CreatePairedObject issue (#32)

Local draft only. No GitHub PR has been published.

## Related Issue
Fixes #32 (the reproduced pairing/count-order defect; full package validation remains blocked below).

## Problem
Assay5 extraction can return cells in storage order. Renaming them immediately to
paired-cell identifiers assigns RNA or ATAC counts to the wrong declared pair.

## Root Cause
CreatePairedObject treated LayerData's cell selection as an ordering guarantee.
The installed SeuratObject implementation sorts matched cell indices, while RNA
metadata and embeddings are explicitly selected in the pairing-table order.

## Changes
Reindex both count matrices by their modality-specific cell names before renaming,
with drop=FALSE. Reject missing/invalid pairing identifiers and duplicate output
names; preserve intentional repeated source cells. Correct the function's argument
documentation and example without changing its signature or defaults.

## Regression Tests
Synthetic dense/sparse, legacy Assay/Assay5, identical/permuted/partially overlapping
cells, repeated sources, singleton extraction, factors and invalid IDs. Real Signac
construction verifies counts, feature identities, genomic ranges, metadata and
embeddings for both Assay5 and ChromatinAssay sources. Original source fails the
new count invariants; the corrected implementation passes.

## Validation
- Focused unit/integration tests: 117 expectations passed, zero failures/errors/skips;
  two SeuratObject dense-input conversion warnings recorded.
- Package build: passed, without vignettes/manual.
- Package installation: blocked by 12 missing Imports.
- Built-tarball R CMD check: 1 ERROR, missing dependencies and undeclared Seurat.
- Compatibility: exercised only on the recorded R 4.5.3/SeuratObject 5.4.0/
  Signac 1.17.1/Matrix 1.7-5 environment; no broad supported matrix claimed.

## Scientific Impact
Corrects silent assignment of count vectors to paired cells. Affected paired objects
and downstream analyses should be regenerated. Algorithms, normalization, public
defaults and Signac peak filtering are unchanged.

## Limitations
Focused source tests supply the real Signac constructor in a child environment to
isolate the existing missing namespace import. Installed-package validation and
full CI remain blocked; resolve metadata/dependencies through separate review.
Split-layer merging, real biological datasets and additional dependency versions
were not validated. This draft is not evidence of a successful production rollout.
