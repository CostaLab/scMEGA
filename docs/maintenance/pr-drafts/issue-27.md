# fix: validate coembedding layers and retain legacy reference features (#27)

Related: #27. Local draft; no remote PR.

Split normalized RNA layers can supply only the first batch to TransferData. Report ambiguous/missing data layers before anchor calculation, requiring the caller to deliberately prepare the reference. A tested preparation using JoinLayers preserves named counts and completes the pipeline.

Legacy Assay coembedding also fails because JoinLayers has no legacy method and merging loses the supplied reference variable-feature set. Join only Assay5 output and restore the existing reference selection only when the merged object has none. Do not choose new variable features or automatically combine ambiguous measurements. Document normalized RNA and selected-feature prerequisites.

The real synthetic Seurat/Signac pipeline passes for both assay classes, preserving observed RNA data, ATAC counts, cell identities and the supplied variable-feature selection. Original source fails the regressions. Final: 17 expectations passed, five upstream warnings, no failures/skips. Installed tests also pass with separate metadata proposal ec24fdf in a scratch copy. Install/build pass; R CMD check retains six existing warnings and three notes. Original plant data and old development versions were not validated.

Branch codex/issue-27-coembed, commit 7fbcce6, based on main e6ce4f9. No #25/#32 commits are included. Integrate and validate #32 separately before release; no normalization, matching, PCA/UMAP defaults or public signature changes are included here.
