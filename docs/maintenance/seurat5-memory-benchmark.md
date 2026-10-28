# Seurat v5 modernization benchmark

Environment: macOS arm64, R 4.5.3, Seurat 5.5.0, SeuratObject 5.4.0,
Signac 1.17.1, Matrix 1.7-5. Baseline and modernized runs use seed 1234,
identical ordering, parameters, and numerical libraries.

| Dataset | Cells | Baseline RSS | Modernized RSS | Reduction | Results equivalent |
|---|---:|---:|---:|---:|---|
| Small layer fixture | 20 | not separately measurable | not separately measurable | none claimed | Exact for Assay5, Assay, ChromatinAssay |
| Medium real CoembedData fixture | 1,800 | 1,904,787,456 B | 1,907,638,272 B | -0.15% (noise) | Exact |
| Full Cardiomyocyte | 83,020 before merging | historical 46.31 GB RSS / 96.46 GB footprint | not rerun: CCA code unchanged | none claimed | Full output unavailable |

Medium runtime changed from 4.391 to 4.817 seconds in single runs. This is not
treated as a regression or improvement without repetitions; the source changes do
not alter numerical operations.

Every medium-fixture component is `identical()`: cells, metadata, assay list,
variable features, GeneActivity scale layer, RNA data, PCA, and UMAP. Maximum
absolute numerical difference is zero. Existing pairing/trajectory/GRN regression
tests provide downstream structural and numerical coverage.

## BPCells evaluation

| Operation | Representation | Result |
|---|---|---|
| Raw counts and subsetting | BPCells `MatrixDir`/`RenameDims` | exact values and names |
| NormalizeData | BPCells queued matrix | maximum difference `4.82422e-7` from dgCMatrix |
| GetTrajectory on normalized layer | BPCells queued matrix | maximum difference `1.838324e-5` |
| Legacy Assay construction | `dgCMatrix` | disk backing materialized |
| Signac ChromatinAssay construction | `dgCMatrix` | disk backing materialized |
| Exact CCA | dense matrices/cross-product | dominant allocation unchanged |

BPCells support is rejected for this strict-equivalence release. No PR C is
created, and no dependency declaration changes.
