Thank you for the detailed warnings. The 2,976-versus-4,224 mismatch is explained by the split RNA layers in the report: `LayerData(..., layer="data")` selects only the first matching layer, while the anchors use all reference cells. I reproduced the same error with 60 selected columns versus 120 reference cells. This does not require switching ATAC preprocessing from Signac to ArchR.

If your `counts.*`/`data.*` layers are disjoint batches of the same RNA assay that should be used together, explicitly join them on a copy:

```r
SeuratObject::Layers(RNA[["RNA"]], search = NA)
rna.work <- SeuratObject::JoinLayers(RNA, assay = "RNA")
refdata <- SeuratObject::LayerData(rna.work, assay = "RNA", layer = "data")
stopifnot(setequal(colnames(refdata), SeuratObject::Cells(rna.work[["RNA"]])))
stopifnot(length(SeuratObject::VariableFeatures(rna.work[["RNA"]])) > 0)
```

Inspect layer membership first; do not join alternative measurements/normalizations of the same cells blindly. The RNA reference needs normalized data and your intended variable features. If either is absent, complete the appropriate RNA preprocessing; do not change normalization or feature selection merely to silence an error.

Your example also uses `weight.reduction="umap"` with the default `dims=1:30`. A two-dimensional UMAP does not contain 30 dimensions. Select an appropriate ATAC reduction and dimensions for your analysis. For example, this call shape was tested with a synthetic LSI reduction:

```r
use.dims <- 1:10  # choose suitable LSI components for your own data
stopifnot(max(use.dims) <= ncol(SeuratObject::Embeddings(ATAC, "lsi")))
stopifnot(setequal(colnames(gene.activities), colnames(ATAC)))
gene.activities <- gene.activities[, colnames(ATAC), drop = FALSE]
obj.coembed <- CoembedData(
  obj.rna = rna.work, obj.atac = ATAC, gene.activity = gene.activities,
  reference.assay = "RNA", weight.reduction = "lsi", dims = use.dims
)
```

Verification: joining the intended RNA layers plus normal RNA preprocessing completed real Seurat anchor transfer/PCA/UMAP on 220 synthetic cells. A separate legacy-Assay failure was also reproduced: `JoinLayers()` has no legacy Assay method, and legacy merging drops the reference variable-feature set. A local correction now checks the assay class, restores the existing reference feature selection only when merging lost it, and reports ambiguous RNA layers before anchor calculation. It does not choose a new feature set or automatically combine ambiguous layers.

The corrected function passes 17 expectations, including both assay classes, preservation of observed RNA data/ATAC counts and cell IDs, retained reference variable features, and finite embeddings. Five upstream warnings are recorded. Changes are local on `codex/issue-27-coembed`; no PR is published, and installed-package tests also pass with the separately reviewed metadata proposal. R CMD check still reports six warnings and three notes. The synthetic test verifies this mechanism, not your original plant dataset or every Seurat development version.
