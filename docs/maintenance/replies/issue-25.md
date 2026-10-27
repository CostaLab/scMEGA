Thank you for the report. I reproduced this error in `GetTrajectory()` with synthetic Seurat objects. The expression `data.use[, cell_names]` can become a vector when a bin contains one cell **or when there is only one feature**, even if `cell_names` contains multiple cells. A different bin from the one inspected during debugging can also be the failing bin.

A local correction now retains matrix dimensions, keeps single-bin aggregation two-dimensional (tested with smoothing disabled), and preserves dimensions when `GetCorrelation()` has one shared TF. Empty bins now report an explicit error; their means cannot be defined without choosing a different analysis/binning policy. No counts, normalization formula, correlation method, or defaults were changed.

To inspect the input without changing it:

```r
assay_name <- "RNA"  # also inspect the assay used for chromVAR activity
SeuratObject::Layers(object[[assay_name]])
m <- SeuratObject::LayerData(object, assay = assay_name, layer = "data")
dim(m)
t <- object[["Trajectory"]][, 1]
groupEvery <- 1
breaks <- seq(0, 100, by = groupEvery)
table(cut(t, breaks = breaks, right = TRUE), useNA = "ifany")
setdiff(rownames(object[[]])[!is.na(t)], colnames(m))
```

Use the exact trajectory column and assay names in your object. Do not replace absent cells/TFs or empty bins with zeros. If you choose coarser bins, justify that resolution biologically. `SelectTFs()` currently fixes the smoothing window at 7 bins; coarsening to fewer than 7 bins is not a general workaround. Direct diagnostic `GetTrajectory(..., smoothWindow = NULL)` can help isolate aggregation from smoothing.

Verification: 29 source-level regression expectations pass for legacy Assay and Assay5 inputs, singleton cells/features, a single bin, normal means/normalization, invalid metadata, absent TFs, and one shared TF. The original source fails the new edge-case tests. The fix is local on `codex/issue-25-trajectory`; no PR has been published. With the separate metadata proposal applied in a scratch copy, installation and the installed-package tests also pass. R CMD check still reports six warnings and three notes; this is not a clean CI result.

Your original object was not available, so this establishes a reproducible defect with the same error, not proof of which dimension or bin failed in your data. If the problem persists after the reviewed correction, please provide `sessionInfo()`, `traceback()`, `dim(m)`, assay classes/layers, the occupancy table, and a small anonymized example. The dimension error itself is independent of species; other-species motif and genome support are separate questions.
