Thank you for reporting this. Several distinct issues are involved, so replacing one accessor is not enough to establish compatibility for every workflow.

Current scMEGA source already uses `LayerData()` in `GetTrajectory()`. Local tests verify counts/data access for legacy `Assay` and `Assay5` objects under SeuratObject 5.4.0, and layer access for Signac 1.17.1 `ChromatinAssay` without losing genomic ranges. This is not a claim that all scMEGA workflows or historical v5 development builds pass.

Please inspect the actual input and installed function:

```r
packageVersion("scMEGA")
packageVersion("Seurat")
packageVersion("SeuratObject")
packageVersion("Signac")
body(scMEGA::GetTrajectory)
assay_name <- "RNA"
class(object[[assay_name]])
SeuratObject::Layers(object[[assay_name]], search = NA)
```

Relevant cases:

- `SelectTFs()` asks for normalized `data` from RNA and TF-activity assays. A counts-only RNA assay needs the appropriate normalization; counts must not simply be relabeled as data. chromVAR activity is a separate upstream result.
- Split Assay5 layers require deliberate handling: selecting `layer="data"` may return only the first matching batch. If those layers represent disjoint cells from the same assay and are intended to be used together, explicitly join them and verify cell identities. See [Seurat's layer documentation](https://satijalab.org/seurat/articles/seurat5_essential_commands.html).
- Do not cast a `ChromatinAssay` to `Assay5` just to bypass an error. `SelectTFs()` and `GetTFGeneCorrelation()` need its motif metadata.
- A local #25 correction handles singleton trajectory bins/features; the #32 correction handles paired-count ordering. These are separate defects and are not yet published.
- Missing scMEGA dependency/import declarations, including `LayerData`/`JoinLayers` imports, have a separate local metadata proposal awaiting review.

Please provide your exact error and traceback, `sessionInfo()`, the relevant assay classes/layers, and a small synthetic or anonymized example. That will distinguish missing layers, matrix dimensions, split data and an import problem. No universal Seurat downgrade or wholesale assay conversion is recommended. The final isolated environment now supports package installation. Metadata-corrected scratch copies pass installed tests, but R CMD check retains six warnings and three notes; this does not establish compatibility for your original input.
