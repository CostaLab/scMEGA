# Seurat v5 compatibility

scMEGA requires Seurat >= 5.0.0 and SeuratObject >= 5.0.0. These requirements
are enforced by `DESCRIPTION`; redundant per-function version checks are not used.

## Supported assays

- `Assay5` with one unambiguous requested `counts`, `data`, or `scale.data` layer;
- legacy `Assay` objects that remain supported by SeuratObject v5;
- Signac `ChromatinAssay` objects, retaining ranges, fragments, and motif metadata.

Functions that require a unified matrix now inspect matching layers explicitly.
When multiple `counts.*` or `data.*` layers are present, scMEGA reports the layer
names and asks the caller to join only the biologically intended disjoint layers
with `JoinLayers()`. It does not silently choose the first match or automatically
combine ambiguous measurements. Legacy Assay and ChromatinAssay objects are not
passed blindly to `JoinLayers()`.

The public SeuratObject APIs `Assays()`, `Layers()`, `LayerData()`, `Embeddings()`,
`Reductions()`, and `object[[]]` are used for data, reductions, and metadata. Public
`Signac::GetMotifData()` replaces direct motif-slot access.

## BPCells status

BPCells 0.3.1 at commit `5f2a5d21a21f6e75cf3650dee8b65f6220f50a41`
was evaluated but is not declared as supported by this release:

- integer RNA counts, names, and subsetting remain disk-backed and exact;
- Seurat normalization remains a BPCells queued matrix, but differs from the
  current dgCMatrix path by up to `4.82422e-7` on the test fixture;
- the difference reaches `1.838324e-5` after trajectory aggregation;
- legacy `CreateAssayObject()` and Signac `CreateChromatinAssay()` materialize the
  test matrix as `dgCMatrix`;
- exact CCA still constructs dense scaled matrices and a dense cell cross-product.

Because strict scientific equivalence is not established and the Cardiomyocyte
peak occurs after dense CCA materialization, BPCells is not added to `Suggests` or
used automatically. `scripts/maintenance/evaluate_bpcells.R` reproduces the audit
when BPCells is installed separately.

## Sketching

Seurat v5 sketch workflows could reduce resource needs by analyzing a subset of
cells, but that changes the statistical analysis. Sketching is not enabled or used
for validation. It remains a possible future optional mode requiring its own
scientific contract.
