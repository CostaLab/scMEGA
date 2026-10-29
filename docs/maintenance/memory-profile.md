# scMEGA memory profile

Baseline commit: `3b9e904960982fc189d98bfceaaa313984f43a96`.
Profile date: 2026-10-28. Platform: macOS arm64, R 4.5.3, Seurat 5.5.0,
SeuratObject 5.4.0, Signac 1.17.1, Matrix 1.7-5, and scMEGA 1.1.0.

## Cardiomyocyte call graph

The tutorial calls scMEGA `CoembedData()`, then Seurat Harmony/UMAP/clustering,
scMEGA `PairCells()`, `CreatePairedObject()`, `AddTrajectory()`, `SelectTFs()`,
`SelectGenes()`, `GetTFGeneCorrelation()`, `GetGRN()`, and plotting helpers. The
observed peak occurs in the first call, `CoembedData()`:

```text
CoembedData
  -> NormalizeData / ScaleData (Seurat)
  -> FindTransferAnchors(reduction = "cca") (Seurat)
       -> ScaleData(reference and query)
       -> RunCCA
            -> Standardize two dense matrices
            -> crossprod(reference, query)
            -> irlba
  -> TransferData / merge / PCA / UMAP (Seurat)
```

The full input is 28,933 genes × 63,161 RNA cells, 182,116 peaks × 19,859
ATAC cells, and a 24,919 × 19,859 sparse gene-activity matrix. There are 19,026
common nonzero genes entering CCA. Source files and SHA-256 values are recorded in
`docs/maintenance/tutorial-validation.md` and `scripts/maintenance/tutorial-data.tsv`
on Draft PR #45; this branch does not modify them.

## Ranked bottlenecks

| Rank | Function / allocation | Input dimensions | Size before overhead | Evidence |
|---:|---|---:|---:|---|
| 1 | Seurat `RunCCA.default()` dense cross-product | 63,161 × 19,859 | 9.35 GiB | Source inspection: `crossprod(object1, object2)` |
| 2 | Standardized reference matrix | 19,026 × 63,161 | 8.95 GiB | Seurat `ScaleData` plus `Standardize` |
| 3 | Standardized query matrix | 19,026 × 19,859 | 2.82 GiB | Seurat `ScaleData` plus `Standardize` |
| 4 | IRLBA bases/workspace and CCA merged objects | cell matrices above | variable | Live alongside the preceding dense objects |
| 5 | Redundant scMEGA query scale layer | 19,026 × 19,859 | 2.82 GiB | scMEGA scaled before Seurat scaled the diet query again |

The matrix sizes exclude R headers, dimnames, allocator overhead, temporary
copies, the original Seurat objects, sparse assays, and IRLBA work arrays.

## Measurements

| Function | Input dimensions | Peak RSS | Allocated memory | Runtime | Bottleneck |
|---|---|---:|---:|---:|---|
| Baseline `CoembedData` fixture | 2,000 genes; 1,200 RNA; 600 ATAC; 182,116 peaks | 1,904,787,456 B | 550,075,480 B in allocations ≥1 MB | 4.391 s | CCA plus duplicate query scale |
| Optimized `CoembedData` fixture | identical | 1,891,352,576 B | 540,475,432 B in allocations ≥1 MB | 4.447 s | CCA remains dominant |
| Baseline full preprocessing only | 19,026 genes; 63,161 RNA; 19,859 ATAC | sampling dominated by input load; not used for reduction claim | retained 3,022,698,672-B scale matrix | 24.563 s | redundant pre-anchor scaling |
| Delayed full preprocessing only | identical | sampling dominated by input load; not used for reduction claim | no retained pre-anchor scale matrix | 11.022 s | later reconstructed for output equivalence |

`Rprofmem()` measures R allocations and `/usr/bin/time -l` measures process peak
RSS. The earlier full tutorial retry reported 46,314,422,272 B maximum resident
memory and 96,461,476,840 B peak footprint, with about 47 GB swap. Those values
came from the Phase VII `/usr/bin/time -l` and `sysctl vm.swapusage` logs; they
are historical baseline evidence and were not reproduced uncontrolled here.

## Conclusion

The dominant memory requirement belongs to Seurat's exact dense CCA. Sparse input
alone cannot avoid its mathematically dense cell-by-cell cross-product. Replacing
CCA with RPCA, approximate neighbors, fewer features/cells, or a custom implicit
SVD would change dependency internals or analytical behavior and was rejected for
this maintenance optimization.
