# Seurat v5 whole-pipeline memory profile

Date: 2026-10-29. Baseline commit:
`3b9e904960982fc189d98bfceaaa313984f43a96`. The modernization branch is based
directly on this commit and is independent of Draft PR #46.

## Actual Cardiomyocyte workflow

```text
load RNA, ATAC, gene activity
  -> CoembedData
       -> NormalizeData / ScaleData
       -> FindTransferAnchors(reduction = "cca")
       -> TransferData / merge / PCA / UMAP
  -> Harmony / clustering / marker analysis
  -> PairCells
  -> CreatePairedObject
  -> AddTrajectory
  -> motif matching / chromVAR
  -> SelectTFs / SelectGenes / GetTFGeneCorrelation
  -> GetGRN / visualization
```

The full workflow stops in `FindTransferAnchors()` before later stages exist.
Those downstream full-data stages therefore cannot be profiled honestly on this
host. Successful Myofibroblast and PBMC runs show that they are executable on
smaller inputs, but are not substitutes for Cardiomyocyte measurements.

## Stage evidence

| Stage | Dimensions / objects | Peak RSS / allocation evidence | Runtime | Result |
|---|---|---|---:|---|
| Data inputs | RNA 28,933 × 63,161; ATAC 182,116 × 19,859; activity 24,919 × 19,859 | Three source files total 4.33 GB compressed; in-memory objects are larger | not isolated | Loaded successfully |
| Pre-anchor preparation | 19,026 common genes | A retained query scale layer is 3,022,698,672 B | baseline 24.563 s in #46 stage audit | Optimized separately in #46 |
| Exact CCA | 63,161 RNA × 19,859 ATAC cells | 8.95 GiB reference scale, 2.82 GiB query scale, 9.35 GiB cross-product, plus Standardize/IRLBA/object copies | historical attempt stopped at 384.02 s | Dominant blocker |
| Guarded full attempt (#46) | unchanged full inputs | stopped at 33.80 GiB RSS / 20.23 GiB swap | stopped before anchors returned | Safe stop |
| Pairing through GRN | unavailable because CCA did not complete | not measured on full Cardiomyocyte data | n/a | Not reached |

The historical uncontrolled-boundary evidence was 46.31 GB maximum RSS,
96.46 GB peak footprint, and approximately 47 GB swap. It was not repeated.

## Allocation analysis

Seurat `RunCCA.default()` standardizes two dense matrices and executes
`crossprod(reference, query)` before IRLBA. BPCells and sparse input storage cannot
avoid this mathematically dense intermediate in the current Seurat implementation.
The stopped future was a `SequentialFuture`; worker duplication was not the cause
of the observed peak.

This PR changes access validation and public APIs only. It produces exact output
on the 1,800-cell real-data fixture, but peak RSS is 1,907,638,272 B versus
1,904,787,456 B at baseline, a 0.15% increase within measurement noise. No memory
improvement is claimed for accessor modernization.

The only accepted measured memory reduction remains Draft PR #46: 0.705% on the
representative fixture and removal of one 3.02 GB simultaneously retained full-data
scale layer. Exact CCA remains infeasible on the 64 GB host.
