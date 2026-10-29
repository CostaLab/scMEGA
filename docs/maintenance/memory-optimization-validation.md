# Memory optimization validation

## A. Baseline

- Commit: `3b9e904960982fc189d98bfceaaa313984f43a96`.
- Full input: 63,161 RNA cells, 19,859 ATAC cells, 28,933 RNA features,
  182,116 peaks, and a 24,919 × 19,859 sparse gene-activity matrix.
- Input files:

| File | Bytes | SHA-256 |
|---|---:|---|
| `snRNA.rds` | 1,291,963,378 | `234edb0b3544dcb97084d1de684712b4c9b77f665f2e3ba1fe4ea17a3dd282c2` |
| `snATAC.rds` | 1,782,277,344 | `06603ef6b93a3e9e67352fe0dd65b5891733dbc0960ff4a8740241a8b42c151a` |
| `gene.activity.rds` | 1,251,316,755 | `d2771ddc2188074e3e1d24a016720032710a7cbbd5c06ce1dae6976edb82eaaf` |

All are from `https://costalab.ukaachen.de/open_data/scMEGA/Cardiomyocyte/` and
were reverified after profiling.
- Historical failed full run: 384.02 seconds, 46.31 GB maximum resident memory,
  96.46 GB peak footprint, and approximately 47 GB swap.
- Deterministic correctness fixture and output digests are documented in
  `memory-benchmark.md`; no large biological fixture is committed.

## B. Identified bottlenecks

1. Seurat exact CCA creates a dense 63,161 × 19,859 cross-product.
2. Reference/query scaling and `Standardize()` hold multiple 19,026-feature dense
   matrices at the same time.
3. scMEGA previously retained a query scale layer before Seurat created another
   scaled diet query internally.
4. IRLBA and merged CCA objects add live workspace and copies.

The top two bottlenecks are dependency-algorithm behavior, not an accidental sparse
conversion in scMEGA.

## C. Implemented optimization

File: `R/pair_cells.R`, function `CoembedData()`.

The GeneActivity assay is still normalized before anchoring. Its retained
`scale.data` layer is now calculated after anchor transfer rather than before
`FindTransferAnchors()`. Seurat continues to scale the diet reference/query for
CCA exactly as before, and scMEGA reconstructs the same returned GeneActivity
scale layer. No signature, default, feature, cell, threshold, normalization, CCA,
or downstream operation changes.

## D. Result equivalence

The representative real-data output is exactly identical for cell/feature order,
metadata, assays, variable features, GeneActivity scale values, RNA values, PCA,
and UMAP. Maximum numerical difference is zero. Existing valid joined/scaled/legacy
Seurat paths pass 34 focused expectations, including exact retained scale data.
The complete installed-source suite passes 190 expectations with 0 failures,
0 skips, and 0 warnings; all 22 maintenance security tests pass. `R CMD build`
passes and ordinary `R CMD check` reports 0 errors, 0 warnings, and the same two
pre-existing notes (private ArchR helpers and data-mask/global diagnostics).

Because the optimization acts only before coembedding output exists, downstream
pairing, trajectory, TF/gene ranking, peak-to-gene links, GRN edges, and graph
statistics receive identical inputs. No scientific threshold decision changes.

## E. Full tutorial

One guarded optimized full-data attempt used the original data and arguments with
no subsampling. A monitor terminated it during `FindTransferAnchors()` when swap
reached 20.23 GiB; RSS was 33.80 GiB. The process had printed “Performing data
integration using Seurat...” and had not returned anchors. The 50 GB RSS threshold
was not reached. No outputs were produced.

The full tutorial therefore still cannot complete safely on this 64 GB host. The
remaining blocker is Seurat's exact dense CCA, not scMEGA's redundant scale layer.

## F. Performance

Representative peak RSS improves by 13.43 MB (0.705%) and large R allocations by
9.60 MB (1.745%). On the full input, the delayed duplicate is 3.02 GB, but a safe
whole-workflow percentage cannot be claimed because the workflow remains blocked.
Runtime on the representative fixture changes by +1.28%, within ordinary noise.

## G. Remaining limitations

- Exact CCA remains intrinsically memory intensive at these cell counts.
- An implicit/blockwise CCA SVD could avoid materializing the cross-product, but
  integrating it would require reimplementing or changing Seurat internals and may
  alter convergence/numerics. It is a separate research task, not accepted here.
- RPCA, subsampling, fewer features, approximate algorithms, lower precision, and
  threshold changes were not used.
- Disk-backed inputs would not prevent Seurat's current dense CCA materialization.
- A high-memory runner is still required for the complete tutorial.

## H. Pull request

Draft PR [#46](https://github.com/CostaLab/scMEGA/pull/46) proposes the isolated
`codex/memory-optimization` change. It should remain a Draft until hosted package
checks and maintainer review confirm the modest benefit justifies the sequencing
change. No release, Issue, or other tutorial was modified.
