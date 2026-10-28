# Memory optimization benchmark

All comparisons use commit `3b9e904` as baseline, the same R process environment,
seed 1234, and the same deterministic fixture: the first 2,000 sorted common genes,
1,200 sorted RNA cells, 600 sorted ATAC cells, all 182,116 ATAC peaks, and the same
30 transfer dimensions. The fixture SHA-256 is
`e9b8d650fabe639f359555f0803c28a780f02e7ad4d8b40d4bc80ed8fbd8bf5e`.

| Operation | Original peak RSS | Optimized peak RSS | Reduction | Runtime change | Equivalent |
|---|---:|---:|---:|---:|---|
| Representative `CoembedData()` | 1,904,787,456 B | 1,891,352,576 B | 0.705% (13,434,880 B) | 4.391 → 4.447 s (+1.28%) | Exact |
| R allocations ≥1 MB | 550,075,480 B | 540,475,432 B | 1.745% (9,600,048 B) | n/a | Exact |
| Retained pre-anchor full query scale | 3,022,698,672 B | 0 B during CCA | 100% for this duplicate | reconstructed after transfer | Exact on fixture |
| Whole Cardiomyocyte workflow | 46.31 GB historical RSS; 96.46 GB footprint | guarded run stopped at 33.80 GiB RSS / 20.23 GiB swap | not comparable: neither run completed | not comparable | Not reached |

The whole-workflow row is not presented as a memory reduction because the optimized
run was intentionally stopped at a different safety threshold. It demonstrates
only that the optimization does not make the full CCA feasible on this 64 GB host.

## Exact result comparison

Every saved component is `identical()` between baseline and optimized runs:

- ordered cells and metadata;
- variable features and assay list;
- retained GeneActivity `scale.data` (maximum absolute difference 0);
- RNA data (maximum absolute difference 0);
- PCA and UMAP embeddings (maximum absolute difference 0).

Reference component SHA-256 digests:

| Component | SHA-256 |
|---|---|
| cells | `2e315e8a231eb92d1b6005e2b7d7c58554b9df06505acd41fa40523959c6236c` |
| metadata | `11b30d548ef83261344ce36c0f182a1fa9d4f0db919fa3bc4cafdce8c2916673` |
| variable features | `b78f35634599ad7b72d587403a4e8a949fda84d4847775baee2af850081677f2` |
| assays | `48e10d8d02d91d17affd1a354a1f77c22ccb32391a13fa5f51e4cf039c8c0454` |
| GeneActivity scale | `63d37d7cacb789985ff4c92b9939e4ffcf51443f7f6e372292b7ee6a0b032019` |
| RNA data | `bb6f9645dcf93ffec6dfd026ee4c76a7f492f8bbc5d000678c8bb208b7a1b499` |
| PCA | `0701baedf4a6773bd632040515eb859ad50d7d065f24388d7c909e8e97d15cf2` |
| UMAP | `bbf886bcbb2886634c1da948db83ce5e2249650d520c98c2c484288178e1904e` |

The benchmark was executed once per state because the result was exact and the
measured reduction equals the removed dense allocation. More repetitions are
recommended on a dedicated quiet runner before making broader performance claims.
