# Final tutorial and release validation

Date: 2026-10-29. Version proposed for publication: **1.1.0**. No tag, release,
Issue closure, or website deployment has been performed.

## Tutorials

| Tutorial | Full execution | Data mode | Scientific validation | Status |
|---|---|---|---|---|
| Myofibroblast differentiation | 739 s; 33.5 GB max RSS | Full | 3,319 paired trajectory cells; 50,595-edge GRN; identifiers aligned | Passed |
| Cardiomyocyte remodeling | 1,605.65 s; 38.76 GiB max RSS | Fixed downsample | all strata/groups retained; 5,450 paired trajectory cells; 26,440-edge GRN | Passed |
| CD4 T-cell activation | 362 s; 25.0 GB max RSS | Full | 2,999 trajectory cells; 34,926-edge filtered GRN; paired identifiers aligned | Passed |

All runs used clean R sessions, no knitr cache, verified input checksums, and the
local combined #47 + #46 package candidate. Generated inputs/results remain outside
Git. Myofibroblast and CD4 evidence was reused because the release-candidate source
changes are exact-equivalent at tested package boundaries and neither tutorial
source changed after its complete run.

## Cardiomyocyte sampling

- Original: 63,161 RNA cells and 19,859 ATAC cells.
- Published fixed selection: 40,000 RNA and 8,000 ATAC cells.
- Strata: `patient_group × patient_region_id`.
- Seeds: 42 for RNA and 43 for ATAC.
- Allocation: proportional largest-remainder; original within-modality order retained.
- All 20 RNA and 17 ATAC patient-region strata retained.
- Selected groups: RNA 34,640 myogenic / 2,688 ischemic / 2,672 fibrotic;
  ATAC 5,842 myogenic / 303 ischemic / 1,855 fibrotic.
- GeneActivity columns are exactly identical to the ordered selected ATAC vector.
- RNA cell record SHA-256:
  `1bf32443cddbd95c84a21680c06b64308d690cd6a2a08164f0d4357742a5a5d9`.
- ATAC cell record SHA-256:
  `43d962d39884b00c7a1668823c565e2dbf466e23f2e402d584bf1afc86ecc39a`.
- Runtime: 1,605.65 seconds.
- Peak RSS: 40,647,616 KiB (38.76 GiB).
- Peak system swap: 15,527 MiB (15.16 GiB), approximately 0.49 GiB above the
  pre-run level; the 50 GiB RSS / 20 GiB swap guard did not trigger.

The source contains explicit fixed counts and does not recalibrate on render.
Sampling is applied immediately after load. RNA and ATAC are sampled independently
because the source modalities are unpaired. Source counts and feature values are
not modified.

### Scientific outputs

- composition QC retained 34,202 coembedded cells;
- 6,056 ATAC cells were paired to candidate RNA cells;
- 5,450 cells remained on trajectory clusters 2 → 1 → 0;
- trajectory range: 0.01835–100;
- trajectory groups: 4,632 myogenic / 360 ischemic / 458 fibrotic;
- RNA and ATAC paired column names are identical;
- retained features: 28,933 genes and 182,116 peaks;
- GRN: 26,440 edges, 116 TFs, and 1,709 target genes;
- high-weight predictions include SP3/SERPINE1, IRF9/SERPINE1, SP1/RUNX1,
  and SP1/FN1;
- figures: 20 PNG and one SVG embedded in a 17 MB article.

Downsampling changes cell membership, clustering, trajectory composition, and GRN
weights. The result is scientifically valid but is not claimed to reproduce the
unavailable full-data result numerically. Patient/condition representation, RNA/
ATAC/GeneActivity alignment, complete feature sets, and the original analytical
method are preserved.

## Package candidate

- Default-branch baseline: `3b9e904960982fc189d98bfceaaa313984f43a96`.
- Locally integrated order: #47 → #46 → #45.
- Validated code/tutorial candidate before this report-only commit: `2ddc574`.
- R 4.5.3; Seurat 5.5.0; SeuratObject 5.4.0; Signac 1.17.1;
  ArchR 1.0.3; Matrix 1.7-5; chromVAR 1.32.0.
- Package installation and `library(scMEGA)`: passed.
- Regression suite: 201 passed, 0 failed, 0 skipped, 0 warnings.
- Maintenance security: 22 passed.
- Tutorial static/API/sampling contract: passed.
- `R CMD build`: passed (`scMEGA_1.1.0.tar.gz`).
- Ordinary `R CMD check`: 0 errors, 0 warnings, 2 known notes.
- Notes: four private ArchR helpers and existing data-mask/global diagnostics.

PR #46 and #47 both have green hosted package/security checks at their current
heads. PR #45 must obtain new hosted checks for the final tutorial head before
merge. A final post-merge check on `main` is required before tagging.

## Documentation website

All three articles rendered from source without cache. A local release preview
contains pkgdown home, news, search, all 31 reference pages, installation, and the
three regenerated self-contained articles. Automated checking found zero missing
relative links across the preview. Nothing was deployed.

## Release plan and current state

Already merged maintenance PRs: #35, #36, #37, #38, #41, #42, #43, and #44.

Pending approval:

1. [#47](https://github.com/CostaLab/scMEGA/pull/47) — Seurat v5 API/layer modernization.
2. [#46](https://github.com/CostaLab/scMEGA/pull/46) — exact-equivalent memory lifetime improvement.
3. [#45](https://github.com/CostaLab/scMEGA/pull/45) — tutorials, fixed Cardiomyocyte sampling, validation CI, and 1.1.0 notes.

After approval: merge sequentially, rebuild/check updated `main`, verify hosted CI,
tag `1.1.0`, publish GitHub Release `scMEGA 1.1.0`, verify installation from the
tag, then deploy the reviewed site through the repository's normal Pages process.

## Issue recommendations after release

- #32: close as resolved by the pairing correction, with recomputation guidance.
- #25: close for the reproduced trajectory-dimension defects; invite a new report
  if the unavailable original object still fails.
- #27: close for split/legacy layer handling, noting original plant data were not
  available.
- #29: eligible for closure as resolved usage/documentation guidance.
- #31 and #21: close only if maintainers consider the verified guidance sufficient.
- #23 and #24: leave open awaiting reproducible inputs.
- #17, #18, #30, and #33: leave open as feature/scientific-policy work.

## Known limitations

- ArchR remains a GitHub-only strong dependency and four private helpers are used.
- The strict CRAN incoming warning is not suppressed; this is a GitHub maintenance
  release, not a CRAN submission.
- Full-data Cardiomyocyte exact CCA still exceeds the available machine's safe
  memory budget; the published tutorial is explicitly downsampled.
- BPCells and sketching are not enabled.
- Custom genomes and species-aware/composite motif mapping remain unresolved.
