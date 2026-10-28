# Tutorial validation

Date: 2026-10-28. Integrated scMEGA commit:
`3b9e904960982fc189d98bfceaaa313984f43a96` on `main`.
Tutorial review: Draft PR [#45](https://github.com/CostaLab/scMEGA/pull/45).

| Tutorial | Execution | Figures | Scientific consistency | Status |
|---|---|---|---|---|
| Myofibroblast differentiation | Complete, 739 s | 20 embedded PNG and 1 SVG regenerated | Pair identities/count features preserved; GRN completed | **Passed** |
| Cardiomyocyte remodeling | Failed during full-data CCA | Not regenerated | No result comparison possible; workflow was not downsampled | **Blocked by memory** |
| CD4 T-cell activation | Complete, 362 s | 11 embedded PNG and 1 SVG regenerated | Paired RNA/ATAC cell names identical; trajectory and GRN completed | **Passed** |

No website was deployed and no release was created. Generated HTML and analysis
objects remain in isolated `/private/tmp/scmega-tutorial-runs` directories pending
review; large input data are not in Git.

## Integrated package gate

Fresh detached checkout: `3b9e904`. Platform: macOS arm64, R 4.5.3. Important
versions are scMEGA 1.1.0, Seurat 5.5.0, SeuratObject 5.4.0, Signac 1.17.1,
ArchR 1.0.3, Matrix 1.7-5, irlba 2.3.7, chromVAR 1.32.0, MOJITOO 1.0 at
`1ee6af174fc8ba31641f82baaec589ec6ebb4576`, and BSgenome hg38 1.4.5. The full
package inventory and loaded session are in `tutorial-environment.tsv` and
`tutorial-session-info.txt`.

Validation results:

- isolated `R CMD INSTALL`: passed; `library(scMEGA)` passed;
- testthat: 187 passed, 0 failed, 0 skipped, 0 warnings;
- maintenance security: 22 passed;
- `R CMD build .`: passed, producing `scMEGA_1.1.0.tar.gz`;
- ordinary `R CMD check`: 0 errors, 0 warnings, 2 notes;
- online `R CMD check --as-cran`: remote CRAN archive metadata failed TLS before
  package checking; the attempt is not reported as a package pass;
- `--as-cran` with only remote incoming checks disabled: 0 errors, 1 warning,
  5 notes. The warning is missing local `qpdf`; notes cover `Remotes`, clock
  verification, four private ArchR helpers, existing data-mask diagnostics, and
  old HTML Tidy;
- final hosted run
  [40696699280](https://github.com/CostaLab/scMEGA/actions/runs/40696699280): passed on
  `3b9e904`.

## Inputs and provenance

All ten HTTPS URLs returned HTTP 200 with certificate verification enabled. Exact
sizes, SHA-256 hashes, destinations, and URLs are recorded in
`scripts/maintenance/tutorial-data.tsv`. No source object was modified.

Input summaries:

- Fibroblast RNA: Seurat/legacy Assay, 28,933 features × 45,515 cells; ATAC:
  Seurat/ChromatinAssay, 173,843 peaks × 6,481 cells; gene activity:
  `dgCMatrix`, 24,919 × 6,481.
- Cardiomyocyte RNA: Seurat/legacy Assay, 28,933 features × 63,161 cells; ATAC:
  Seurat/ChromatinAssay, 182,116 peaks × 19,859 cells; gene activity:
  `dgCMatrix`, 24,919 × 19,859.
- PBMC 10x HDF5: 36,601 genes and 143,887 peaks across 11,898 cells. The label
  reference updates to 29,732 features × 11,909 cells.

The article downloads no longer use `--no-check-certificate`; downloads are
conditional and use verified HTTPS. The release-validation workflow verifies file
sizes and SHA-256 hashes before execution.

## Source compatibility corrections

Changes are limited to tutorial and validation infrastructure:

- preserve real `ChromatinAssay` objects instead of coercing them to `Assay5`;
- fix the cardiomyocyte `reduction.ke` typo to `reduction.key`;
- set a deterministic seed;
- replace removed Signac 1.17 `RunChromVAR()` calls with its Signac 1.16-equivalent
  public chromVAR sequence: motif matrix, nonzero peaks, GC bias, background peaks,
  deviation z-scores, and a Seurat data assay;
- pin MOJITOO for website/release validation.

The chromVAR change restores an upstream wrapper removed by Signac; it does not
change scMEGA algorithms or select a different motif/background method.

## Myofibroblast differentiation

Source: `vignettes/articles/myofibroblast-GRN.Rmd`.

The first attempt reached TF selection and failed because Signac 1.17 no longer
exports `RunChromVAR()`. After the public-API correction, a clean rerun completed
all 71 knit steps in 739 seconds. Maximum resident memory was 33.5 GB (reported
peak footprint 46.1 GB).

Verified outputs:

- 3,319 cells in the saved trajectory object;
- RNA and ATAC output column names are identical (SHA-256 of the ordered names:
  `1c8de0e1fdf3deabbde04fe55edfcb0eacbf61baedda92a09ca3c2a627ee8353`);
- 28,933 RNA features and 173,843 ATAC features retain stable ordered hashes;
- GRN: 50,595 edges, 131 TFs, and 2,050 target genes;
- saved assays: RNA, ATAC, chromvar, and target;
- regenerated standalone HTML: approximately 20 MB.

No old serialized paired/trajectory/GRN objects are published, so edge-by-edge and
embedding-coordinate comparison to the old website is impossible. The article's
biological interpretation and selected workflow are unchanged. Differences are
classified as expected compatibility/output regeneration, not evidence of a new
biological conclusion.

## Cardiomyocyte remodeling

Source: `vignettes/articles/cardiomyocyte-GRN.Rmd`.

The normal run stopped in `CoembedData()`/Seurat CCA at R's 64 GB vector limit
after 129 seconds. One controlled clean retry used `R_MAX_VSIZE=100Gb`; it reached
46.3 GB maximum resident memory and about 96.5 GB peak memory footprint, consumed
approximately 47 GB swap, and was manually interrupted to prevent host failure.

No downsampling, altered features, or changed integration method was used to force
completion. Figures and biological results were therefore not regenerated. The
current checked-in cardiomyocyte HTML also references 21 missing relative figure
files, an existing documentation defect that cannot be safely repaired without a
successful full render. A host with substantially more than 64 GB RAM (or a
maintainer-approved scientifically equivalent memory strategy) is required.

## CD4 T-cell activation

Source: `vignettes/articles/pbmc_10x_multiome.Rmd`.

A clean run completed all 59 knit steps in 362 seconds. Maximum resident memory was
25.0 GB (reported peak footprint 27.6 GB).

Verified outputs:

- 2,999 trajectory cells, all with nonmissing pseudotime spanning 0.0333–100;
- RNA and ATAC ordered cell names are identical (SHA-256:
  `aaeafc7b0eed23752f922d9db1f6d6fffad8209fe35c1baeb84aa659e781e325`);
- 71 selected TFs, 2,233 selected genes, and 21,932 peak-to-gene links;
- raw GRN: 83,128 rows; correlation-filtered visualization GRN: 34,926 edges,
  70 TFs, and 1,990 target genes;
- regenerated standalone HTML: approximately 10 MB.

The published repository contains no serialized PBMC result baseline, so exact
edge/ranking equivalence cannot be established. Cell correspondence, dimensions,
identifiers, full workflow completion, and the existing biological narrative are
preserved. Stochastic steps use a documented seed; numerical differences from old
unseeded figures are classified as stochastic unless a future baseline proves
otherwise.

## Website validation

pkgdown successfully rebuilt the home page, all 31 function-reference pages, news,
and search into `/private/tmp/scmega-site-preview`. The two newly rendered,
self-contained articles have zero missing relative links. A complete
`pkgdown::build_site()` cannot pass because it necessarily reruns the blocked
cardiomyocyte article. The copied historical cardiomyocyte page has 21 missing
figure links, while navigation and regenerated function reference files exist.
Nothing was deployed to GitHub Pages.

## Sustainable CI

`.github/workflows/tutorial-validation.yml` adds:

- a lightweight PR job that parses all three Rmd files, rejects disabled TLS,
  verifies referenced scMEGA exports, validates the manifest, and checks pkgdown
  navigation;
- a manual six-hour full-data job for one tutorial at a time, with checksum-keyed
  data caching, installed checked-out scMEGA, rendered artifact upload, read-only
  permissions, and no deployment.

Local static validation passes for all three sources. The full workflow remains
manual so routine PRs do not download 1–5 GB biological inputs.

## Conclusion

The integrated package and two tutorials are validated. The maintenance release is
**not ready for publication** because the cardiomyocyte tutorial cannot be
regenerated on the available 64 GB host and its existing website figures are
missing. Maintainer review should decide whether to provision a larger runner or
approve a scientifically validated memory-reduction strategy before merging the
tutorial PR or publishing a release.
