# Tutorial validation

Updated: 2026-10-29. Default-branch baseline:
`3b9e904960982fc189d98bfceaaa313984f43a96`. Tutorial execution used the local
combined #47 + #46 candidate `c81406f` plus this Draft's tutorial source.
Tutorial review: Draft PR [#45](https://github.com/CostaLab/scMEGA/pull/45).

| Tutorial | Execution | Figures | Scientific consistency | Status |
|---|---|---|---|---|
| Myofibroblast differentiation | Complete, 739 s | 20 embedded PNG and 1 SVG regenerated | Pair identities/count features preserved; GRN completed | **Passed** |
| Cardiomyocyte remodeling | Complete, 1,606 s | 20 embedded PNG and 1 SVG regenerated | All strata retained; pairing, trajectory, and GRN completed | **Passed (downsampled)** |
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
- testthat: 198 passed, 0 failed, 0 skipped, 0 warnings on the combined candidate;
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
- apply fixed patient-group/patient-region stratified Cardiomyocyte sampling before
  integration, with exact selected-cell records and GeneActivity alignment;
- update Cardiomyocyte pairing to the current `CreatePairedObject()` interface,
  use public Signac motif access, and correct inherited fibroblast prose;
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

The source objects contain 63,161 RNA and 19,859 ATAC cells and cannot complete
exact CCA safely on the 64 GB host. The published example now uses a fixed sample
of 40,000 RNA and 8,000 ATAC cells. Sampling is proportional within
`patient_group × patient_region_id`, uses seeds 42/43, and preserves original
within-modality order. All 20 RNA and 17 ATAC strata and all three patient groups
are retained. GeneActivity columns are exactly the selected ATAC vector.

Selection records:

- RNA text SHA-256: `1bf32443cddbd95c84a21680c06b64308d690cd6a2a08164f0d4357742a5a5d9`;
- ATAC text SHA-256: `43d962d39884b00c7a1668823c565e2dbf466e23f2e402d584bf1afc86ecc39a`.

The final clean run completed all 69 knit steps in 1,605.65 seconds. Peak RSS was
38.76 GiB. Peak system swap was 15.16 GiB, only about 0.49 GiB above the pre-run
level; the 50 GiB RSS / 20 GiB swap guard did not trigger.

Verified outputs:

- 34,202 cells remain after explicit patient/modality composition QC;
- 5,450 paired trajectory cells with identical RNA/ATAC output identifiers;
- trajectory range 0.0183–100 across sampled clusters 2 → 1 → 0;
- trajectory cells retain myogenic, ischemic, and fibrotic groups;
- 28,933 RNA features and 182,116 ATAC peaks retain their identifiers;
- filtered GRN: 26,440 edges, 116 TFs, and 1,709 target genes;
- highly weighted predictions include SP3/SERPINE1, IRF9/SERPINE1,
  SP1/RUNX1, and SP1/FN1;
- regenerated standalone HTML: approximately 17 MB.

Downsampling necessarily changes clustering, selected cells, edge membership, and
weights. These results are not claimed numerically equivalent to the unavailable
full-data baseline. They preserve the complete method and major cardiomyocyte
remodeling/fibrotic genes while making the example executable.

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

All three articles now render successfully from source. pkgdown home, all 31
function-reference pages, news, and search were rebuilt locally; final combined
site/link validation is recorded in `final-tutorial-release-validation.md`.
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

The integrated package and all three tutorials are validated. The release remains
unpublished pending maintainer approval of PRs #47, #46, and #45, the downsampled
Cardiomyocyte scientific impact, version 1.1.0, and the final merge/tag plan.
