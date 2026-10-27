# Phase II validation and review readiness

Dates: 2026-10-27–28. Reporting base dde1d7c; upstream main e6ce4f9. All 12 open issues were processed sequentially and refreshed again at the end of issue investigation. No remote write or branch merge occurred.

## Environment

R 4.5.3; Seurat 5.5.0; SeuratObject 5.4.0; Signac 1.17.1; Matrix 1.7-5; primary irlba 2.3.7; testthat 3.3.2; ArchR 1.0.3. The alternative irlba 2.3.5.1 test build is isolated separately. All original 12 missing direct dependencies are now available. The exact 342-package dependency closure, metadata/shared-library checksums, source archive hashes and GitHub SHAs, system-library package URLs and compiler settings are archived in phase2-evidence/.

Targeted repairs addressed observed failures: missing CMake, archived smoother/TTR, missing GSL/libgit2, duplicate OpenMP linkage (isolated pthreads BLAS selected), and C++14 required by installed RcppArmadillo headers. No duplicate-OpenMP bypass, algorithm edits or existing site-library upgrades were used. The library is reusable locally and checked for drift; portable restoration, a coherent release-only Bioconductor matrix and hosted Linux CI remain unverified.

## Actual commands/results

Commands were run from the reporting repository with `R_LIBS_USER=/tmp/scmega-r-lib` unless another prefix is shown.

| Command / target | Actual result |
|---|---|
| `gh issue list --repo CostaLab/scMEGA --state open --limit 100 --json number,title,body,labels,url`, then full `gh issue view` for every issue | 12 issues plus all returned comments; #21 screenshot inspected, #29 immutable screenshot investigation preserved; final numbers/update timestamps unchanged |
| `Rscript scripts/maintenance/run_pairing_tests.R .` with `SCMEGA_REQUIRE_INTEGRATION=true` | 117 expectations pass; 2 warnings; 0 failures/errors/skips |
| `Rscript docs/maintenance/examples/run-issue25-tests.R /tmp/scmega-issue-25` | 29 expectations pass; 0 warnings/failures/errors/skips; original source had 4 test errors |
| `Rscript docs/maintenance/examples/run-issue27-tests.R /tmp/scmega-issue-27` | 17 expectations pass; 5 upstream warnings; 0 failures/errors/skips; original source had 2 test errors |
| `Rscript docs/maintenance/examples/issue-31-matrix-irlba.R` | Current irlba sparse SVD agrees with dense reference; Signac TF-IDF/SVD pipeline passes |
| `R_LIBS_USER=/tmp/scmega-irlba2351-lib:/tmp/scmega-r-lib Rscript docs/maintenance/examples/issue-31-matrix-irlba.R` | Same tests pass after rebuilding archived irlba 2.3.5.1 against installed Matrix; original Ubuntu/R 4.3.1 binary unavailable |
| `Rscript docs/maintenance/examples/issue-24-seurat-layers.R` | Assay/Assay5 values, split-layer behavior/join, absent data and ChromatinAssay range preservation verified |
| `Rscript docs/maintenance/examples/issue-27-coembed.R .` | Diagnostic original-source run reproduces split-data mismatch; properly preprocessed/joined Assay5 succeeds on 220 cells; original legacy Assay fails at JoinLayers. Caught errors are diagnostic outcomes, not passing workflows |
| `Rscript docs/maintenance/examples/issue-33-motif-identifiers.R` | Real Signac motif metadata, exact-label loss, duplicate/composite ambiguity and mm10 annotation contract verified; duplicate-name warning recorded. Initial invalid fixture lacked PWM metadata and was corrected before making conclusions |
| `Rscript docs/maintenance/examples/issue-23-missing-layers.R` | Missing RNA data failure reproduced; appropriate normalization restores trajectory; ordinary data-only assay holds signed TF activity; 2 relevant layer warnings recorded |
| `Rscript docs/maintenance/examples/issue-21-network-labels.R` | Disconnected graph fails as expected; selected weak component renders with RUNX1 label; original graph unchanged; 2 guide warnings recorded |
| `Rscript docs/maintenance/examples/issue-17-genome-contract.R` | Unsupported rn6 rejected; 3 mm10 links match actual C++/direct named-row correlations after gene permutation |
| `Rscript docs/maintenance/examples/issue-30-grn-meaning.R` | GetGRN joins/counts supplied support; missing/negative correlations remain; no trajectory-free inference implied |
| `python3 -m unittest discover -s tests/maintenance -v` | 22 passed; no infrastructure changes or weakened checks |
| R parsing of all example/reply R blocks; Python compileall | Passed; parsing is not claimed as execution of user-data-dependent examples |
| `Rscript docs/maintenance/examples/check-environment.R` | R and all 342 recorded versions/DESCRIPTION/shared-library fingerprints match |

## Installed-package checks

`docs/maintenance/examples/check-package.py` runs actual `R CMD INSTALL`, `R CMD build --no-build-vignettes --no-manual`, and built-tarball `R CMD check --no-manual --no-build-vignettes`, recording each command/exit code. Source-test switches are removed from its environment. Full biological vignette rendering was not attempted.

Three scratch copies were produced with `prepare-validation-copy.py`: each independent issue branch plus metadata commit ec24fdf. This is a filesystem-only validation composition, not a Git merge, and does not change any branch. A fourth scratch copy contains metadata alone.

| Exact candidate | Install | Build | Built-tarball check | Installed testthat output |
|---|---|---|---|---|
| Reporting base / existing #32, original metadata | Exit 0 | Exit 0 | Exit 1: undeclared Seurat | Dependency-stage check stops before tests |
| Metadata ec24fdf alone | Exit 0 | Exit 0 | Exit 0; **6 warnings, 3 notes**, no errors | No baseline test suite on main |
| #25 98ff9ad + metadata | Exit 0 | Exit 0 | Exit 0; **6 warnings, 3 notes**, no errors | FAIL 0 / WARN 0 / SKIP 0 / PASS 29 |
| #27 7fbcce6 + metadata | Exit 0 | Exit 0 | Exit 0; **6 warnings, 3 notes**, no errors | FAIL 0 / WARN 5 / SKIP 0 / PASS 17 |
| Existing #32 branch + metadata | Exit 0 | Exit 0 | Exit 0; **6 warnings, 3 notes**, no errors | FAIL 0 / WARN 2 / SKIP 0 / PASS 117 |

The five `final-*-results.json` files contain exact paths/commands and exit codes; corresponding logs and installed test output are archived. Earlier failed installation attempts are retained as historical evidence and are superseded by these results.

Warnings already present on the metadata-only baseline concern missing object/data documentation, stale/missing Rd argument documentation, LazyData compression, serialized data compression, and vignette registration/build output. Notes concern private ArchR APIs, code-analysis globals and vignette recognition. These have not been hidden or converted into a clean result. Code-analysis findings need contextual review; do not blindly add imports/global declarations to silence data-mask symbols or real undefined variables.

## Readiness

#25/#27/#32 are ready for human review as locally verified corrections. Installed-package verification depends on the separate metadata proposal; each issue branch stays independent. Local PR descriptions are in pr-drafts/. Remaining package warnings still fail the existing CI warning threshold, and branch integration needs review/testing. Metadata is reviewable with successful load/unload/check execution, not a claim that every existing function or supported R release is validated.

Replies #21/#29/#30 are supported by source/examples; #31 has a tested rebuild mechanism; #24/#23 need original/current reproductions for blanket resolution; #33/#18/#17 need scientific mappings or annotation/interface decisions and data. A reply draft does not resolve or close an issue. No comments/PRs were published. Hosted authorization, sandbox/prompt-injection, Codex generation, artifact transfer and real draft PR creation remain untested. Automatic publication remains disabled.
