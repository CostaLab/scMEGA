# Validation report

Date: 2026-10-27. Baseline `e6ce4f9`. Local infrastructure commit `a9bea1e` precedes
pilot source edits. No remote writes, pushes, labels, settings or PRs occurred.
The implementation is locally reviewable; **production rollout is not complete**.

## Commands and actual outcomes

| Command | Outcome |
|---|---|
| `git status --short`; `git symbolic-ref refs/remotes/origin/HEAD` | Initial tree clean; main default, later confirmed by GitHub API. |
| `gh auth status`; `gh repo view CostaLab/scMEGA --json defaultBranchRef,nameWithOwner,viewerPermission` | Network-authorized retry succeeded; ADMIN viewer, main. Initial sandbox attempt failed. |
| `gh issue list --repo CostaLab/scMEGA --state open --limit 1000 --json number,title,body,labels,createdAt,updatedAt,url` plus `gh issue view N --json number,title,body,comments,labels,state,url,createdAt,updatedAt` for every issue | 12 open issues and full comments retrieved and classified before source edits. |
| `python3 scripts/maintenance/fetch_issues.py --repo CostaLab/scMEGA --output /tmp/scmega-live-audit.json` | Live REST pagination succeeded; all titles, bodies and comments match the 12-issue audit snapshot. Read-only. |
| `R CMD INSTALL --library=/tmp/scmega-baseline/lib .` | Failed: 13 missing Imports. |
| `R CMD build --no-build-vignettes --no-manual .` (baseline) | Passed. |
| `R CMD check --no-manual --no-build-vignettes --output=/tmp/scmega-baseline scMEGA_1.1.0.tar.gz` | Failed: 1 ERROR, missing Imports/testthat and undeclared Seurat. Network retry removed index-access warnings but retained the same dependency error. |
| Source-only `CreatePairedObject(df.pair=data.frame(), object=NULL, use.assay1="RNA", use.assay2="ATAC")` inside tryCatch | Reproduces the unused-arguments error visible in #29 screenshot; documentation mismatch is separate from #32. |
| `Rscript scripts/maintenance/reproduce_32.R` | Passed diagnostic assertions: original Assay5 extraction/renaming corrupts pairing order; explicit name reindexing preserves it. |
| `R_LIBS_USER=/tmp/scmega-r-lib SCMEGA_REQUIRE_INTEGRATION=true Rscript scripts/maintenance/run_pairing_tests.R /tmp/scmega-original .` | Pristine source failed: 14 failed expectations, three test errors, two dense-input warnings; real Signac enabled, zero skips. Five then-present cases; an additional aligned/partial-overlap case was subsequently added. |
| `R_LIBS_USER=/tmp/scmega-r-lib SCMEGA_REQUIRE_INTEGRATION=true Rscript scripts/maintenance/run_pairing_tests.R .` | Corrected source: 117 expectations pass in six cases, zero failures/errors/skips, two dense-to-dgCMatrix warnings. Real Signac integration included. |
| `R_LIBS_USER=/tmp/scmega-r-lib R CMD INSTALL --library=/tmp/scmega-final/lib .` | Failed: 12 missing Imports. Signac installation supplied GenomicRanges, reducing initial missing count by one. |
| `R_LIBS_USER=/tmp/scmega-r-lib R CMD build --no-build-vignettes --no-manual .` | Passed. |
| `R_LIBS_USER=/tmp/scmega-r-lib R CMD check --no-manual --no-build-vignettes --output=/tmp/scmega-final /tmp/scmega-final/scMEGA_1.1.0.tar.gz` | Failed: 1 ERROR, 12 missing Imports and undeclared Seurat. Package tests cannot execute through R CMD check. |
| `python3 scripts/maintenance/validate_patch.py --repo <fresh-checkout> --patch <pilot.diff> --base a9bea1ece9f9712e5d6debf7942528d888456548 --output docs/maintenance/evidence/pilot-patch-validation.json`, then trusted pairing runner against that checkout | Actual pilot patch passed policy/application; 117 expectations passed after independent application. Patch SHA-256 and resulting tree archived. This is a real local patch round trip, not a Codex/GitHub hosted run. |
| `python3 -m unittest discover -s tests/maintenance -v` | 22 tests passed; includes many adversarial path/role subcases. No real PR API writes. |
| `/tmp/scmega-actionlint/actionlint -shellcheck= -pyflakes= .github/workflows/codex-issue-fix.yml .github/workflows/R-CMD-check.yml .github/workflows/render-rmarkdown.yaml` | Passed actionlint 1.7.12. Download checksum verified. Shellcheck/pyflakes unavailable and explicitly disabled; this is workflow syntax/semantic lint, not hosted execution. |
| `python3 -m compileall -q scripts/maintenance tests/maintenance` | Passed. |
| Parse all nine R source files and `tools::parse_Rd("man/CreatePairedObject.Rd")` | Passed. |
| `git diff --check` | Passed after trimming trailing whitespace in captured logs (log content otherwise retained). |

Tests use synthetic data only. testthat and Signac were installed into
`/tmp/scmega-r-lib`; BiocManager used `update=FALSE`, with no repository dependency
upgrade or change to the existing library. The exact loaded versions are archived.
Local existing S4/Bioconductor versions are not a newly certified compatibility matrix.

## Dry-run coverage versus hosted verification

| Required scenario | Local evidence | Hosted status |
|---|---|---|
| Issue retrieval | Live read-only REST script and exact audit comparison pass; mock pagination spans multiple pages. | Workflow run not performed. |
| Authorization | API response test doubles cover write/maintain/admin, read/triage, spoofed issue-author role, unauthorized rerun, removed labels and non-default manual refs. | Real event/GITHUB_TOKEN behavior pending. |
| Read-only audit | Retrieval script uses only GET operations; no writes occurred. | Manual audit workflow pending deployment. |
| Synthetic patch generation | Known synthetic textual fix and actual pilot patch pass static validation and Git application; fresh-checkout scientific tests pass. | Actual Codex generation not run; API secret not configured by this task. |
| Draft PR creation | GitHub API double validates draft=true, issue branch, expected tree and no merge/force operations. | No remote PR: explicit authorization absent. |
| Unsafe patch rejection | Traversal, malformed/mismatched headers, binary/empty/oversized patches, protected paths, duplicate diffs, symlink ancestors/leaves, executable and gitlink modes rejected. | Hosted artifact transfer and runner tests pending. |
| Failure recovery | Missing-dependency/failed-test/unresolved reports produce diagnostic-only output; speculative patch reports rejected; actual R installation/check failures recorded. | Hosted job-gating and orphan-branch recovery pending. |
| Scientific regression | Real SeuratObject/Signac synthetic count and metadata invariants pass after fix, fail before. | Linux CI and alternate compatible package versions pending. |
| Prompt injection | Adversarial text preserved as JSON data, protected workflow path rejected, generated PR text cannot mention users/close unrelated issues. | Actual model compliance, network/sudo sandbox adversarial tests NOT run. |
| Artifact provenance | Wrong repository/run/attempt/base/issue/workflow, replay and wrong event rejected in unit tests. | GitHub artifact metadata and transfer end-to-end pending. |

An initial security test incorrectly expected Git to reject a displaced hunk;
Git legitimately applies matching context with offsets. The fixture was corrected
to alter the required original context, and now verifies an actual non-applicable
patch is rejected. No security enforcement was weakened to make that test pass.

## Acceptance status

Complete locally: repository/issue audits, confirmed pilot investigation and minimal
correction, regression fixtures, project instructions, R CI definition, gated Codex
workflow and trusted publication tooling, local security tests, setup/security/
troubleshooting documentation and PR draft.

Not verified/blocked: green package-wide installation/check, supported version
matrix, real hosted authorization/sandbox/injection behavior, model-generated patch
round trip, real artifact provenance transfer, draft PR creation and PR CI triggering.
Maintain both opt-in variables off in production and require review of the named
publication environment. Do not mark the project deployed or the full package
validated based on this session.
