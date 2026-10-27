# Troubleshooting

| Symptom | Meaning / safe response |
|---|---|
| `gh auth status` fails in restricted execution | Verify network access separately. This session succeeded after sandbox network authorization. Never print raw tokens. |
| No generation after labeling | Confirm exact `codex-fix` label, current write/maintain/admin rights, default-branch workflow, open issue, and CODEX_FIX_ENABLED. Audit mode intentionally generates no patch. |
| Authorization API 403/404 | Check repository access and token metadata permissions; fail closed. Do not trust author_association or add `allow-users: '*'`. |
| R setup/install/check fails | Inspect logs and environment versions. Current package has undeclared Seurat and missing imports/dependencies. Repair metadata in a separately reviewed change; do not disable checks. |
| Matrix `as_cholmod_sparse` / irlba error | Capture sessionInfo and library paths; reproduce with a coherent compiled dependency stack. Rebuild linked dependencies as needed; do not change SVD math to bypass the error. |
| Multiple counts/data layers warning | Inspect Layers and cell names. Do not assume first layer contains all cells or discard warnings. Joining layers needs scientific review. |
| Codex network/sudo denied | Expected sandbox controls; provision trusted dependencies before Codex. Never change to unsafe/danger-full-access. |
| Missing/invalid JSON result | Inspect generate logs, output size limits and schema. A failed diagnostic artifact is expected, with no PR. |
| Protected path / symlink / malformed patch | Do not weaken the allowlist. Use a separate maintainer-reviewed process for metadata/config changes. |
| Artifact provenance/digest mismatch | Stop, preserve run/artifact IDs, investigate substitution/replay. Do not manually override identity fields. |
| Tests skip Signac | Focused extraction tests do not validate ChromatinAssay construction. Install missing dependency into an isolated library and rerun; do not count a skip as a pass. |
| Publication skipped | Confirm independent validation succeeded, CODEX_PUBLISH_ENABLED, and protected environment approval. Never bypass failed validation. |
| Existing branch/PR | Inspect the existing `codex/issue-N`; no overwrite or force-push is attempted. Reconcile through maintainer review before a new run. |
| Default branch advanced | Rerun against current default SHA; never silently rebase a validated patch and claim the same test result. |
| Draft PR has no automatic CI event | GITHUB_TOKEN-created PRs may not trigger pull_request workflows. Use manual R-CMD-check on the branch; retain pre-publication validation evidence. |

Local commands:

```sh
python3 scripts/maintenance/fetch_issues.py --repo CostaLab/scMEGA --output /tmp/issues.json
python3 -m unittest discover -s tests/maintenance -v
Rscript scripts/maintenance/reproduce_32.R
Rscript scripts/maintenance/run_pairing_tests.R .
R CMD INSTALL .
R CMD build --no-build-vignettes --no-manual .
R CMD check --no-manual --no-build-vignettes scMEGA_1.1.0.tar.gz
```

For this session's isolated packages, prefix R commands with
`R_LIBS_USER=/tmp/scmega-r-lib`. Do not assume that path exists on another machine.
Record commands, exit codes, skip counts and sessionInfo. Logs under evidence/
show what actually ran. No dataset, secret or repository-wide dependency upgrade is
needed to run the synthetic accessor reproduction.
