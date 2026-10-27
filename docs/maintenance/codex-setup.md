# Codex maintenance setup

## Current state

The local implementation is **not enabled** and has not run on GitHub. No secret,
label, environment, setting, branch or PR has been created remotely. Local security
and synthetic tests are not evidence of a successful hosted end-to-end run.
See [validation-report.md](validation-report.md) and [repository-audit.md](repository-audit.md).
Package installation/check currently has pre-existing dependency/metadata failures;
these block automatic publication by design.

## Architecture and operation

There are two security boundaries, implemented as five jobs to avoid executing
anything after Codex on its runner:

1. **authorize** (read-only): validate current permissions for both original actor
   and rerun actor through GitHub's collaborator API (write/maintain/admin), verify
   an open issue and current `codex-fix` label, record immutable default-branch SHA,
   retrieve all comments, and archive issue evidence. Manual `audit` mode needs no
   label/key and performs only this read-only work.
2. **generate** (contents:read): trusted checkout, no persisted Git credentials,
   minimal R dependencies installed before Codex, restricted `:workspace` profile,
   `drop-sudo` on disposable Ubuntu 24.04, API key supplied only to the official
   action's proxy input. Codex is the last step. Its bounded JSON output contains
   the patch and diagnostic sections. No repository write token is available.
3. **bundle** (fresh read-only runner): trusted script serializes the untrusted JSON
   as immutable patch/report artifacts and stamps run/repository/base/attempt
   identity. Invalid/empty agent output becomes a failed diagnostic artifact.
4. **validate** (fresh read-only runner, no secrets): same-run artifact provenance,
   digest and path checks; apply to a clean pinned checkout; run trusted baseline
   pairing tests against candidate source, candidate tests, installation and package
   checks. No passing application-only shortcut. No dependency cache is saved.
5. **publish** (fresh runner, protected environment): requires independent validation
   success and publication opt-in. Rechecks the artifact and patch against a clean
   base. It uses only trusted Python/git inspection and GitHub's Git API to create
   blobs, a commit, a new issue branch and a **draft** PR. It never loads R, generated
   scripts or test files. Existing branches/PRs and moved default branch fail closed.

Concurrency serializes by issue number and publication refuses an existing issue
branch/PR (including closed PRs). There are no force updates, merges or issue closes.
Repeated label events may perform redundant diagnostics but cannot overwrite a PR.
Branch: `codex/issue-N`. Failure leaves logs/artifacts, not a speculative PR. A failed
PR API call may leave a branch; maintainers reconcile it before retrying.

## Maintainer configuration, after reviewing the local changes

1. Merge the infrastructure and reviewed pilot through normal review. The issue
   workflow only executes from the default branch. Use a sandbox repository first.
2. Repair declared dependencies/imports in a separately reviewed change. Establish
   a green supported baseline and an environment manifest. Existing README/remotes
   use historical branches; neither successful installation nor broad compatibility
   is assumed by these workflows.
3. Create repository labels: `codex-fix`, `codex-investigating`, `codex-pr-created`,
   `codex-needs-info`, `codex-failed`. Only `codex-fix` triggers generation. The current
   implementation deliberately does not mutate status labels; apply them manually
   from diagnostics. No issues:write token is granted.
4. Add repository Actions secret **OPENAI_API_KEY** with a dedicated project budget
   and least necessary API access. Never put it in YAML, a repository file, general
   job env, or the publication job. No PAT is needed.
5. Create GitHub Environment **codex-publication**, require trusted reviewers,
   prevent self-review where available, restrict deployment branches to the default
   branch, and prohibit administrator bypass where your plan supports it. Merely
   naming an environment in YAML does not configure required reviewers.
6. Restrict default-branch changes and protect workflow/control files through branch
   rules and CODEOWNERS in a separate maintainer change. Set Actions defaults to
   read-only. Review allowed actions and fork-PR rules. Permit GitHub Actions to
   create PRs in repository settings only when ready.
7. Leave `CODEX_FIX_ENABLED` and `CODEX_PUBLISH_ENABLED` absent/false initially.
   Run manual `audit` on an open issue and verify the archived title/body/comments.
   Complete the dry-run checklist below in a sandbox repository before enabling
   unattended generation. Temporarily opt into supervised sandbox generation only
   when the key and reviewers are ready. Production publication remains off.
8. After live security testing and green R validation, set `CODEX_FIX_ENABLED=true`.
   A trusted maintainer can apply `codex-fix` or dispatch `fix` for an already-labeled
   issue. Set `CODEX_PUBLISH_ENABLED=true` only after publication testing; the protected
   environment still requires reviewer approval. Never disable validation to get a PR.

When auto-PRs use GITHUB_TOKEN, their creation normally does not trigger another
pull_request workflow. The required pre-publication validate job supplies the CI
checks; explicitly run R-CMD-check via workflow_dispatch on the PR branch (or use a
separately reviewed GitHub App trigger design). Do not treat missing PR checks as
passed or bypass branch protection. No personal token workaround is embedded here.

## Versions and compatible environments

Official action `openai/codex-action@v1` resolves to commit
`bdf19a4a223ec2549a3e2274a0cf61556bc07675` (v1.13) on audit date. Its annotated tag
object is **not** the commit SHA. Inputs verified against action.yml/README/security
at that commit: `openai-api-key`, `permission-profile`, `safety-strategy`,
`prompt-file`, `output-schema-file`. The action's inputs/security guidance were
reviewed; this is not a complete third-party code audit. CLI/proxy are both pinned to npm version `0.164.0` (resolved on audit date).
Their installation and all transitive action dependencies still require supply-chain review.

r-lib/actions v2.14.0 is pinned to `f9a764fea8d5c63df6ef9a5c7795bf7deb5d7e05`;
checkout/upload/download actions are pinned in YAML. The single candidate CI
configuration is Linux R 4.5 with Seurat v5-oriented dependencies. R >=4.3 is the
package's declared constraint, not proof all versions work. Bioconductor must match
R; do not combine arbitrary Bioconductor/Matrix/Seurat versions. Expand only after
a baseline passes. Dependencies are provisioned without requested upgrades and
without reusable cache; the current resolver/remotes are not a reproducibility lock.
No large biological data are downloaded by the unit fixtures.

## Dry-run acceptance before production enablement

* Authenticated read-only audit matches live issue/comment pagination.
* Reader/triage/outsider and unauthorized rerun actors cannot reach Codex/key access;
  writer/maintainer label events can. Check real API failures as well as roles.
* Synthetic issue produces a valid tested patch; unresolved/missing-dependency and
  failing-test cases produce diagnostics and **no** publication.
* Protected path, traversal, quoted/Windows paths, symlinks/gitlinks/executable modes,
  malformed patches, artifact replay and wrong workflow identity are rejected.
* Adversarial issue content requests workflow/secret/unrelated changes; verify actual
  model behavior and runtime network/sudo restrictions. The offline text-boundary
  test alone does not establish model-level injection resistance.
* Valid patch creates a draft PR on a dedicated branch; retries do not overwrite it.
* Cell-order/count-integrity test passes on a fresh supported R environment.
* Record run URLs, artifact IDs/digests, exact environment and outcomes in the
  validation report. Until these checks pass, keep unattended production fixes off.

Sources: [official Codex Action](https://github.com/openai/codex-action),
[security guidance](https://github.com/openai/codex-action/blob/bdf19a4a223ec2549a3e2274a0cf61556bc07675/docs/security.md),
[official Codex documentation](https://developers.openai.com/codex/github-action),
[r-lib actions](https://github.com/r-lib/actions),
[GitHub workflow trigger limitations](https://docs.github.com/en/actions/how-tos/write-workflows/choose-when-workflows-run/trigger-a-workflow).
