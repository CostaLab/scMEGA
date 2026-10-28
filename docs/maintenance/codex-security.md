# Security review

## Trust boundaries

Issue title/body/comments/images/links and all generated code/reports are untrusted.
JSON encoding prevents shell interpolation; the prompt identifies issue text as
non-authoritative. Neither delimiters nor AGENTS.md are security boundaries.
Authorization comes from current GitHub API permissions of the label actor and
rerun actor, not issue authorship. Only `codex-fix` authorizes modifications. Manual
fix mode needs the same label and permissions, and only runs from the default branch.

Codex receives contents:read, no persisted checkout credential and no write-capable
repository credential. The API key is only supplied through the official action's
proxy input. `:workspace` restricts commands' filesystem/network access;
`drop-sudo` provides process privilege isolation on a disposable runner. Codex runs
last, so no later trusted code relies on its mutable host state. A fresh job consumes
only bounded structured output, not its workspace or altered Git metadata.

Publication uses a fresh runner and trusted scripts checked out from the pinned
default-branch revision. Artifact content is treated as data. No R/test/build/install
command, git hook, generated script, action or configuration is executed in this
write-capable job. Git API publication avoids commit hooks and unsafe Git pushes.
R tests execute only in a separate read-only job without API secrets/write tokens.

## Patch/artifact enforcement

* Only same-run, same-attempt immutable artifact ID is downloaded. API metadata must
  match repository, workflow path/event, run ID, attempt and source revision. Locally
  stamped base/issue/workflow-ref and SHA-256 patch/report digests must match trusted
  job context. Agent claims are never used as provenance.
* Strict positive allowlist: direct `R/*.R`, direct `man/*.Rd`, direct
  `tests/testthat/test-*.R`. No other paths; automatic changes cannot edit the validator.
* Thus all `.github/workflows`, `.github/actions`, CODEOWNERS, AGENTS/overrides,
  DESCRIPTION, NAMESPACE, scripts, credentials/dotfiles, hooks, deployment or
  permission settings are excluded. Sensitive filename patterns are also rejected.
* Noncanonical path aliases (absolute, traversal, empty/dot components, backslash,
  control characters, quotes and Unicode) are rejected before filesystem operations.
  Symlink leaves and ancestors, gitlinks and executable modes are rejected. Only
  regular 100644 UTF-8 files, additions/edits, no deletions/renames/binary patches.
* Size and file-count bounds; a regression-test change is required. Header paths
  must match and must agree with actual staged changed paths. Git apply checks
  context against a clean exact base; matching resulting Git tree is checked at
  publication. Passing static policy is not proof of a correct scientific fix.
* Candidate runs include trusted baseline pairing tests as well as candidate tests
  and package check. Installation/check failure blocks publish. The publisher checks
  the issue remains open/labeled, base has not advanced, and no issue branch/PR exists.
* PRs are draft. There is no force-push, merge, issue-close or status-label mutation.
  Protected-area requests require a separate maintainer-reviewed change.

## Remaining risks and limitations

No hosted sandbox, real unauthorized event, model-level injection test or live PR
creation has been exercised. Publication is disabled; the environment's reviewer
rules must be configured externally. Offline mock tests cannot establish platform
permissions, token behavior or runtime isolation. The official action is pinned,
but CLI/npm supply chain and transitive actions/dependency remotes require review.
R package dependency installation is executable third-party code on a read-only
runner; do not add secrets or writable caches there. Source test coverage is limited.

Allowed R/test files are executable code. A malicious or erroneous patch within
allowed paths can still change scientific results, weaken its own tests, contain
network calls or exfiltrate data if later merged. Trusted baseline tests reduce but
do not eliminate this risk. Human review, protected branches and a green validated
baseline remain necessary. No claim is made that prompt injection is solved.

GitHub job output size limits can truncate/drop large or secret-looking model
outputs; this becomes failed diagnostics, not publication. Patch size is capped.
Artifact digests bind bytes to the run, not truthfulness of the report. Diagnostics
are public repository artifacts and must contain synthetic/public data only.
A default-branch movement after the final API check is a residual race; the PR is
based on its explicitly validated commit and requires review/revalidation before
merge. Git API PR failure can leave an orphan branch, which is never force-updated.

The previous website workflow's automatic commit/push was removed and made manual
read-only because it would otherwise execute on automation branches. Other future
workflows must preserve this boundary. Branch protection/CODEOWNERS are recommended
external controls, not assumed to exist. No repository settings were changed.
