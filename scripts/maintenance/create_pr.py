"""Publish only validated text through the GitHub Git API; never run candidate code."""
import argparse
import base64
import hashlib
import json
import os
import re
import subprocess
from pathlib import Path
from bundle_result import FIELDS
from validate_patch import git, provenance, safe_path


def api(method, endpoint, payload=None):
    command = ["gh", "api", "--method", method, endpoint]
    data = None
    if payload is not None:
        command += ["--input", "-"]
        data = json.dumps(payload)
    return json.loads(subprocess.check_output(command, input=data, text=True))


def pr_body(report, number):
    if report["status"] != "fixed" or report["reproduced"] is not True:
        raise ValueError("Diagnostic only: no fix PR")
    sections = ["## Related Issue\nFixes #" + str(number)]
    for key in FIELDS:
        text = report[key]
        if not isinstance(text, str) or not text.strip() or len(text) > 12000:
            raise ValueError("Invalid report section")
        # Do not let generated text mention people, close other issues, or render HTML.
        text = text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace("@", "&#64;").replace("#", "&#35;")
        sections.append("## " + key.replace("_", " ").title() + "\n" + text)
    sections.append("## Independent validation\nRequired unprivileged unit tests, installation and package check passed in the source run. Maintainer review is required; no automatic merge.")
    return "\n\n".join(sections)


def publish(repo, checkout, manifest, validation, report, default, request=api):
    if not re.fullmatch(r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+", repo):
        raise ValueError("Invalid repository")
    n = int(manifest["issue"])
    if n < 1:
        raise ValueError("Invalid issue")
    base = manifest["base"]
    if validation["base"] != base or validation["sha256"] != manifest["sha256"]:
        raise ValueError("Validation does not match patch")
    body = pr_body(report, n)
    branch = f"codex/issue-{n}"
    # Idempotent fail-closed behavior: do not overwrite any existing branch or PR.
    pulls = request("GET", f"repos/{repo}/pulls?state=all&head={repo.split('/')[0]}:{branch}")
    refs = request("GET", f"repos/{repo}/git/matching-refs/heads/{branch}")
    if pulls or any(r["ref"] == "refs/heads/" + branch for r in refs):
        raise ValueError("Existing branch/PR: maintainer intervention required")
    current = request("GET", f"repos/{repo}/commits/{default}")
    if current["sha"] != base:
        raise ValueError("Default branch advanced; rerun against new base")
    issue = request("GET", f"repos/{repo}/issues/{n}")
    if issue["state"] != "open" or "pull_request" in issue or "codex-fix" not in {x["name"] for x in issue["labels"]}:
        raise ValueError("Issue authorization was withdrawn")
    paths = validation["paths"]
    entries = []
    for name in paths:
        safe_path(name)
        # Read from the validated Git index, not a potentially replaced filesystem file.
        content = git(checkout, "show", ":" + name)
        blob = request("POST", f"repos/{repo}/git/blobs", {"content": base64.b64encode(content).decode(), "encoding": "base64"})
        entries.append({"path": name, "mode": "100644", "type": "blob", "sha": blob["sha"]})
    tree = request("POST", f"repos/{repo}/git/trees", {"base_tree": current["commit"]["tree"]["sha"], "tree": entries})
    if tree["sha"] != validation["tree"]:
        raise ValueError("Published tree differs from independently validated tree")
    commit = request("POST", f"repos/{repo}/git/commits", {"message": f"fix: resolve issue #{n}", "tree": tree["sha"], "parents": [base]})
    request("POST", f"repos/{repo}/git/refs", {"ref": "refs/heads/" + branch, "sha": commit["sha"]})
    return request("POST", f"repos/{repo}/pulls", {"title": f"fix: resolve issue #{n}",
        "head": branch, "base": default, "body": body, "draft": True})


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--checkout", required=True)
    p.add_argument("--bundle", required=True)
    p.add_argument("--validation", required=True)
    p.add_argument("--default", required=True)
    a = p.parse_args()
    folder = Path(a.bundle)
    manifest = json.loads((folder / "manifest.json").read_text())
    expected = {k: os.environ[k.upper()] for k in ("repository", "run_id", "run_attempt", "workflow_ref", "base", "issue")}
    provenance(manifest, expected)
    for filename, key in (("fix.patch", "sha256"), ("report.json", "report_sha256")):
        if hashlib.sha256((folder / filename).read_bytes()).hexdigest() != manifest[key]:
            raise ValueError("Artifact digest mismatch")
    print(json.dumps(publish(expected["repository"], a.checkout, manifest,
        json.loads(Path(a.validation).read_text()), json.loads((folder / "report.json").read_text()), a.default)))
