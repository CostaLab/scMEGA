"""Verify immutable artifact identity using GitHub's metadata and local digests."""
import hashlib
import json
import os
from pathlib import Path
from fetch_issues import gh_json
from validate_patch import provenance


def verify(metadata, run, manifest, expected, artifact_id):
    provenance(manifest, expected)
    if str(metadata["id"]) != artifact_id or metadata["expired"]:
        raise ValueError("Wrong/expired artifact")
    if metadata["name"] != f"fix-{expected['run_id']}-{expected['run_attempt']}":
        raise ValueError("Unexpected artifact name")
    if str(metadata["workflow_run"]["id"]) != expected["run_id"]:
        raise ValueError("Artifact from another run")
    if (str(run["id"]) != expected["run_id"] or str(run["run_attempt"]) != expected["run_attempt"]
        or run["repository"]["full_name"] != expected["repository"]
        or run["path"] != ".github/workflows/codex-issue-fix.yml"
        or run["event"] not in {"issues", "workflow_dispatch"}):
        raise ValueError("Unexpected source workflow")
    if run["head_sha"] != metadata["workflow_run"]["head_sha"]:
        raise ValueError("Artifact source revision mismatch")


if __name__ == "__main__":
    expected = {k: os.environ[k.upper()] for k in ("repository", "run_id", "run_attempt", "workflow_ref", "base", "issue")}
    folder = Path("bundle")
    if {p.name for p in folder.iterdir()} != {"fix.patch", "report.json", "manifest.json"}:
        raise ValueError("Unexpected artifact files")
    for p in folder.iterdir():
        if p.is_symlink() or not p.is_file() or p.stat().st_size > 350000:
            raise ValueError("Invalid artifact member")
    manifest = json.loads((folder / "manifest.json").read_text())
    artifact = gh_json("api", f"repos/{expected['repository']}/actions/artifacts/{os.environ['ARTIFACT_ID']}")
    run = gh_json("api", f"repos/{expected['repository']}/actions/runs/{expected['run_id']}")
    verify(artifact, run, manifest, expected, os.environ["ARTIFACT_ID"])
    for name, key in (("fix.patch", "sha256"), ("report.json", "report_sha256")):
        if hashlib.sha256((folder / name).read_bytes()).hexdigest() != manifest[key]:
            raise ValueError("Artifact content digest mismatch")
