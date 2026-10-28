"""Convert bounded untrusted model output to an artifact; no code evaluation."""
import hashlib
import json
import os
from pathlib import Path

FIELDS = ("problem", "root_cause", "changes", "regression_tests", "validation", "scientific_impact", "limitations")


def bundle(raw, out, identity):
    if len(raw.encode()) > 350000:
        raise ValueError("Result too large")
    result = json.loads(raw)
    if result.get("status") not in {"fixed", "needs-info", "failed"}:
        raise ValueError("Invalid result status")
    if type(result.get("reproduced")) is not bool:
        raise ValueError("Missing reproduction status")
    for key in (*FIELDS, "patch"):
        if not isinstance(result.get(key), str):
            raise ValueError("Missing report field: " + key)
    if result["status"] == "fixed" and (not result["reproduced"] or not result["patch"]):
        raise ValueError("Unverified fix claim")
    if result["status"] != "fixed" and result["patch"]:
        raise ValueError("Diagnostic results cannot carry speculative patches")
    out.mkdir(exist_ok=True)
    patch = result.pop("patch").encode()
    (out / "fix.patch").write_bytes(patch)
    (out / "report.json").write_text(json.dumps(result, indent=2) + "\n")
    manifest = dict(identity, sha256=hashlib.sha256(patch).hexdigest())
    manifest["report_sha256"] = hashlib.sha256((out / "report.json").read_bytes()).hexdigest()
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return result["status"]


if __name__ == "__main__":
    identity = {key: os.environ[key.upper()] for key in
        ("repository", "run_id", "run_attempt", "workflow_ref", "base", "issue")}
    try:
        status = bundle(os.environ["CODEX_RESULT"], Path("bundle"), identity)
    except (ValueError, KeyError, TypeError) as exc:
        diagnostic = {key: "Not verified: generation did not produce a valid result." for key in FIELDS}
        diagnostic.update(status="failed", reproduced=False, patch="", limitations=str(exc))
        status = bundle(json.dumps(diagnostic), Path("bundle"), identity)
    with open(os.environ["GITHUB_OUTPUT"], "a") as f:
        f.write(f"status={status}\n")
