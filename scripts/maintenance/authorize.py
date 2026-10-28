"""Authorize using live collaborator permissions, never issue author association."""
import json
import os
from pathlib import Path
from fetch_issues import gh_json

TRUSTED = {"admin", "maintain", "write"}


def allowed(permission):
    return permission in TRUSTED


def authorize(event, event_name, actor, triggering_actor, api):
    repo = event["repository"]["full_name"]
    default = api(f"repos/{repo}")["default_branch"]
    if event_name == "issues":
        if event.get("action") != "labeled" or event.get("label", {}).get("name") != "codex-fix":
            raise ValueError("Only codex-fix label events authorize generation")
        number = event["issue"]["number"]
        mode = "fix"
    elif event_name == "workflow_dispatch":
        if os.environ["GITHUB_REF"] != f"refs/heads/{default}":
            raise ValueError("Manual runs must use the default branch")
        number = int(event["inputs"]["issue"])
        mode = event["inputs"].get("mode", "audit")
        if mode not in {"audit", "fix"}:
            raise ValueError("Invalid mode")
    else:
        raise ValueError("Unexpected event")
    if type(number) is not int or number < 1:
        raise ValueError("Invalid issue")
    for login in {actor, triggering_actor}:
        permission = api(f"repos/{repo}/collaborators/{login}/permission")["permission"]
        if not allowed(permission):
            raise PermissionError("Trigger requires current write, maintain or admin permission")
    issue = api(f"repos/{repo}/issues/{number}")
    if issue["state"] != "open" or "pull_request" in issue:
        raise ValueError("Target must be an open issue")
    if mode == "fix" and "codex-fix" not in {x["name"] for x in issue["labels"]}:
        raise ValueError("Fix mode requires the current codex-fix label")
    sha = api(f"repos/{repo}/commits/{default}")["sha"]
    return dict(issue=str(number), base=sha, default=default, mode=mode)


if __name__ == "__main__":
    event = json.loads(Path(os.environ["GITHUB_EVENT_PATH"]).read_text())
    result = authorize(event, os.environ["GITHUB_EVENT_NAME"], os.environ["GITHUB_ACTOR"],
                       os.environ["GITHUB_TRIGGERING_ACTOR"], lambda p: gh_json("api", p))
    with open(os.environ["GITHUB_OUTPUT"], "a") as out:
        for key, value in result.items():
            if "\n" in value or "\r" in value:
                raise ValueError("Invalid output")
            out.write(f"{key}={value}\n")
