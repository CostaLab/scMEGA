"""Read-only, paginated GitHub issue evidence. Python standard library + gh."""
import argparse
import datetime
import json
import re
import subprocess
from pathlib import Path


def gh_json(*args):
    return json.loads(subprocess.check_output(["gh", *args], text=True))


def fetch(repo, number=None):
    if not re.fullmatch(r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+", repo):
        raise ValueError("Invalid repository")
    if number is not None and (type(number) is not int or number < 1):
        raise ValueError("Invalid issue number")
    if number is None:
        pages = gh_json("api", "--paginate", "--slurp", f"repos/{repo}/issues?state=open&per_page=100")
        numbers = [i["number"] for page in pages for i in page if "pull_request" not in i]
    else:
        numbers = [number]
    issues = []
    for n in numbers:
        issue = gh_json("api", f"repos/{repo}/issues/{n}")
        if "pull_request" in issue:
            raise ValueError("Pull requests are not issue-fix targets")
        comments = gh_json("api", "--paginate", "--slurp", f"repos/{repo}/issues/{n}/comments?per_page=100")
        issues.append({k: issue[k] for k in ("number", "title", "body", "state", "labels", "html_url", "created_at", "updated_at")})
        issues[-1]["comments"] = [c for page in comments for c in page]
    return {"repository": repo, "retrieved_at": datetime.datetime.now(datetime.timezone.utc).isoformat(), "issues": issues}


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--repo", default="CostaLab/scMEGA")
    p.add_argument("--issue", type=int)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    args.output.write_text(json.dumps(fetch(args.repo, args.issue), indent=2) + "\n")
