"""Fail-closed text patch policy. Never imports or executes candidate code."""
import argparse
import hashlib
import json
import os
import re
import subprocess
from pathlib import Path

MAX_PATCH = 262144
MAX_FILES = 20


def git(repo, *args):
    env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
    env.update(GIT_CONFIG_NOSYSTEM="1", GIT_CONFIG_GLOBAL=os.devnull)
    return subprocess.check_output(["git", "-c", "core.hooksPath=/dev/null", "-c",
        "core.attributesFile=/dev/null", "-C", str(repo), *args], env=env)


def safe_path(name):
    # Reject aliases before normalization: no dot components, escapes, quoting, Unicode,
    # control characters, Windows separators, absolute paths, or hidden components.
    if not re.fullmatch(r"[A-Za-z0-9_/-]+(?:\.[A-Za-z0-9_-]+)*", name):
        raise ValueError("Noncanonical path")
    parts = name.split("/")
    if any(p in {"", ".", ".."} or p.startswith(".") for p in parts):
        raise ValueError("Unsafe path components")
    low = name.lower()
    if any(word in low for word in ("secret", "credential", "token", "password", "private-key")):
        raise ValueError("Sensitive path")
    # Narrow positive allowlist also protects DESCRIPTION/NAMESPACE, nested AGENTS,
    # workflow/action files, hooks, credentials, deployment and maintenance scripts.
    permitted = (re.fullmatch(r"R/[A-Za-z0-9_-]+\.R", name)
        or re.fullmatch(r"man/[A-Za-z0-9_-]+\.Rd", name)
        or re.fullmatch(r"tests/testthat/test-[A-Za-z0-9_-]+\.R", name))
    if not permitted:
        raise ValueError("Protected or unsupported path: " + name)
    return name


def patch_paths(data):
    if not data or len(data) > MAX_PATCH or b"\x00" in data or b"\r" in data:
        raise ValueError("Empty, oversized or non-text patch")
    text = data.decode("utf-8", errors="strict")
    if not text.startswith("diff --git ") or not text.endswith("\n"):
        raise ValueError("Expected canonical git diff")
    paths = []
    for section in text.split("diff --git ")[1:]:
        lines = section.splitlines()
        match = re.fullmatch(r"a/(\S+) b/(\S+)", lines[0])
        if not match or match[1] != match[2]:
            raise ValueError("Renames/quoted/malformed paths forbidden")
        name = safe_path(match[1])
        if name in paths:
            raise ValueError("Duplicate diff path")
        paths.append(name)
        # Inspect only metadata before the first hunk; source lines are untrusted data.
        header = []
        for line in lines[1:]:
            if line.startswith("@@ "):
                break
            header.append(line)
        else:
            raise ValueError("No textual hunk")
        for line in header:
            if line.startswith("index "):
                if not re.fullmatch(r"index [0-9a-f]+\.\.[0-9a-f]+(?: 100644)?", line):
                    raise ValueError("Invalid index/mode")
            elif line not in {"new file mode 100644", "--- /dev/null", "--- a/" + name, "+++ b/" + name}:
                raise ValueError("Deletion, binary, symlink, executable or unsupported header")
        if "+++ b/" + name not in header or not ({"--- /dev/null", "--- a/" + name} & set(header)):
            raise ValueError("Mismatched patch headers")
    if not 1 <= len(paths) <= MAX_FILES:
        raise ValueError("File limit")
    if not any(p.startswith("tests/testthat/test-") for p in paths):
        raise ValueError("A regression test change is required")
    return paths


def validate(repo, patch, base):
    repo = Path(repo).resolve()
    if not re.fullmatch(r"[0-9a-f]{40}", base):
        raise ValueError("Invalid base SHA")
    if git(repo, "rev-parse", "HEAD").decode().strip() != base:
        raise ValueError("Base SHA mismatch")
    if git(repo, "status", "--porcelain", "--untracked-files=all"):
        raise ValueError("Checkout must be clean")
    data = Path(patch).read_bytes()
    paths = patch_paths(data)
    for name in paths:
        current = repo
        for part in name.split("/"):
            current = current / part
            if current.is_symlink():
                raise ValueError("Symlink path component")
        entry = git(repo, "ls-files", "--stage", "--", name).decode()
        if entry and not entry.startswith("100644 "):
            raise ValueError("Unsafe original mode")
    patch = str(Path(patch).resolve())
    git(repo, "apply", "--check", "--index", "--whitespace=error-all", patch)
    git(repo, "apply", "--index", "--whitespace=error-all", patch)
    changed = git(repo, "diff", "--cached", "--name-only", "-z").decode().split("\0")[:-1]
    if set(changed) != set(paths):
        raise ValueError("Parsed paths differ from actual changes")
    for name in changed:
        entry = git(repo, "ls-files", "--stage", "--", name).decode()
        if not entry.startswith("100644 ") or not (repo / name).is_file():
            raise ValueError("Only regular non-executable text files may be published")
        content = (repo / name).read_bytes()
        content.decode("utf-8", errors="strict")
        if len(content) > MAX_PATCH or b"\0" in content:
            raise ValueError("Unsafe file content/size")
    return {"base": base, "paths": sorted(changed), "sha256": hashlib.sha256(data).hexdigest(),
            "tree": git(repo, "write-tree").decode().strip()}


def provenance(manifest, expected):
    for key, value in expected.items():
        if manifest.get(key) != value:
            raise ValueError("Artifact provenance mismatch: " + key)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--repo", required=True)
    p.add_argument("--patch", required=True)
    p.add_argument("--base", required=True)
    p.add_argument("--output", required=True)
    a = p.parse_args()
    Path(a.output).write_text(json.dumps(validate(a.repo, a.patch, a.base), indent=2) + "\n")
