"""The version string a run records as having been produced by.

`VERSION` holds the release number, but every commit between releases carries
the same number, so on its own it cannot tell a release from development code.
The string is the bare number only when HEAD is the commit tagged `v<number>`
with no uncommitted changes; anything else gets the commit appended as SemVer
build metadata (`0.1.0+e81487c`, `0.1.0+e81487c.dirty`). Without git or a
checkout the number is returned unchanged, since nothing better is known.
"""
import os
import subprocess


def _git(repo_root, args):
    try:
        return subprocess.run(
            ["git", "-C", repo_root] + args,
            capture_output=True, text=True, check=True, timeout=10,
        ).stdout.strip()
    except Exception:
        return None


def workflow_version(repo_root):
    with open(os.path.join(repo_root, "VERSION")) as f:
        number = f.read().strip()

    sha = _git(repo_root, ["rev-parse", "--short=7", "HEAD"])
    if sha is None:
        return number

    tags = (_git(repo_root, ["tag", "--points-at", "HEAD"]) or "").split()
    dirty = bool(_git(repo_root, ["status", "--porcelain"]))
    if f"v{number}" in tags and not dirty:
        return number
    return f"{number}+{sha}" + (".dirty" if dirty else "")
