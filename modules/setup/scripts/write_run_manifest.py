"""Write RUN_DIR/manifest.json for reproducibility.

Captures the resolved merged config, git SHA of the workflow, timestamps,
SHA256 of sample_annotation, and hashes of every .cfg file referenced by a
`*_config` key in the resolved config. Used as a point-in-time snapshot of
everything that defined the run.

Invoked via Snakemake's `script:` directive — has access to `snakemake.input`,
`snakemake.output`, `snakemake.params`, `snakemake.config`.
"""
import hashlib
import json
import os
import subprocess
from datetime import datetime, timezone


def sha256_of_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def run_git(repo_root, args):
    try:
        result = subprocess.run(
            ["git", "-C", repo_root] + args,
            capture_output=True,
            text=True,
            check=True,
            timeout=10,
        )
        return result.stdout.strip()
    except Exception:
        return None


def json_safe(value):
    """Best-effort conversion of arbitrary config values to JSON-serializable form."""
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


repo_root = snakemake.params.repo_root
output_path = snakemake.output.manifest
config = dict(snakemake.config)
sample_annotation_path = snakemake.input.sample_annotation

git_sha = run_git(repo_root, ["rev-parse", "HEAD"])
git_dirty_status = run_git(repo_root, ["status", "--porcelain"])
git_dirty = bool(git_dirty_status) if git_dirty_status is not None else None

cfg_snapshot = {}
for key, value in config.items():
    if isinstance(value, str) and value.endswith(".cfg"):
        p_abs = value if os.path.isabs(value) else os.path.join(repo_root, value)
        entry = {"path": value}
        try:
            entry["sha256"] = sha256_of_file(p_abs)
        except FileNotFoundError:
            entry["sha256"] = None
            entry["error"] = "file not found"
        cfg_snapshot[key] = entry

containers = config.get("containers", {}) or {}
container_info = {name: {"tag": image} for name, image in containers.items()}

manifest = {
    "workflow": {
        "git_sha": git_sha,
        "git_dirty": git_dirty,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    },
    "run": {
        "experiment": config.get("experiment"),
        "run_name": config.get("run_name"),
        "search_space_method": config.get("search_space_method"),
        "resolved_config": json_safe(config),
    },
    "inputs": {
        "sample_annotation_path": sample_annotation_path,
        "sample_annotation_sha256": sha256_of_file(sample_annotation_path),
    },
    "cfg_files": cfg_snapshot,
    "containers": container_info,
}

os.makedirs(os.path.dirname(output_path), exist_ok=True)
with open(output_path, "w") as f:
    json.dump(manifest, f, indent=2, sort_keys=True)
