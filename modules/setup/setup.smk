import os

EXPERIMENT_DIR = config["experiment_dir"]
RUN_DIR = config["run_dir"]

localrules: setup_diann_spectral_library_config, setup_diann_run_config, setup_peptidotyping_infinidia_config, setup_peptidotyping_standard_config, setup_hapid_infinidia_config, setup_diann_library_search_base_config, write_run_manifest

################################################################################
# Configuration Setup Rules
################################################################################
# Snapshot the DIA-NN spectral-library base config into RUN_DIR/config/
# for reproducibility. Per-rule digest flags are added inline by the
# consuming rules; see diann_spectral_library_base.cfg for context.
rule setup_diann_spectral_library_config:
    output:
        output_config_file = os.path.join(RUN_DIR,"config/diann_spectral_library_base.cfg")
    params:
        selected_config = config.get("diann_spectral_library_base_config")
    log: os.path.join(RUN_DIR,"logs/setup/setup_diann_spectral_library_config.log")
    run:
        import shutil
        shutil.copy(params.selected_config, output.output_config_file)

# Handle DIANN run config
rule setup_diann_run_config:
    output:
        diann_run_config_file = os.path.join(RUN_DIR,"config/run_diann.cfg")
    params:
        selected_config = (
            "config/run_diann_infinidia.cfg"
            if config.get("diann_search_mode") == "infinidia"
            else config.get("run_diann_config")
        )
    log: os.path.join(RUN_DIR,"logs/setup/setup_diann_run_config.log")
    run:
        import shutil
        shutil.copy(params.selected_config, output.diann_run_config_file)

# Snapshot the peptidotyping InfiniDIA DIA-NN cfg into RUN_DIR/config/.
# Source path is config["peptidotyping_infinidia_config"]; experiments can
# override that key to swap in a variant (e.g. proteoform-mode) cfg.
rule setup_peptidotyping_infinidia_config:
    output:
        output_config_file = os.path.join(RUN_DIR,"config/peptidotyping_infinidia.cfg")
    params:
        selected_config = config.get("peptidotyping_infinidia_config")
    log: os.path.join(RUN_DIR,"logs/setup/setup_peptidotyping_infinidia_config.log")
    run:
        import shutil
        shutil.copy(params.selected_config, output.output_config_file)

# Snapshot the peptidotyping standard-mode DIA-NN cfg into RUN_DIR/config/.
# Source path is config["peptidotyping_standard_config"]; same override story as
# the InfiniDIA cfg above.
rule setup_peptidotyping_standard_config:
    output:
        output_config_file = os.path.join(RUN_DIR,"config/peptidotyping_standard.cfg")
    params:
        selected_config = config.get("peptidotyping_standard_config")
    log: os.path.join(RUN_DIR,"logs/setup/setup_peptidotyping_standard_config.log")
    run:
        import shutil
        shutil.copy(params.selected_config, output.output_config_file)

# Snapshot the HAPiID InfiniDIA DIA-NN cfg into RUN_DIR/config/.
rule setup_hapid_infinidia_config:
    output:
        output_config_file = os.path.join(RUN_DIR,"config/hapid_infinidia.cfg")
    params:
        selected_config = config.get("hapid_infinidia_config")
    log: os.path.join(RUN_DIR,"logs/setup/setup_hapid_infinidia_config.log")
    run:
        import shutil
        shutil.copy(params.selected_config, output.output_config_file)

# Snapshot the DIA-NN library-search base cfg into RUN_DIR/config/.
rule setup_diann_library_search_base_config:
    output:
        output_config_file = os.path.join(RUN_DIR,"config/diann_library_search_base.cfg")
    params:
        selected_config = config.get("diann_library_search_base_config")
    log: os.path.join(RUN_DIR,"logs/setup/setup_diann_library_search_base_config.log")
    run:
        import shutil
        shutil.copy(params.selected_config, output.output_config_file)

# Point-in-time reproducibility snapshot: resolved merged config + git SHA +
# sample_annotation hash + cfg-file hashes + container tags.
#
# Runs inline on the host (no container) — Snakemake's `run:` directive
# executes Python in the Snakemake process itself, which already has all
# the stdlib we need and `git` on PATH. Using `script:` here forced a
# container invocation; the conduitr container has `python3` but not
# `python`, which Snakemake's script runner expects.
rule write_run_manifest:
    input:
        sample_annotation = os.path.join(EXPERIMENT_DIR, config["sample_annotation"])
    output:
        manifest = os.path.join(RUN_DIR, "manifest.json")
    params:
        repo_root = workflow.basedir
    log: os.path.join(RUN_DIR, "logs/setup/write_run_manifest.log")
    run:
        import hashlib
        import json
        import subprocess
        from datetime import datetime, timezone

        def _sha256(path):
            h = hashlib.sha256()
            with open(path, "rb") as f:
                for chunk in iter(lambda: f.read(1 << 16), b""):
                    h.update(chunk)
            return h.hexdigest()

        def _git(args):
            try:
                return subprocess.run(
                    ["git", "-C", params.repo_root] + args,
                    capture_output=True, text=True, check=True, timeout=10,
                ).stdout.strip()
            except Exception:
                return None

        def _safe(v):
            if isinstance(v, dict):
                return {str(k): _safe(x) for k, x in v.items()}
            if isinstance(v, (list, tuple)):
                return [_safe(x) for x in v]
            if isinstance(v, (str, int, float, bool)) or v is None:
                return v
            return str(v)

        _git_sha = _git(["rev-parse", "HEAD"])
        _git_dirty_status = _git(["status", "--porcelain"])
        _git_dirty = bool(_git_dirty_status) if _git_dirty_status is not None else None

        _cfg_snapshot = {}
        for _k, _v in config.items():
            if isinstance(_v, str) and _v.endswith(".cfg"):
                _p = _v if os.path.isabs(_v) else os.path.join(params.repo_root, _v)
                _entry = {"path": _v}
                try:
                    _entry["sha256"] = _sha256(_p)
                except FileNotFoundError:
                    _entry["sha256"] = None
                    _entry["error"] = "file not found"
                _cfg_snapshot[_k] = _entry

        _containers = config.get("containers", {}) or {}
        _container_info = {_n: {"tag": _i} for _n, _i in _containers.items()}

        _manifest = {
            "workflow": {
                "git_sha": _git_sha,
                "git_dirty": _git_dirty,
                "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            },
            "run": {
                "experiment": config.get("experiment"),
                "run_name": config.get("run_name"),
                "search_space_method": config.get("search_space_method"),
                "resolved_config": _safe(dict(config)),
            },
            "inputs": {
                "sample_annotation_path": input.sample_annotation,
                "sample_annotation_sha256": _sha256(input.sample_annotation),
            },
            "cfg_files": _cfg_snapshot,
            "containers": _container_info,
        }

        os.makedirs(os.path.dirname(output.manifest), exist_ok=True)
        with open(output.manifest, "w") as _f:
            json.dump(_manifest, _f, indent=2, sort_keys=True)