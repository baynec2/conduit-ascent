import os
import sys

EXPERIMENT_DIR = config["experiment_dir"]
RUN_DIR = config["run_dir"]

# The DIA-NN preflight helper (mode resolution + probe parsing). The Snakefile
# already put modules/_shared on sys.path; re-add defensively so this module
# also works when included directly (e.g. from a test harness). Note that
# workflow.basedir is the top-level Snakefile's directory even inside an
# included module, so the path is built from the repo root, not from here.
sys.path.insert(0, os.path.join(workflow.basedir, "modules", "_shared"))
import diann_env
from workflow_version import workflow_version

DIANN_COMPAT_JSON = os.path.join(RUN_DIR, "logs/diann/diann_compatibility.json")
DIANN_PROBE_RAW = os.path.join(RUN_DIR, "logs/diann/diann_compat_probe_raw.txt")

# Every DIA-NN rule in every module reads a cfg snapshot out of RUN_DIR/config/,
# so gating the six snapshot rules below on the compatibility check is enough to
# guarantee it runs before any search — no need to touch 35 DIA-NN rules.
DIANN_COMPAT_GATE = [DIANN_COMPAT_JSON] if config.get("diann_compat_gate") else []

localrules: setup_diann_spectral_library_config, setup_diann_run_config, setup_peptidotyping_infinidia_config, setup_peptidotyping_standard_config, setup_hapid_infinidia_config, setup_diann_library_search_base_config, check_diann_compatibility, write_run_manifest

################################################################################
# DIA-NN compatibility check
################################################################################
# The parse-time preflight (modules/_shared/diann_env.py) only confirms the
# user's DIA-NN exists and has the right shape. These two rules confirm it is a
# DIA-NN this workflow can actually drive.
#
# Why it matters: DIA-NN does not fail on an unknown flag — it prints
# "WARNING: unrecognised option [--flag]" and carries on. So an older build
# silently ignores --pre-search/--pre-filter and runs every InfiniDIA search as
# a plain library search: no error, no crash, wrong science. Verified against
# DIA-NN 2.1.0, which does exactly that.
#
# Split in two because the probe must run inside the DIA-NN container (which
# has no Python) while the parsing must run on the host (which does):
#   probe_diann_cli         -> raw DIA-NN banner + warnings, in the container
#   check_diann_compatibility -> parse, verdict, JSON + human-readable report
rule probe_diann_cli:
    output:
        raw = DIANN_PROBE_RAW
    params:
        probe_cmd = config.get("diann_compat_probe_cmd", ""),
        tmpdir = config.get("diann_compat_probe_tmpdir", ""),
        # Not used by the shell command — carried so that swapping or upgrading
        # the DIA-NN install re-triggers the probe (params is a default
        # rerun trigger); otherwise a stale PASS would stand.
        diann_fingerprint = (config.get("diann") or {}).get("fingerprint", ""),
    log: os.path.join(RUN_DIR, "logs/diann/probe_diann_cli.log")
    # An incompatible DIA-NN fails identically every time; the profile's blanket
    # `retries: 2` (a net for transient network blips) would just print the same
    # verdict three times.
    retries: 0
    container:
        config["containers"]["diann"]
    shell:
        # DIA-NN is handed every flag the workflow uses but no input files, so
        # it reports "0 files will be processed" and exits in a few seconds.
        # `|| true`: a non-zero exit is itself data (the parser reports it),
        # and must not abort the run before we can explain why.
        """
        mkdir -p $(dirname {output.raw}) {params.tmpdir}
        {params.probe_cmd} > {output.raw} 2>&1 || true
        rm -rf {params.tmpdir}
        cp {output.raw} {log}
        """

rule check_diann_compatibility:
    input:
        raw = DIANN_PROBE_RAW
    output:
        compat = DIANN_COMPAT_JSON
    log: os.path.join(RUN_DIR, "logs/diann/check_diann_compatibility.log")
    retries: 0
    run:
        with open(input.raw) as _f:
            _probe_output = _f.read()

        _verdict = diann_env.evaluate_compat_probe(_probe_output)
        _info = config.get("diann", {})
        _txt = diann_env.write_compat_artifacts(output.compat, _verdict, _info)

        _report = diann_env.render_compat_report(_verdict, _info)
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)
        with open(log[0], "w") as _f:
            _f.write(_report)
        print(_report, file=sys.stderr)

        if not _verdict["ok"] and config.get("diann_compat_check") == "strict":
            raise ValueError(
                _report
                + "\nSet diann_compat_check: warn (or off) to proceed anyway, or "
                  "install a DIA-NN version this workflow supports "
                  f"(tested: {', '.join(diann_env.DIANN_TESTED_VERSIONS)}).\n"
                  f"Full report: {_txt}"
            )

################################################################################
# Configuration Setup Rules
################################################################################
# Snapshot the DIA-NN spectral-library base config into RUN_DIR/config/
# for reproducibility. Per-rule digest flags are added inline by the
# consuming rules; see diann_spectral_library_base.cfg for context.
rule setup_diann_spectral_library_config:
    input:
        DIANN_COMPAT_GATE
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
    input:
        DIANN_COMPAT_GATE
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
    input:
        DIANN_COMPAT_GATE
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
    input:
        DIANN_COMPAT_GATE
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
    input:
        DIANN_COMPAT_GATE
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
    input:
        DIANN_COMPAT_GATE
    output:
        output_config_file = os.path.join(RUN_DIR,"config/diann_library_search_base.cfg")
    params:
        selected_config = config.get("diann_library_search_base_config")
    log: os.path.join(RUN_DIR,"logs/setup/setup_diann_library_search_base_config.log")
    run:
        import shutil
        shutil.copy(params.selected_config, output.output_config_file)

# Point-in-time reproducibility snapshot: resolved merged config + git SHA +
# workflow version + sample_annotation hash + cfg-file hashes + container tags.
#
# Runs inline on the host (no container) — Snakemake's `run:` directive
# executes Python in the Snakemake process itself, which already has all
# the stdlib we need and `git` on PATH. Using `script:` here forced a
# container invocation; the conduitr container has `python3` but not
# `python`, which Snakemake's script runner expects.
rule write_run_manifest:
    input:
        sample_annotation = os.path.join(EXPERIMENT_DIR, config["sample_annotation"]),
        # DIA-NN is no longer pinned by an image tag — it is whatever the user
        # installed — so the manifest has to capture the version that actually
        # ran. The compatibility check already parsed it out of DIA-NN's banner;
        # depend on its JSON rather than invoking DIA-NN a second time.
        diann_compat = DIANN_COMPAT_GATE
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

        # DIA-NN provenance: how it was resolved, and which build actually ran.
        _diann_info = dict(config.get("diann", {}) or {})
        _diann_info["compat_check"] = config.get("diann_compat_check")
        _diann_info["version"] = None
        if input.diann_compat:
            try:
                with open(input.diann_compat[0]) as _f:
                    _compat = json.load(_f)
                _diann_info["version"] = _compat.get("version")
                _diann_info["compat_ok"] = _compat.get("ok")
                _diann_info["unsupported_flags"] = _compat.get("unsupported_flags")
            except Exception as _e:  # never let provenance capture fail the run
                _diann_info["compat_error"] = str(_e)

        _manifest = {
            "workflow": {
                "version": workflow_version(params.repo_root),
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
            "diann": _diann_info,
        }

        os.makedirs(os.path.dirname(output.manifest), exist_ok=True)
        with open(output.manifest, "w") as _f:
            json.dump(_manifest, _f, indent=2, sort_keys=True)