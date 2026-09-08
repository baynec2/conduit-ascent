"""Resolve the user-supplied DIA-NN installation and check it against the workflow.

Why this module exists
----------------------
DIA-NN's licence (LICENSE.txt §3, §5) allows one backup copy and forbids
renting/leasing/lending/sublicensing. Publishing a public registry image that
bundles the binary is neither, so conduit-ascent stopped shipping DIA-NN
(issue #63). Each user now obtains DIA-NN themselves, accepts its terms, and
points the workflow at their copy with a single config key::

    diann_path: resources/diann/diann-2.5.0

One key, three autodetected shapes:

===== ============================== ==================================
Mode  Detected when                  Behaviour
===== ============================== ==================================
A     directory containing           run inside the ``diann_runtime``
      ``diann-linux``                image with the directory bind-mounted
B     an executable file             run directly on the host, no container
C     ``docker://…``/``…​.sif``       use as the container, exactly as before
===== ============================== ==================================

Mode A is the recommended shape: the user unzips the vendor archive and stops
there, and the .NET 8 / libgomp dependencies still come from an image so a bare
host works. Mode B suits ``DIA-NN.AppImage`` (self-contained) or a host that
already has .NET 8. Mode C is the migration path for anyone who already built
their own private full image — nobody is broken by this change.

The module has two jobs, split by cost:

``resolve_diann_environment``
    Parse-time, filesystem-only. Classifies the mode, validates the shape, and
    writes ``config["diann_cmd"]`` / ``config["containers"]["diann"]`` /
    ``config["diann"]`` so every rule can stay written as it was. Fails fast
    with a message that tells the user what to download and where to put it.

``build_compat_probe_argv`` / ``evaluate_compat_probe``
    Run-time, executed by the ``check_diann_compatibility`` rule (which runs
    inside whatever container the resolved mode selected). Confirms the user's
    DIA-NN is actually one this workflow can drive: its version is in range and
    it recognises every CLI flag the workflow passes. The report-column half of
    the contract is checked separately, after a search, by
    ``modules/diann/scripts/check_diann_report_columns.R``.
"""

import json
import os
import re
import sys

# ---------------------------------------------------------------------------
# What the user has to download
# ---------------------------------------------------------------------------
DIANN_DOWNLOAD_URL = "https://github.com/vdemichev/DiaNN/releases"
DIANN_LICENCE_URL = "https://github.com/vdemichev/DiaNN/blob/master/LICENSE.txt"
DEFAULT_DIANN_PATH = "resources/diann/diann-2.5.0"
DIANN_PATH_ENV_VAR = "CONDUIT_DIANN_PATH"

#: The executable inside an extracted DIA-NN directory (mode A).
DIANN_LINUX_BINARY = "diann-linux"

#: Container-URI prefixes that mean "this is mode C, an image".
_IMAGE_PREFIXES = ("docker://", "shub://", "library://", "oras://", "http://", "https://")


# ---------------------------------------------------------------------------
# The compatibility contract: what this workflow needs DIA-NN to support
# ---------------------------------------------------------------------------
# Recorded in the compatibility JSON so an old artifact can be told apart from
# one produced by the current probe. Bump when the flag list or the verdict
# rules change meaningfully.
COMPAT_PROBE_VERSION = 1

#: Every DIA-NN CLI flag the workflow passes, with an inert value where the
#: flag takes one. Sourced from the .cfg files in config/ and the inline flags
#: in the DIA-NN rules across modules/. DIA-NN has no --help and no --version;
#: what it *does* have is a hard guarantee that it prints
#: "WARNING: unrecognised option [--flag]" for anything it does not know. So
#: the probe hands it this whole list with no input files ("0 files will be
#: processed"), and any flag echoed back as unrecognised is a flag this DIA-NN
#: build cannot honour. This is how we catch, e.g., DIA-NN 2.1.0 — which
#: silently ignores --pre-search/--pre-filter and would therefore run every
#: InfiniDIA search as a plain library search.
DIANN_REQUIRED_FLAGS = [
    # --- quantification / FDR (config/run_diann.cfg, *_base.cfg) ---
    ("--qvalue", "0.01"),
    ("--min-corr", "2.0"),
    ("--corr-diff", "1.0"),
    ("--time-corr-only", "true"),
    ("--unimod4", "true"),
    ("--mass-acc", "10"),
    ("--mass-acc-ms1", "4"),
    # --- digest (inline per rule; values here are only placeholders) ---
    ("--met-excision", "true"),
    ("--cut", "K*,R*"),
    ("--missed-cleavages", "1"),
    ("--min-pep-len", "7"),
    ("--max-pep-len", "30"),
    # --- library generation (config/diann_spectral_library_base.cfg) ---
    ("--predictor", None),
    ("--fasta-search", None),
    ("--gen-spec-lib", None),
    ("--min-fr-mz", "200"),
    ("--max-fr-mz", "1800"),
    ("--min-pr-mz", "300"),
    ("--max-pr-mz", "1800"),
    ("--min-pr-charge", "2"),
    ("--max-pr-charge", "3"),
    # --- search behaviour ---
    ("--rt-profiling", None),
    ("--species-ids", None),
    # --- InfiniDIA: the peptidotyping / HAPiID broad-database strategy ---
    ("--pre-search", None),
    ("--pre-filter", None),
    # --- taxonomic FDR needs unfiltered PSMs and decoys ---
    ("--no-prot-inf", None),
    ("--report-decoys", None),
    # --- plumbing ---
    ("--threads", "1"),
    ("--verbose", "1"),
]

#: A flag DIA-NN will never implement. The probe passes it too: if DIA-NN does
#: NOT report this one as unrecognised, the detection mechanism itself is
#: broken (message reworded upstream, output swallowed) and a clean probe run
#: proves nothing. Fail loudly rather than issue a false PASS.
COMPAT_SENTINEL_FLAG = "--conduit-ascent-sentinel-not-a-real-flag"

#: DIA-NN versions this workflow is known to drive. Anything below the minimum
#: is rejected; anything at or above it but outside `tested` is a warning.
DIANN_MIN_VERSION = (2, 2)
DIANN_TESTED_VERSIONS = ("2.3.0", "2.5.0")

#: Columns each consumer needs in a DIA-NN report. Checked after a search by
#: modules/diann/scripts/check_diann_report_columns.R — the CLI probe, which
#: runs no files, cannot see these.
DIANN_REQUIRED_REPORT_COLUMNS = {
    # conduitR::diann_to_qfeatures + modules/diann/scripts/extract_detected_proteins.R
    "main": [
        "Run",
        "Precursor.Id",
        "Modified.Sequence",
        "Stripped.Sequence",
        "Precursor.Charge",
        "Precursor.Mz",
        "Proteotypic",
        "Protein.Ids",
        "Protein.Group",
        "Protein.Names",
        "Genes",
        "Q.Value",
        "Global.Q.Value",
        "Lib.Q.Value",
        "PG.Q.Value",
        "Lib.PG.Q.Value",
        "Global.PG.Q.Value",
        "Precursor.Normalised",
        "PG.MaxLFQ",
    ],
    # presence_lib.R / infer_{family,species_strain}_presence.R
    "peptidotyping": [
        "PEP",
        "Q.Value",
        "Stripped.Sequence",
        "Decoy",
        "Protein.Ids",
    ],
    # build_genome_spectrum_mapping.py / build_taxon_spectrum_mapping.py
    "hapiid": [
        "Run",
        "Precursor.Id",
        "Protein.Ids",
        "Protein.Group",
        "Proteotypic",
    ],
}


class DiannEnvError(Exception):
    """Raised for an unusable DIA-NN configuration; message is user-facing."""


# ---------------------------------------------------------------------------
# Parse-time resolution
# ---------------------------------------------------------------------------
def _install_instructions(path, source):
    """The message a user sees when their diann_path does not resolve."""
    return f"""
DIA-NN was not found.

  configured path : {path}
  came from       : {source}

conduit-ascent does not distribute DIA-NN. Its licence permits one backup copy
and forbids sublicensing, so the binary cannot ship in a public image; you
download it yourself and accept its terms directly. DIA-NN is free for academic
use; commercial use requires a licence from the authors.

  download : {DIANN_DOWNLOAD_URL}
             (DIA-NN-<version>-Academia-Linux.zip)
  licence  : {DIANN_LICENCE_URL}

Then unzip it and point the workflow at it, in any one of these shapes:

  A. an extracted directory (recommended) — contains {DIANN_LINUX_BINARY}

       mkdir -p resources/diann
       unzip DIA-NN-2.5.0-Academia-Linux.zip -d resources/diann
       chmod +x resources/diann/diann-2.5.0/{DIANN_LINUX_BINARY}

     giving the default layout:

       resources/diann/diann-2.5.0/
       ├── {DIANN_LINUX_BINARY}
       ├── libtorch_cpu.so, libc10.so, libtimsdata.so, ...
       └── models/

     .NET 8 and libgomp come from the diann_runtime image; the directory is
     bind-mounted in. Keeping it under resources/ means no extra bind is
     needed — that path is inside the working directory and already visible
     to the container.

  B. a single executable — e.g. the self-contained DIA-NN.AppImage, or
     diann-linux on a host that already has .NET 8. Runs on the host, no
     container.

  C. a container image you built yourself — "docker://myrepo/diann:2.5.0" or
     "/path/to/diann.sif". The vendor ships a Dockerfile in the zip.

Select it with any of (first wins):

  snakemake --config diann_path=/path/to/diann-2.5.0 ...
  export {DIANN_PATH_ENV_VAR}=/path/to/diann-2.5.0
  diann_path: /path/to/diann-2.5.0      # in config/snakemake.yaml
""".rstrip()


def _resolve_configured_path(config, workflow):
    """Return (path, human-readable source) following CLI > env > config file."""
    cli = (getattr(workflow.config_settings, "overwrite_config", None) or {}).get("diann_path")
    if cli:
        return str(cli), "--config diann_path="

    env = os.environ.get(DIANN_PATH_ENV_VAR)
    if env:
        return env, f"${DIANN_PATH_ENV_VAR}"

    cfg = config.get("diann_path")
    if cfg:
        return str(cfg), "config file (diann_path)"

    return DEFAULT_DIANN_PATH, "built-in default"


def _classify(path):
    """Map a configured path to one of the three modes without touching config."""
    if path.startswith(_IMAGE_PREFIXES) or path.endswith(".sif"):
        return "image"
    if os.path.isdir(path):
        return "directory"
    if os.path.isfile(path):
        return "executable"
    return "missing"


def _is_within(child, parent):
    """True if `child` lives under `parent` (both resolved through symlinks)."""
    child = os.path.realpath(child)
    parent = os.path.realpath(parent)
    return child == parent or child.startswith(parent.rstrip(os.sep) + os.sep)


def _fingerprint(mode, cmd, container):
    """Identity of the DIA-NN that will run, for Snakemake's rerun triggers.

    The compatibility probe's output is a normal rule output, so without this
    Snakemake would happily keep a PASS recorded against a DIA-NN the user has
    since replaced. Feeding the fingerprint through the probe rule's params
    (a default rerun trigger) makes swapping the install re-run the check —
    including an in-place upgrade, which leaves the path unchanged but not the
    binary's size or mtime.
    """
    if mode in ("directory", "executable") and cmd and os.path.isfile(cmd):
        stat = os.stat(cmd)
        return f"{cmd}:{stat.st_size}:{int(stat.st_mtime)}"
    return str(container or cmd or "")


def resolve_diann_environment(config, workflow, strict=True):
    """Resolve the DIA-NN installation and wire it into `config`.

    Sets, for the rules to consume:

    ``config["diann_cmd"]``
        The command that invokes DIA-NN, used in shell blocks as
        ``{config[diann_cmd]}``. ``diann`` in mode C (on the image's PATH),
        an absolute path to the binary in modes A and B.
    ``config["containers"]["diann"]``
        The container the DIA-NN rules run in: the ``diann_runtime`` image in
        mode A, the user's own image in mode C, ``None`` in mode B (Snakemake
        reads ``container: None`` as "run this rule on the host").
    ``config["diann"]``
        Provenance: mode, resolved path, where the setting came from, and the
        bind added. Lands in manifest.json via the resolved-config snapshot.

    With ``strict=False`` (dry runs) a broken installation is reported as a
    warning instead of an exception, so `tests/run_dry_runs.sh` and CI can
    still build the DAG on a machine that has no DIA-NN.

    Returns the ``config["diann"]`` dict.
    """
    path, source = _resolve_configured_path(config, workflow)
    mode = _classify(path)

    container = None
    cmd = None
    bind = None
    problem = None

    if mode == "image":
        # Mode C — unchanged behaviour: DIA-NN is on the image's PATH.
        container = path
        cmd = "diann"

    elif mode == "directory":
        # Mode A — the extracted vendor archive, run inside diann_runtime.
        binary = os.path.join(os.path.abspath(path), DIANN_LINUX_BINARY)
        if not os.path.isfile(binary):
            problem = (
                f"{path} is a directory but contains no {DIANN_LINUX_BINARY}.\n"
                f"Point diann_path at the directory the zip extracts to "
                f"(the one holding {DIANN_LINUX_BINARY}), not its parent."
            )
        elif not os.access(binary, os.X_OK):
            problem = (
                f"{binary} is not executable.\n"
                f"The vendor zip does not always preserve the execute bit. Fix it with:\n"
                f"    chmod +x {binary}"
            )
        else:
            runtime = (config.get("containers") or {}).get("diann_runtime")
            if not runtime:
                problem = (
                    "config['containers']['diann_runtime'] is not set, so there is no "
                    "runtime image to bind the DIA-NN directory into."
                )
            else:
                container = runtime
                cmd = binary
                # Apptainer binds are global (--singularity-args), not per-rule.
                # Anything under the working directory is auto-bound already;
                # an out-of-tree install needs an explicit bind appended to the
                # deployment settings at parse time.
                install_dir = os.path.realpath(path)
                if not _is_within(install_dir, os.getcwd()):
                    bind = install_dir

    elif mode == "executable":
        # Mode B — run on the host; the file carries its own runtime.
        if not os.access(path, os.X_OK):
            problem = (
                f"{path} is not executable.\n"
                f"    chmod +x {path}"
            )
        else:
            container = None  # `container: None` → Snakemake runs it on the host
            cmd = os.path.abspath(path)

    else:
        problem = None  # missing: the full install message below says it all

    if problem is not None or mode == "missing":
        message = _install_instructions(path, source)
        if problem:
            message = f"{problem}\n{message}"
        if strict:
            raise DiannEnvError(message)
        print(
            "WARNING: DIA-NN preflight failed; continuing because this is a dry run.\n"
            + message,
            file=sys.stderr,
        )
        # Keep the DAG buildable: a placeholder command that cannot silently
        # succeed if a rule were somehow executed.
        cmd = "diann"
        container = (config.get("containers") or {}).get("diann_runtime")

    config.setdefault("containers", {})["diann"] = container
    config["diann_cmd"] = cmd

    if bind:
        args = workflow.deployment_settings.apptainer_args or ""
        if bind not in args:
            workflow.deployment_settings.apptainer_args = (
                f"{args} --bind {bind}".strip()
            )

    info = {
        "mode": mode,
        "configured_path": path,
        "source": source,
        "resolved_command": cmd,
        "container": container,
        "bind": bind,
        "fingerprint": _fingerprint(mode, cmd, container),
        "usable": problem is None and mode != "missing",
    }
    config["diann"] = info

    if info["usable"]:
        detail = f"mode {mode} via {container}" if container else f"mode {mode} on the host"
        print(f"DIA-NN: {cmd} ({detail}; from {source})", file=sys.stderr)

    return info


# ---------------------------------------------------------------------------
# Run-time compatibility probe
# ---------------------------------------------------------------------------
def build_compat_probe_argv(diann_cmd, tmpdir):
    """Argv that exercises every flag the workflow uses, processing no files.

    DIA-NN prints its banner (version included) and then, for each flag it does
    not know, ``WARNING: unrecognised option [--flag]``. With no ``--f``/``--dir``
    it reports "0 files will be processed" and exits in a couple of seconds, so
    this is a cheap, complete check of the CLI half of the contract.
    """
    argv = [diann_cmd]
    for flag, value in DIANN_REQUIRED_FLAGS:
        argv.append(flag)
        if value is not None:
            argv.append(value)
    argv.append(COMPAT_SENTINEL_FLAG)
    argv += ["--temp", tmpdir, "--out", os.path.join(tmpdir, "probe_report")]
    return argv


_VERSION_RE = re.compile(r"DIA-NN\s+(\d+(?:\.\d+){0,2})")
# DIA-NN 2.x writes "unrecognised"; accept the US spelling too so a future
# rewording of the message does not quietly turn every check into a pass.
_UNRECOGNISED_RE = re.compile(r"unrecogni[sz]ed option \[([^\]]+)\]")


def _parse_version(text):
    match = _VERSION_RE.search(text)
    return match.group(1) if match else None


def _version_tuple(version):
    return tuple(int(part) for part in version.split("."))


def evaluate_compat_probe(output, min_version=DIANN_MIN_VERSION,
                          tested_versions=DIANN_TESTED_VERSIONS):
    """Turn probe output into a verdict dict. Pure function — easy to unit-test.

    Returns ``{"ok", "version", "unsupported_flags", "errors", "warnings"}``.
    """
    errors = []
    warnings = []

    version = _parse_version(output)
    if version is None:
        errors.append(
            "Could not read a version banner from DIA-NN's output. The command "
            "ran but does not look like DIA-NN."
        )
    else:
        try:
            parsed = _version_tuple(version)
        except ValueError:
            parsed = None
            warnings.append(f"Could not compare version string {version!r} numerically.")
        if parsed is not None:
            padded = parsed + (0,) * (len(min_version) - len(parsed))
            if padded < min_version:
                errors.append(
                    f"DIA-NN {version} is older than the minimum this workflow "
                    f"supports ({'.'.join(str(p) for p in min_version)}). Earlier "
                    f"builds lack flags the pipeline depends on."
                )
            elif version not in tested_versions:
                warnings.append(
                    f"DIA-NN {version} has not been tested with this workflow "
                    f"(tested: {', '.join(tested_versions)}). Proceeding because it "
                    f"accepts every flag the pipeline uses."
                )

    unrecognised = set(_UNRECOGNISED_RE.findall(output))

    # Self-test: if the deliberately bogus flag was not reported, DIA-NN is not
    # telling us about unknown flags at all and a clean result means nothing.
    if COMPAT_SENTINEL_FLAG not in unrecognised:
        errors.append(
            "The compatibility probe could not verify itself: DIA-NN did not "
            f"report the deliberately invalid flag {COMPAT_SENTINEL_FLAG} as "
            "unrecognised. Its diagnostics have changed shape, so flag support "
            "cannot be confirmed. Treat this DIA-NN build as unverified."
        )
    unrecognised.discard(COMPAT_SENTINEL_FLAG)

    required = {flag for flag, _ in DIANN_REQUIRED_FLAGS}
    unsupported = sorted(unrecognised & required)
    if unsupported:
        errors.append(
            "This DIA-NN build does not support flags the workflow passes: "
            + ", ".join(unsupported)
            + ".\nDIA-NN ignores unknown flags rather than failing, so a run would "
            "silently do the wrong thing (e.g. without --pre-search/--pre-filter "
            "every InfiniDIA search degrades to a plain library search)."
        )

    return {
        "ok": not errors,
        "version": version,
        "unsupported_flags": unsupported,
        "errors": errors,
        "warnings": warnings,
    }


def render_compat_report(verdict, diann_info):
    """Human-readable summary written next to the machine-readable JSON."""
    lines = [
        "DIA-NN compatibility check",
        "==========================",
        f"mode       : {diann_info.get('mode')}",
        f"path       : {diann_info.get('configured_path')}",
        f"command    : {diann_info.get('resolved_command')}",
        f"container  : {diann_info.get('container') or '(host)'}",
        f"version    : {verdict.get('version') or 'unknown'}",
        f"result     : {'PASS' if verdict.get('ok') else 'FAIL'}",
    ]
    for warning in verdict.get("warnings", []):
        lines += ["", "WARNING: " + warning]
    for error in verdict.get("errors", []):
        lines += ["", "ERROR: " + error]
    return "\n".join(lines) + "\n"


def write_compat_artifacts(json_path, verdict, diann_info):
    os.makedirs(os.path.dirname(json_path) or ".", exist_ok=True)
    payload = dict(verdict)
    payload["diann"] = diann_info
    payload["probe_version"] = COMPAT_PROBE_VERSION
    with open(json_path, "w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
    txt_path = os.path.splitext(json_path)[0] + ".txt"
    with open(txt_path, "w") as handle:
        handle.write(render_compat_report(verdict, diann_info))
    return txt_path
