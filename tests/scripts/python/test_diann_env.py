"""Unit tests for the DIA-NN environment resolution and compatibility probe.

These cover the pure logic in modules/_shared/diann_env.py — mode detection,
precedence, and probe-output parsing — without needing DIA-NN, a container, or
snakemake. The parsing tests are the important ones: DIA-NN never exits
non-zero for an unknown flag, so the *only* signal that a user's build cannot
honour a flag is the wording of a warning line, and that wording is what these
assertions pin down.
"""

import os
import stat
import types

import pytest

import diann_env


# ---------------------------------------------------------------------------
# Mode detection
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "path,expected",
    [
        ("docker://baynec2/diann:abc123", "image"),
        ("docker://myrepo/diann:2.5.0", "image"),
        ("/opt/images/diann.sif", "image"),
        ("library://someone/diann", "image"),
    ],
)
def test_classify_image_uris(path, expected):
    assert diann_env._classify(path) == expected


def test_classify_directory_and_executable(tmp_path):
    install = tmp_path / "diann-2.5.0"
    install.mkdir()
    binary = install / diann_env.DIANN_LINUX_BINARY
    binary.write_text("#!/bin/sh\n")

    assert diann_env._classify(str(install)) == "directory"
    assert diann_env._classify(str(binary)) == "executable"
    assert diann_env._classify(str(tmp_path / "nope")) == "missing"


# ---------------------------------------------------------------------------
# Resolution and precedence
# ---------------------------------------------------------------------------
def _fake_workflow(cli_config=None):
    """Minimal stand-in for snakemake's Workflow: the two attrs we touch."""
    return types.SimpleNamespace(
        config_settings=types.SimpleNamespace(overwrite_config=cli_config or {}),
        deployment_settings=types.SimpleNamespace(apptainer_args=""),
    )


def _make_install(tmp_path, name="diann-2.5.0", executable=True):
    install = tmp_path / name
    install.mkdir()
    binary = install / diann_env.DIANN_LINUX_BINARY
    binary.write_text("#!/bin/sh\necho DIA-NN\n")
    if executable:
        binary.chmod(binary.stat().st_mode | stat.S_IXUSR)
    return install, binary


def test_directory_mode_uses_runtime_image_and_binary(tmp_path):
    install, binary = _make_install(tmp_path)
    config = {
        "diann_path": str(install),
        "containers": {"diann_runtime": "docker://baynec2/diann_runtime:abc123"},
    }
    info = diann_env.resolve_diann_environment(config, _fake_workflow())

    assert info["mode"] == "directory"
    assert info["usable"]
    assert config["diann_cmd"] == str(binary)
    assert config["containers"]["diann"] == "docker://baynec2/diann_runtime:abc123"


def test_directory_mode_binds_out_of_tree_install(tmp_path):
    install, _ = _make_install(tmp_path)
    workflow = _fake_workflow()
    config = {
        "diann_path": str(install),
        "containers": {"diann_runtime": "docker://baynec2/diann_runtime:abc123"},
    }
    diann_env.resolve_diann_environment(config, workflow)

    # tmp_path is outside the working directory, so a bind must be appended.
    assert f"--bind {os.path.realpath(install)}" in workflow.deployment_settings.apptainer_args


def test_executable_mode_runs_on_host(tmp_path):
    _, binary = _make_install(tmp_path)
    config = {"diann_path": str(binary), "containers": {}}
    info = diann_env.resolve_diann_environment(config, _fake_workflow())

    assert info["mode"] == "executable"
    # `container: None` is how a rule opts out of the workflow's container.
    assert config["containers"]["diann"] is None
    assert config["diann_cmd"] == str(binary)


def test_image_mode_is_unchanged_behaviour():
    config = {"diann_path": "docker://myrepo/diann:2.5.0", "containers": {}}
    info = diann_env.resolve_diann_environment(config, _fake_workflow())

    assert info["mode"] == "image"
    assert config["containers"]["diann"] == "docker://myrepo/diann:2.5.0"
    assert config["diann_cmd"] == "diann"


def test_non_executable_binary_is_rejected_with_chmod_hint(tmp_path):
    install, binary = _make_install(tmp_path, executable=False)
    config = {
        "diann_path": str(install),
        "containers": {"diann_runtime": "docker://baynec2/diann_runtime:abc123"},
    }
    with pytest.raises(diann_env.DiannEnvError) as excinfo:
        diann_env.resolve_diann_environment(config, _fake_workflow())
    assert "chmod +x" in str(excinfo.value)


def test_missing_install_is_a_warning_not_an_error_when_not_strict(tmp_path, capsys):
    config = {"diann_path": str(tmp_path / "absent"), "containers": {"diann_runtime": "img"}}
    info = diann_env.resolve_diann_environment(config, _fake_workflow(), strict=False)
    assert not info["usable"]
    assert "download" in capsys.readouterr().err


def test_cli_config_beats_env_var(tmp_path, monkeypatch):
    monkeypatch.setenv(diann_env.DIANN_PATH_ENV_VAR, "docker://from-env/diann:1")
    config = {"diann_path": "docker://from-file/diann:1", "containers": {}}
    workflow = _fake_workflow({"diann_path": "docker://from-cli/diann:1"})
    info = diann_env.resolve_diann_environment(config, workflow)
    assert info["configured_path"] == "docker://from-cli/diann:1"


def test_env_var_beats_config_file(tmp_path, monkeypatch):
    monkeypatch.setenv(diann_env.DIANN_PATH_ENV_VAR, "docker://from-env/diann:1")
    config = {"diann_path": "docker://from-file/diann:1", "containers": {}}
    info = diann_env.resolve_diann_environment(config, _fake_workflow())
    assert info["configured_path"] == "docker://from-env/diann:1"


# ---------------------------------------------------------------------------
# Probe construction and parsing
# ---------------------------------------------------------------------------
def test_probe_argv_carries_every_required_flag_and_the_sentinel():
    argv = diann_env.build_compat_probe_argv("/opt/diann/diann-linux", "/tmp/probe")
    assert argv[0] == "/opt/diann/diann-linux"
    for flag, _ in diann_env.DIANN_REQUIRED_FLAGS:
        assert flag in argv
    assert diann_env.COMPAT_SENTINEL_FLAG in argv


def _probe_output(version="2.5.0", unrecognised=(), include_sentinel=True):
    lines = [
        "",
        f"DIA-NN {version} Academia  (Data-Independent Acquisition by Neural Networks)",
        "Compiled on Sep 26 2025 02:56:25",
    ]
    flags = list(unrecognised)
    if include_sentinel:
        flags.append(diann_env.COMPAT_SENTINEL_FLAG)
    lines += [f"WARNING: unrecognised option [{flag}]" for flag in flags]
    lines += ["0 files will be processed", "Finished"]
    return "\n".join(lines)


def test_supported_version_with_all_flags_passes():
    verdict = diann_env.evaluate_compat_probe(_probe_output("2.5.0"))
    assert verdict["ok"]
    assert verdict["version"] == "2.5.0"
    assert verdict["unsupported_flags"] == []


def test_unsupported_infinidia_flags_fail():
    # The real DIA-NN 2.1.0 failure: --pre-search/--pre-filter silently ignored,
    # which would downgrade every InfiniDIA search to a plain library search.
    verdict = diann_env.evaluate_compat_probe(
        _probe_output("2.1.0", ["--pre-search", "--pre-filter"])
    )
    assert not verdict["ok"]
    assert verdict["unsupported_flags"] == ["--pre-filter", "--pre-search"]
    assert any("older than the minimum" in e for e in verdict["errors"])


def test_untested_but_new_enough_version_only_warns():
    verdict = diann_env.evaluate_compat_probe(_probe_output("2.9.0"))
    assert verdict["ok"]
    assert any("has not been tested" in w for w in verdict["warnings"])


def test_missing_sentinel_fails_the_self_test():
    # A DIA-NN that reports nothing as unrecognised gives us no evidence either
    # way, so a clean probe must not be read as a pass.
    verdict = diann_env.evaluate_compat_probe(
        _probe_output("2.5.0", include_sentinel=False)
    )
    assert not verdict["ok"]
    assert any("could not verify itself" in e for e in verdict["errors"])


def test_unparseable_output_fails():
    verdict = diann_env.evaluate_compat_probe("bash: diann-linux: No such file or directory")
    assert not verdict["ok"]
    assert verdict["version"] is None


def test_us_spelling_of_unrecognized_is_also_matched():
    output = _probe_output("2.5.0", include_sentinel=False).replace(
        "0 files will be processed",
        f"WARNING: unrecognized option [{diann_env.COMPAT_SENTINEL_FLAG}]\n"
        "WARNING: unrecognized option [--pre-search]\n0 files will be processed",
    )
    verdict = diann_env.evaluate_compat_probe(output)
    assert verdict["unsupported_flags"] == ["--pre-search"]


def test_required_report_columns_cover_every_consumer_profile():
    assert set(diann_env.DIANN_REQUIRED_REPORT_COLUMNS) == {"main", "peptidotyping", "hapiid"}
    for columns in diann_env.DIANN_REQUIRED_REPORT_COLUMNS.values():
        assert columns, "a profile with no required columns checks nothing"
        assert len(columns) == len(set(columns))
