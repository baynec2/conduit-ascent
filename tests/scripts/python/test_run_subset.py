"""Unit tests for MS run-subset selection (modules/_shared/run_subset.py)."""

import os

import pandas as pd
import pytest

import run_subset


def _annotation():
    rows = []
    for i in range(12):
        rows.append({
            "file": f"s{i:02d}",
            "treatment": ["untreated", "metro", "amp"][i % 3],
            "day": [0, 21][i // 6],
        })
    return pd.DataFrame(rows)


SAMPLES = [f"s{i:02d}" for i in range(12)]


# ---------------------------------------------------------------------------
# Spec validation
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("spec", [None, False, "null"])
def test_null_spec_means_no_subset(spec):
    assert run_subset.parse_subset_spec(spec, "k") is None


@pytest.mark.parametrize(
    "spec,match",
    [
        ({}, "exactly one"),
        ({"fraction": 0.5, "samples": ["s00"]}, "exactly one"),
        ({"fraction": 0}, r"\(0, 1\]"),
        ({"fraction": 1.5}, r"\(0, 1\]"),
        ({"fraction": "0.5"}, "number"),
        ({"fraction": 0.5, "seed": "x"}, "integer"),
        ({"fraction": 0.5, "nope": 1}, "unknown field"),
        ({"samples": []}, "non-empty"),
        ({"samples": "s00"}, "non-empty"),
        ({"samples": ["s00"], "stratify_by": ["day"]}, "only applies"),
        (0.25, "mapping"),
    ],
)
def test_invalid_specs_raise(spec, match):
    with pytest.raises(ValueError, match=match):
        run_subset.parse_subset_spec(spec, "k")


def test_stratify_by_string_is_normalised_to_list():
    spec = run_subset.parse_subset_spec({"fraction": 0.5, "stratify_by": "day"}, "k")
    assert spec == {"fraction": 0.5, "stratify_by": ["day"], "seed": 1}


# ---------------------------------------------------------------------------
# Selection
# ---------------------------------------------------------------------------
def _select(**spec):
    return run_subset.select_run_subset(
        _annotation(), SAMPLES, run_subset.parse_subset_spec(spec, "k"), "k"
    )


def test_unstratified_fraction_takes_ceil_of_total():
    picked = _select(fraction=0.25)
    assert len(picked) == 3
    assert {label for _, label in picked} == {"all"}


def test_stratified_takes_ceil_per_stratum_and_covers_every_stratum():
    picked = _select(fraction=0.25, stratify_by=["treatment", "day"])
    # 6 strata of 2 runs each: ceil(0.25 * 2) = 1 per stratum.
    assert len(picked) == 6
    assert len({label for _, label in picked}) == 6


def test_small_strata_do_not_inflate_the_total():
    # CB024 shape: 24 strata of 5 runs. Per-stratum ceil would take 2 each (48).
    df = pd.DataFrame({"file": [f"r{i:03d}" for i in range(120)],
                       "grp": [i // 5 for i in range(120)]})
    spec = run_subset.parse_subset_spec({"fraction": 0.25, "stratify_by": ["grp"]}, "k")
    picked = run_subset.select_run_subset(df, df["file"].tolist(), spec, "k")
    assert len(picked) == 30
    per = pd.Series([label for _, label in picked]).value_counts()
    assert len(per) == 24 and set(per) == {1, 2}


def test_allocation_is_proportional_to_stratum_size():
    quotas = run_subset._allocate({"big": 30, "small": 10}, 0.25, 1)
    assert quotas == {"big": 7, "small": 3}


def test_every_stratum_keeps_at_least_one_run():
    picked = _select(fraction=0.01, stratify_by=["treatment"])
    assert sorted(label for _, label in picked) == [
        "treatment=amp", "treatment=metro", "treatment=untreated"
    ]


def test_fraction_one_selects_everything():
    assert [s for s, _ in _select(fraction=1, stratify_by=["day"])] == SAMPLES


def test_selection_is_deterministic_and_seed_dependent():
    a = _select(fraction=0.5, stratify_by=["day"], seed=1)
    assert a == _select(fraction=0.5, stratify_by=["day"], seed=1)
    draws = {tuple(_select(fraction=0.5, stratify_by=["day"], seed=s)) for s in range(10)}
    assert len(draws) > 1


def test_adding_a_stratum_leaves_other_strata_draws_unchanged():
    before = _select(fraction=0.5, stratify_by=["treatment"])
    df = pd.concat([_annotation(), pd.DataFrame([
        {"file": "x0", "treatment": "vanc", "day": 0},
        {"file": "x1", "treatment": "vanc", "day": 0},
    ])])
    spec = run_subset.parse_subset_spec({"fraction": 0.5, "stratify_by": ["treatment"]}, "k")
    after = run_subset.select_run_subset(df, SAMPLES + ["x0", "x1"], spec, "k")
    assert [p for p in after if not p[1].endswith("vanc")] == before


def test_missing_stratify_values_form_their_own_stratum():
    df = _annotation()
    df.loc[0, "treatment"] = None
    picked = run_subset.select_run_subset(
        df, SAMPLES, run_subset.parse_subset_spec({"fraction": 0.1, "stratify_by": ["treatment"]}, "k"), "k"
    )
    assert ("s00", "treatment=NA") in picked


def test_unknown_stratify_column_raises():
    with pytest.raises(ValueError, match="not in sample_annotation"):
        _select(fraction=0.5, stratify_by=["cage"])


def test_annotation_rows_without_ms_files_are_ignored():
    df = pd.concat([_annotation(), pd.DataFrame([{"file": "ghost", "treatment": "amp", "day": 0}])])
    spec = run_subset.parse_subset_spec({"fraction": 1}, "k")
    assert "ghost" not in [s for s, _ in run_subset.select_run_subset(df, SAMPLES, spec, "k")]


def test_explicit_samples():
    assert _select(samples=["s03", "s01", "s01"]) == [("s01", "explicit"), ("s03", "explicit")]


def test_explicit_samples_must_exist():
    with pytest.raises(ValueError, match="no MS file: zz"):
        _select(samples=["s01", "zz"])


# ---------------------------------------------------------------------------
# Resolution into rule inputs
# ---------------------------------------------------------------------------
def test_no_subset_renders_original_inputs(tmp_path):
    raws = [str(tmp_path / f"exp/input/ms_files/{s}.raw") for s in SAMPLES]
    config = {}
    config["_run_subsets"] = run_subset.resolve_run_subsets(config, _annotation(), SAMPLES, raws)
    assert config["_run_subsets"] == {"diann_empirical_lib": None, "search_space": None}
    ms_dir, out_raws, out_samples = run_subset.run_subset_inputs(
        config, "search_space", str(tmp_path / "exp"), "run", raws, SAMPLES
    )
    # Must be byte-identical to what the rules used before subsets existed.
    assert ms_dir == os.path.join(str(tmp_path / "exp"), "input/ms_files")
    assert out_raws is raws and out_samples is SAMPLES


def test_subset_points_rules_at_staged_dir(tmp_path):
    raws = [str(tmp_path / f"ms/{s}.mzML") for s in SAMPLES]
    config = {"search_space_run_subset": {"samples": ["s02", "s05"]}}
    config["_run_subsets"] = run_subset.resolve_run_subsets(config, _annotation(), SAMPLES, raws)
    assert config["_run_subsets"]["diann_empirical_lib"] is None
    assert config["_run_subsets"]["search_space"]["raws"] == [raws[2], raws[5]]
    ms_dir, out_raws, out_samples = run_subset.run_subset_inputs(
        config, "search_space", "exp", "runs/r1", raws, SAMPLES
    )
    assert ms_dir == "runs/r1/ms_subsets/search_space"
    assert out_raws == [f"{ms_dir}/s02.mzML", f"{ms_dir}/s05.mzML"]
    assert out_samples == ["s02", "s05"]


# ---------------------------------------------------------------------------
# Staging
# ---------------------------------------------------------------------------
def test_stage_runs_hard_links_and_resolves_symlinked_inputs(tmp_path):
    real = tmp_path / "data" / "a.raw"
    real.parent.mkdir()
    real.write_bytes(b"raw")
    ms_files = tmp_path / "ms_files"
    ms_files.mkdir()
    (ms_files / "a.raw").symlink_to(real)  # e.g. integration test ms_files
    out = tmp_path / "subset"
    run_subset.stage_runs([str(ms_files / "a.raw")], str(out))
    staged = out / "a.raw"
    # A hard link to the real file, not a symlink: DIA-NN cannot read .raw via symlinks.
    assert not staged.is_symlink()
    assert os.path.samefile(staged, real)


def _cross_device(*_args, **_kwargs):
    import errno
    raise OSError(errno.EXDEV, "Invalid cross-device link")


def test_stage_runs_raw_across_filesystems_raises(tmp_path, monkeypatch):
    raw = tmp_path / "a.raw"
    raw.write_bytes(b"raw")
    monkeypatch.setattr(run_subset.os, "link", _cross_device)
    with pytest.raises(OSError, match="different filesystems"):
        run_subset.stage_runs([str(raw)], str(tmp_path / "subset"))


def test_stage_runs_mzml_across_filesystems_falls_back_to_symlink(tmp_path, monkeypatch):
    mzml = tmp_path / "a.mzML"
    mzml.write_bytes(b"mzml")
    monkeypatch.setattr(run_subset.os, "link", _cross_device)
    run_subset.stage_runs([str(mzml)], str(tmp_path / "subset"))
    assert (tmp_path / "subset" / "a.mzML").is_symlink()
