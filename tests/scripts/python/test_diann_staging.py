"""Unit tests for the staged DIA-NN helpers (modules/_shared/diann_staging.py)."""

import os

import diann_staging


def test_stage3_links_basename_quants_for_quant_ori_names(tmp_path):
    # Combine runs with --quant-ori-names, so DIA-NN looks for <basename>.quant
    # whatever directory the raw sits in (input/ms_files or an ms_subsets dir).
    raws = [
        "experiments/e/input/ms_files/CB_Stool-P1A1.raw",
        "experiments/e/runs/r/ms_subsets/search_space/sample2.mzML",
    ]
    cmds = diann_staging.stage3_symlink_commands("tmp/combine", raws, str(tmp_path))
    assert cmds.split(" && ") == [
        f"ln -sf {os.path.join(tmp_path, 'CB_Stool-P1A1.quant')} tmp/combine/CB_Stool-P1A1.quant",
        f"ln -sf {os.path.join(tmp_path, 'sample2.quant')} tmp/combine/sample2.quant",
    ]


def test_stage3_with_no_raws_is_a_noop():
    assert diann_staging.stage3_symlink_commands("tmp", [], "q") == "true"
