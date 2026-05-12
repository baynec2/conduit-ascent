"""Helpers for the staged DIA-NN search rules.

The DIA-NN search modules are split into three stages: build empirical library,
per-raw search (MBR off, conventional engine), combine via --use-quant. These
helpers handle (a) enumerating MS files for the per-raw fan-out and (b) bridging
between Snakemake's canonical {sample}.quant output paths and DIA-NN's
mangled-path .quant filenames.

DIA-NN names .quant files by mangling the absolute path of the raw file:
non-alphanumeric chars become '_', then '.quant' is appended.

Verified against baynec2/diann:f5f961d (DIA-NN 2.5.0) on 2026-05-08:
    /work/raw/sample1.mzML  ->  _work_raw_sample1_mzML.quant
"""
import glob
import os
import re

_QUANT_SUFFIX = ".quant"


def list_raw_files(experiment_dir):
    return sorted(
        glob.glob(os.path.join(experiment_dir, "input/ms_files/*.raw"))
        + glob.glob(os.path.join(experiment_dir, "input/ms_files/*.mzML"))
    )


def list_samples(experiment_dir):
    return sorted({
        os.path.splitext(os.path.basename(p))[0]
        for p in list_raw_files(experiment_dir)
    })


def raw_path_for_sample(experiment_dir, sample):
    for p in list_raw_files(experiment_dir):
        if os.path.splitext(os.path.basename(p))[0] == sample:
            return p
    raise KeyError(
        f"no raw file for sample {sample!r} in {experiment_dir}/input/ms_files"
    )


def diann_quant_filename(abs_raw_path):
    return re.sub(r"[^A-Za-z0-9]", "_", abs_raw_path) + _QUANT_SUFFIX


def stage3_symlink_commands(combine_temp_dir, raws, quant_files_dir):
    """`&&`-joined `ln -sf` chain: stage canonical {sample}.quant files into
    DIA-NN's expected mangled-path layout in combine_temp_dir for --use-quant."""
    cmds = []
    for raw in raws:
        sample = os.path.splitext(os.path.basename(raw))[0]
        src = os.path.abspath(os.path.join(quant_files_dir, sample + _QUANT_SUFFIX))
        dst = os.path.join(combine_temp_dir, diann_quant_filename(os.path.abspath(raw)))
        cmds.append(f"ln -sf {src} {dst}")
    return " && ".join(cmds) if cmds else "true"
