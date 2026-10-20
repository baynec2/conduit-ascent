"""Helpers for the staged DIA-NN search rules.

The DIA-NN search modules are split into three stages: build empirical library,
per-raw search (MBR off, conventional engine), combine via --use-quant. These
helpers handle (a) enumerating MS files for the per-raw fan-out and (b) staging
the per-raw {sample}.quant outputs where the combine stage's DIA-NN looks for
them.

The combine stage passes --quant-ori-names, under which DIA-NN looks for each
run's .quant as <temp>/<raw basename without extension>.quant — the same name
the per-raw rules give their output. Without it DIA-NN mangles the raw path *as
passed* (relative, '/' and '.' -> '_', '-' kept), which is easy to mispredict:
a mismatch is not an error, DIA-NN just silently re-searches the run.

Verified against DIA-NN 2.5.0 on 2026-10-19: a .quant written without the flag
is byte-identical to one written with it, and combine reuses it (no "Main
search") whether --dir is input/ms_files or a directory of symlinks to it.
"""
import glob
import os

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


def stage3_symlink_commands(combine_temp_dir, raws, quant_files_dir):
    """`&&`-joined `ln -sf` chain: link each canonical {sample}.quant into
    combine_temp_dir for --use-quant --quant-ori-names."""
    cmds = []
    for raw in raws:
        sample = os.path.splitext(os.path.basename(raw))[0]
        src = os.path.abspath(os.path.join(quant_files_dir, sample + _QUANT_SUFFIX))
        dst = os.path.join(combine_temp_dir, sample + _QUANT_SUFFIX)
        cmds.append(f"ln -sf {src} {dst}")
    return " && ".join(cmds) if cmds else "true"
