"""Select a subset of MS runs for DIA-NN steps that don't need every run.

Two kinds of step only use the runs to *learn* something, not to quantify them:

  * ``diann_empirical_lib`` — stage 1 of the standard main search builds the
    empirical spectral library. Stages 2 and 3 still search and quantify every
    run against that library, so no sample is dropped from the final matrix.
  * ``search_space`` — the DIA-NN passes of the organism-selection methods
    (unipept_peptidotyping, unipept_hapiid, hapiid, genome_peptidotyping). Their
    only output is the set of organisms that goes into the database.

Both steps' memory grows with the number of runs (DIA-NN keeps per-run state
for cross-run quantification), which is what OOMs large cohorts. Each step has
its own config key because they tolerate dropped runs differently: a peptide
missing from the library costs a few IDs, while an organism missed during
selection never enters the database at all.

Config shape (either key; ``null`` = use every run, the default)::

    diann_empirical_lib_subset:
      fraction: 0.25                            # 0 < fraction <= 1
      stratify_by: [antibiotic_treatment, day]  # sample_annotation columns
      seed: 1
    # or pick the runs explicitly:
    search_space_run_subset:
      samples: [CB_Stool_P1A1, CB_Stool_P1B1]

With ``fraction``, runs are grouped by the ``stratify_by`` columns (one group if
omitted) and about ``fraction`` of the runs are drawn, spread across groups in
proportion to their size with at least one per group, so no condition drops out
(see ``_allocate``). Draws are seeded, so the same config always selects the
same runs — Snakemake re-parses on every invocation,
and an unseeded draw would invalidate the library each time.
"""
import math
import os
import random

SUBSET_KEYS = {
    "diann_empirical_lib": "diann_empirical_lib_subset",
    "search_space": "search_space_run_subset",
}
_SPEC_FIELDS = {"fraction", "stratify_by", "seed", "samples"}
_DEFAULT_SEED = 1
# Annotation cells that are empty (NaN) are grouped under this label.
_NA_LABEL = "NA"


def parse_subset_spec(spec, key):
    """Validate a subset spec from config; return a normalised dict or None."""
    if spec in (None, False, "null", "None"):
        return None
    if not isinstance(spec, dict):
        raise ValueError(f"{key} must be a mapping or null (got {spec!r})")
    unknown = set(spec) - _SPEC_FIELDS
    if unknown:
        raise ValueError(
            f"{key}: unknown field(s) {sorted(unknown)}; allowed: {sorted(_SPEC_FIELDS)}"
        )
    fraction = spec.get("fraction")
    samples = spec.get("samples")
    stratify_by = spec.get("stratify_by")
    if (fraction is None) == (samples is None):
        raise ValueError(f"{key}: set exactly one of 'fraction' or 'samples'")
    if samples is not None:
        if stratify_by is not None:
            raise ValueError(f"{key}: 'stratify_by' only applies with 'fraction'")
        if isinstance(samples, str) or not samples:
            raise ValueError(f"{key}: 'samples' must be a non-empty list")
        return {"samples": [str(s) for s in samples]}
    if isinstance(fraction, bool) or not isinstance(fraction, (int, float)):
        raise ValueError(f"{key}: 'fraction' must be a number (got {fraction!r})")
    if not 0 < fraction <= 1:
        raise ValueError(f"{key}: 'fraction' must be in (0, 1] (got {fraction})")
    if stratify_by is None:
        stratify_by = []
    elif isinstance(stratify_by, str):
        stratify_by = [stratify_by]
    seed = spec.get("seed", _DEFAULT_SEED)
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise ValueError(f"{key}: 'seed' must be an integer (got {seed!r})")
    return {"fraction": float(fraction), "stratify_by": list(stratify_by), "seed": seed}


def select_run_subset(sample_df, samples, spec, key):
    """Return ``[(sample, stratum_label), ...]`` for the selected runs, sorted by sample.

    ``sample_df`` is the sample annotation (needs a ``file`` column); ``samples``
    are the sample names that have an MS file. ``spec`` is a parsed spec.
    """
    available = sorted(set(samples))
    if "samples" in spec:
        missing = sorted(set(spec["samples"]) - set(available))
        if missing:
            raise ValueError(f"{key}: samples with no MS file: {', '.join(missing)}")
        return [(s, "explicit") for s in sorted(set(spec["samples"]))]

    stratify_by = spec["stratify_by"]
    absent = [c for c in stratify_by if c not in sample_df.columns]
    if absent:
        raise ValueError(
            f"{key}: stratify_by column(s) not in sample_annotation: {', '.join(absent)}"
        )
    df = sample_df.assign(file=sample_df["file"].astype(str))
    df = df[df["file"].isin(available)].drop_duplicates("file")

    def _label(row):
        if not stratify_by:
            return "all"
        return "|".join(
            f"{c}={_NA_LABEL if _is_na(row[c]) else str(row[c]).strip()}"
            for c in stratify_by
        )

    strata = {}
    for _, row in df.iterrows():
        strata.setdefault(_label(row), []).append(row["file"])

    quotas = _allocate(
        {label: len(members) for label, members in strata.items()},
        spec["fraction"], spec["seed"],
    )
    selected = []
    for label, members in strata.items():
        # Shuffle with a per-stratum seed and take a prefix: adding a stratum or
        # changing one quota leaves every other stratum's draw unchanged.
        order = sorted(members)
        random.Random(f"{spec['seed']}|{label}").shuffle(order)
        selected.extend((s, label) for s in order[:quotas[label]])
    return sorted(selected)


def _allocate(sizes, fraction, seed):
    """Per-stratum run counts totalling ~``ceil(fraction * N)``, at least 1 each.

    Each stratum gets ``floor(fraction * size)`` (min 1); leftover slots up to the
    target go to the strata with the largest rounding remainder, ties broken by a
    seeded shuffle. Rounding up within each stratum instead would inflate small
    strata — 24 strata of 5 at 0.25 would take 48 of 120 runs, not 30.
    """
    target = max(len(sizes), math.ceil(fraction * sum(sizes.values())))
    quotas = {label: max(1, math.floor(fraction * n)) for label, n in sizes.items()}
    tie_order = sorted(sizes)
    random.Random(f"{seed}|allocation").shuffle(tie_order)
    rank = {label: i for i, label in enumerate(tie_order)}
    candidates = sorted(
        (label for label in sizes if quotas[label] < sizes[label]),
        key=lambda label: (-(fraction * sizes[label] - math.floor(fraction * sizes[label])), rank[label]),
    )
    for label in candidates[: max(0, target - sum(quotas.values()))]:
        quotas[label] += 1
    return quotas


def _is_na(v):
    return v is None or (isinstance(v, float) and math.isnan(v))


def resolve_run_subsets(config, sample_df, samples, raw_filepaths):
    """Resolve every subset key in config into the form the rules consume.

    Returns ``{name: None | {"samples": [...], "strata": [...], "raws": [...],
    "spec": {...}}}`` for each name in SUBSET_KEYS. ``raws`` are the absolute
    paths of the selected MS files, in sample order.
    """
    raw_by_sample = {os.path.splitext(os.path.basename(p))[0]: p for p in raw_filepaths}
    resolved = {}
    for name, key in SUBSET_KEYS.items():
        spec = parse_subset_spec(config.get(key), key)
        if spec is None:
            resolved[name] = None
            continue
        picked = select_run_subset(sample_df, samples, spec, key)
        resolved[name] = {
            "spec": spec,
            "samples": [s for s, _ in picked],
            "strata": [label for _, label in picked],
            "raws": [os.path.abspath(raw_by_sample[s]) for s, _ in picked],
        }
    return resolved


def subset_ms_dir(run_dir, name):
    """Directory of symlinks to the selected MS files, staged by the setup module."""
    return os.path.join(run_dir, "ms_subsets", name)


def run_subset_inputs(config, name, experiment_dir, run_dir, raw_filepaths, samples):
    """``(ms_dir, raw_filepaths, samples)`` a DIA-NN step should read for ``name``.

    With no subset configured this returns the experiment's own MS directory and
    the inputs unchanged, so rules render exactly as they did before subsetting
    existed (no spurious params/code rerun triggers on completed runs).
    """
    subset = (config.get("_run_subsets") or {}).get(name)
    if subset is None:
        return os.path.join(experiment_dir, "input/ms_files"), raw_filepaths, samples
    ms_dir = subset_ms_dir(run_dir, name)
    staged = [os.path.join(ms_dir, os.path.basename(p)) for p in subset["raws"]]
    return ms_dir, staged, list(subset["samples"])
