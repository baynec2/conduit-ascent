#!/usr/bin/env python3
"""
coverAllSpectra_greedy.py

Greedy algorithm that selects genomes to cover all profiling spectra.
Adapted from HAPiID (https://github.com/mgtools/HAPiID).

Usage:
    python coverAllSpectra_greedy.py genome2spectrum_dic.json output.tsv

Output TSV columns: genome, nSpectraCovered, cumulative_pct
"""

import json
import sys
import copy
import pandas as pd


def getNextBestGenome(genome2remainingSpectrum_dic):
    """Return the genome ID covering the most remaining spectra."""
    genome2nSpectrum_tuple = [(k, len(v)) for k, v in genome2remainingSpectrum_dic.items()]
    genome2nSpectrum_tuple = sorted(genome2nSpectrum_tuple, key=lambda x: x[1], reverse=True)
    return genome2nSpectrum_tuple[0][0]


def updateGenome2remainingSpectrum_dic(genome2remainingSpectrum_dic, remaining_spectra, covered_spectra):
    """Remove covered spectra from remaining and update per-genome intersections."""
    remaining_spectra = list(set(remaining_spectra).difference(covered_spectra))
    for genome in genome2remainingSpectrum_dic:
        genome2remainingSpectrum_dic[genome] = set(
            genome2remainingSpectrum_dic[genome]
        ).intersection(set(remaining_spectra))
    return genome2remainingSpectrum_dic, remaining_spectra


def greedy_cover(genome2spectrum_dic):
    """Run greedy set-cover; return list of (genome, n_spectra_covered_cum) tuples."""
    all_spectra = []
    for spectra in genome2spectrum_dic.values():
        all_spectra.extend(spectra)
    all_spectra = list(set(all_spectra))

    remaining_spectra = copy.deepcopy(all_spectra)
    genome2remainingSpectrum_dic = copy.deepcopy(genome2spectrum_dic)

    selected = []
    covered_so_far = []

    while remaining_spectra and genome2remainingSpectrum_dic:
        nxt = getNextBestGenome(genome2remainingSpectrum_dic)
        covered = genome2remainingSpectrum_dic[nxt]
        # If the best remaining genome covers nothing, no further progress is
        # possible — stop rather than emit zero-coverage rows forever.
        if len(covered) == 0:
            break
        covered_so_far.extend(covered)
        covered_so_far = list(set(covered_so_far))
        selected.append((nxt, len(covered_so_far)))
        genome2remainingSpectrum_dic, remaining_spectra = updateGenome2remainingSpectrum_dic(
            genome2remainingSpectrum_dic, remaining_spectra, covered
        )

    return selected


def main(argv):
    if len(argv) != 3:
        print(
            "Usage: python coverAllSpectra_greedy.py genome2spectrum_dic.json output.tsv"
        )
        return 1

    genome2spectrum_dic_f, out_f = argv[1], argv[2]

    with open(genome2spectrum_dic_f) as in_f:
        genome2spectrum_dic = json.load(in_f)

    selected = greedy_cover(genome2spectrum_dic)

    if not selected:
        # Empty first-pass result (no spectra / no taxa): emit a header-only
        # selection so downstream resolves to an empty (no-detection) conduit
        # instead of crashing on the cumulative-percentage step below.
        with open(out_f, "w") as _out_f:
            _out_f.write("genome\tnSpectraCovered\tcumulative_pct\n")
        return 0

    with open(out_f, "w") as _out_f:
        _out_f.write("genome\tnSpectraCovered\n")
        for genome, n_cum in selected:
            _out_f.write(f"{genome}\t{n_cum}\n")

    # Add cumulative percentage column
    covered_spectra_df = pd.read_csv(out_f, sep="\t")
    total = covered_spectra_df["nSpectraCovered"].iloc[-1]
    covered_spectra_df["cumulative_pct"] = [
        (item * 100) / total for item in covered_spectra_df["nSpectraCovered"]
    ]
    covered_spectra_df.to_csv(out_f, sep="\t", index=False)

    # Optional: generate visualization if plotly is available
    try:
        import plotly.graph_objects as go

        fig = go.Figure()
        fig.add_trace(
            go.Scatter(
                x=covered_spectra_df["genome"],
                y=covered_spectra_df["cumulative_pct"],
                mode="lines+markers",
            )
        )
        fig.update_layout(
            yaxis=dict(range=[0, 100]),
            width=max(400, len(covered_spectra_df) * 5),
            xaxis_title="top N genomes",
            yaxis_title="percentage of cumulative spectra covered",
            height=800,
            title={
                "text": "Cumulative % of spectra covered by genomes (greedy approach)",
                "y": 0.96,
                "x": 0.02,
            },
        )
        out_dir = out_f.rsplit("/", 1)[0] + "/"
        out_fname = out_f.rsplit("/", 1)[-1].replace(".tsv", "")
        fig.write_image(out_dir + out_fname + "_cumulativePercentSpecsCovered.png")
    except Exception:
        pass  # Visualization is optional

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
