#!/usr/bin/env python3
"""
coverAllSpectra_greedy.py

Greedy set-cover selection of genomes (or taxa) over identified spectra.

Repeatedly take whichever genome explains the most spectra not yet explained,
until nothing further can be added. This is the standard greedy approximation
to set cover; the method is the one described in HAPiID:

    Stamboulian M, Li S, Ye Y. Using high-abundance proteins as guides for fast
    and effective peptide/protein identification from human gut metaproteomic
    data. Microbiome 9, 80 (2021). doi:10.1186/s40168-021-01035-8

The implementation here is our own, written against the algorithm and this
module's tests rather than derived from HAPiID's source. See
THIRD_PARTY_NOTICES.md.

Usage:
    python coverAllSpectra_greedy.py genome2spectrum_dic.json output.tsv

Output TSV columns: genome, nSpectraCovered, cumulative_pct
"""

import json
import sys


def greedy_cover(genome2spectrum_dic):
    """Select genomes greedily by uncovered-spectrum count.

    Takes {genome_id: [spectrum_id, ...]} and returns a list of
    (genome_id, cumulative_spectra_covered) in selection order.

    Spectrum lists are collapsed to sets on entry, so a genome listing the same
    spectrum twice gains nothing from the repeat. Both callers already
    deduplicate (build_taxon_spectrum_mapping.py, build_genome_spectrum_mapping.py),
    but relying on that silently would make a duplicate leak into a future
    caller show up as a quietly worse genome selection rather than as an error.

    Ties are broken by iteration order — the first genome reaching the maximum
    wins — which makes the output deterministic for a given input file, since
    json.load preserves the order the mapping was written in.
    """
    uncovered_by_genome = {g: set(spectra) for g, spectra in genome2spectrum_dic.items()}

    selected = []
    total_covered = 0

    while True:
        best = None
        best_size = 0
        for genome, uncovered in uncovered_by_genome.items():
            if len(uncovered) > best_size:
                best, best_size = genome, len(uncovered)

        # No genome adds anything further: either everything is covered, or the
        # remainder is unreachable. Either way, stop rather than emit rows that
        # advance coverage by zero.
        if best is None:
            return selected

        newly_covered = uncovered_by_genome.pop(best)
        total_covered += len(newly_covered)
        selected.append((best, total_covered))

        for uncovered in uncovered_by_genome.values():
            uncovered -= newly_covered


def write_selection(selected, out_f):
    """Write the selection table, including cumulative percentage of coverage."""
    with open(out_f, "w") as fh:
        fh.write("genome\tnSpectraCovered\tcumulative_pct\n")
        if not selected:
            # Empty first-pass result (no spectra / no taxa). A header-only file
            # lets downstream resolve to an empty (no-detection) conduit instead
            # of dividing by a total of zero here.
            return
        total = selected[-1][1]
        for genome, n_cumulative in selected:
            fh.write(f"{genome}\t{n_cumulative}\t{n_cumulative * 100 / total}\n")


def main(argv):
    if len(argv) != 3:
        print("Usage: python coverAllSpectra_greedy.py genome2spectrum_dic.json output.tsv")
        return 1

    in_f, out_f = argv[1], argv[2]

    with open(in_f) as fh:
        genome2spectrum_dic = json.load(fh)

    write_selection(greedy_cover(genome2spectrum_dic), out_f)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
