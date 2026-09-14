"""Tests for the greedy set-cover helper used by the HAPiID-style genome
selection (modules/search_space/_shared/scripts/coverAllSpectra_greedy.py).
Tests the pure `greedy_cover` function — file I/O and CLI are tested by the
hapid integration test."""

import pytest

import coverAllSpectra_greedy as cag


def test_greedy_cover_single_genome_covers_all():
    out = cag.greedy_cover({"G1": ["s1", "s2", "s3"]})
    assert out == [("G1", 3)]


def test_greedy_cover_picks_highest_coverage_first():
    # G2 covers 3, G1 covers 2 — G2 should be selected first.
    g2s = {
        "G1": ["s1", "s2"],
        "G2": ["s1", "s2", "s3"],
    }
    out = cag.greedy_cover(g2s)
    assert out[0][0] == "G2"
    # Two-element output: G2 covers 3 of 3; G1 contributes nothing new, so
    # the loop should stop after G2 rather than picking redundant genomes.
    assert out == [("G2", 3)]


def test_greedy_cover_covers_disjoint_genomes_iteratively():
    g2s = {
        "G1": ["s1", "s2"],
        "G2": ["s3", "s4"],
        "G3": ["s5"],
    }
    out = cag.greedy_cover(g2s)
    genomes_selected = [g for g, _ in out]
    # All three are required since they're disjoint; the cumulative count
    # must reach the total (5) exactly.
    assert set(genomes_selected) == {"G1", "G2", "G3"}
    assert out[-1][1] == 5
    # Cumulative counts must be monotonically non-decreasing.
    cums = [n for _, n in out]
    assert cums == sorted(cums)


def test_greedy_cover_handles_overlapping_genomes():
    # G1 fully contains G2's spectra. G1 alone covers all → G2 not selected.
    g2s = {
        "G1": ["s1", "s2", "s3"],
        "G2": ["s1", "s2"],
    }
    out = cag.greedy_cover(g2s)
    assert out == [("G1", 3)]


def test_greedy_cover_empty_input():
    assert cag.greedy_cover({}) == []


def test_greedy_cover_genome_with_no_spectra_is_not_selected():
    # An empty genome should never advance coverage; loop should terminate.
    g2s = {
        "G1": ["s1", "s2"],
        "G2": [],
    }
    out = cag.greedy_cover(g2s)
    assert [g for g, _ in out] == ["G1"]


def test_greedy_cover_counts_duplicate_spectra_once():
    # A genome listing the same spectrum repeatedly must not appear to cover
    # more than it does. G1 has 4 entries but only 2 distinct spectra, so G2's
    # 3 distinct spectra should win the first pick.
    g2s = {
        "G1": ["s1", "s1", "s2", "s2"],
        "G2": ["s3", "s4", "s5"],
    }
    out = cag.greedy_cover(g2s)
    assert out[0][0] == "G2"
    # Total coverage is 5 distinct spectra, not the 7 raw list entries.
    assert out[-1][1] == 5


def test_greedy_cover_ties_resolve_by_iteration_order():
    # Equal coverage must resolve deterministically to the first genome seen,
    # so a given input file always yields the same selection.
    tie = {"G1": ["s1", "s2"], "G2": ["s3", "s4"]}
    assert cag.greedy_cover(tie)[0][0] == "G1"
    assert cag.greedy_cover({"G2": ["s3", "s4"], "G1": ["s1", "s2"]})[0][0] == "G2"


def test_greedy_cover_cumulative_counts_are_monotonic_and_complete():
    g2s = {
        "G1": ["s1", "s2", "s3"],
        "G2": ["s3", "s4"],
        "G3": ["s5"],
        "G4": ["s1"],
    }
    out = cag.greedy_cover(g2s)
    cums = [n for _, n in out]
    assert cums == sorted(cums)
    # Every reachable spectrum is accounted for, and no genome is picked twice.
    assert out[-1][1] == len({s for spectra in g2s.values() for s in spectra})
    assert len({g for g, _ in out}) == len(out)


def test_write_selection_empty_emits_header_only(tmp_path):
    out_f = tmp_path / "sel.tsv"
    cag.write_selection([], str(out_f))
    assert out_f.read_text() == "genome\tnSpectraCovered\tcumulative_pct\n"


def test_write_selection_cumulative_pct_reaches_100(tmp_path):
    out_f = tmp_path / "sel.tsv"
    cag.write_selection([("G1", 3), ("G2", 5)], str(out_f))
    rows = [ln.split("\t") for ln in out_f.read_text().strip().split("\n")[1:]]
    assert [r[0] for r in rows] == ["G1", "G2"]
    assert float(rows[0][2]) == 60.0
    assert float(rows[-1][2]) == 100.0
