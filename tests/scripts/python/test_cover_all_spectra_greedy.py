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
