"""Tests for modules/search_space/unipept_hapid/scripts/build_taxon_spectrum_mapping.py.
Targets the pure transform `build_taxon_spectrum_mapping(df)`."""

import pandas as pd
import pytest

import build_taxon_spectrum_mapping as btsm


def make_diann_df(rows):
    """Build a synthetic DIA-NN parquet-shape df. Each row dict supplies
    Run, Precursor.Id, Protein.Group, Proteotypic."""
    return pd.DataFrame(rows)


def test_groups_unique_spectra_per_taxon():
    df = make_diann_df([
        {"Run": "r1", "Precursor.Id": "p1", "Protein.Group": "umgap|abc|562", "Proteotypic": 1},
        {"Run": "r1", "Precursor.Id": "p2", "Protein.Group": "umgap|def|562", "Proteotypic": 1},
        {"Run": "r1", "Precursor.Id": "p3", "Protein.Group": "umgap|ghi|1280", "Proteotypic": 1},
    ])
    out = btsm.build_taxon_spectrum_mapping(df)
    assert set(out.keys()) == {"562", "1280"}
    assert out["562"] == ["r1||p1", "r1||p2"]
    assert out["1280"] == ["r1||p3"]


def test_drops_non_proteotypic_precursors():
    df = make_diann_df([
        {"Run": "r1", "Precursor.Id": "p1", "Protein.Group": "umgap|abc|562", "Proteotypic": 1},
        {"Run": "r1", "Precursor.Id": "p2", "Protein.Group": "umgap|def|562", "Proteotypic": 0},
    ])
    out = btsm.build_taxon_spectrum_mapping(df)
    assert out == {"562": ["r1||p1"]}


def test_extracts_taxon_as_last_pipe_segment():
    # umgap|seq_id|taxid — rsplit("|", n=1) returns the trailing taxid even
    # when seq_id contains additional pipes (shouldn't happen, but defensive).
    df = make_diann_df([
        {"Run": "r1", "Precursor.Id": "p1", "Protein.Group": "umgap|abc|562", "Proteotypic": 1},
        {"Run": "r1", "Precursor.Id": "p2", "Protein.Group": "umgap|with|pipes|999", "Proteotypic": 1},
    ])
    out = btsm.build_taxon_spectrum_mapping(df)
    assert "562" in out
    assert "999" in out


def test_deduplicates_spectrum_ids_within_taxon():
    # Same Run+Precursor.Id appearing twice (e.g., two FASTA proteins with
    # same LCA taxid) should produce one spectrum_id, not two.
    df = make_diann_df([
        {"Run": "r1", "Precursor.Id": "p1", "Protein.Group": "umgap|abc|562", "Proteotypic": 1},
        {"Run": "r1", "Precursor.Id": "p1", "Protein.Group": "umgap|def|562", "Proteotypic": 1},
    ])
    out = btsm.build_taxon_spectrum_mapping(df)
    assert out["562"] == ["r1||p1"]


def test_spectrum_ids_are_sorted():
    # Output contract: per-taxon spectrum lists are sorted (greedy cover
    # downstream is order-sensitive on deterministic iteration order).
    df = make_diann_df([
        {"Run": "r2", "Precursor.Id": "p2", "Protein.Group": "umgap|abc|562", "Proteotypic": 1},
        {"Run": "r1", "Precursor.Id": "p1", "Protein.Group": "umgap|def|562", "Proteotypic": 1},
    ])
    out = btsm.build_taxon_spectrum_mapping(df)
    assert out["562"] == sorted(out["562"])


def test_empty_after_proteotypic_filter_returns_empty_dict():
    df = make_diann_df([
        {"Run": "r1", "Precursor.Id": "p1", "Protein.Group": "umgap|abc|562", "Proteotypic": 0},
    ])
    out = btsm.build_taxon_spectrum_mapping(df)
    assert out == {}
