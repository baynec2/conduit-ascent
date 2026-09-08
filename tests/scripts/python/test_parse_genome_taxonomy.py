"""Tests for modules/search_space/genomes/scripts/parse_genome_taxonomy.py.
Targets the pure transforms: filter_to_selected_genomes, assign_organism_ids,
fill_taxonomy_defaults, and the build_genome_taxonomy end-to-end pipe."""

import pandas as pd
import pytest

import parse_genome_taxonomy as pmt


def make_tax_df(rows):
    """rows is a list of dicts; each must include 'genome'. Missing taxonomy
    columns are intentionally allowed — tests cover the fallback path."""
    return pd.DataFrame(rows)


def test_assign_organism_ids_sorts_alphabetically():
    df = make_tax_df([
        {"genome": "B_subtilis"},
        {"genome": "A_baumannii"},
        {"genome": "E_coli"},
    ])
    out = pmt.assign_organism_ids(df)
    assert list(out["genome"]) == ["A_baumannii", "B_subtilis", "E_coli"]
    assert list(out["organism_id"]) == [1, 2, 3]


def test_filter_to_selected_genomes_keeps_only_listed():
    df = make_tax_df([
        {"genome": "G1"},
        {"genome": "G2"},
        {"genome": "G3"},
    ])
    out = pmt.filter_to_selected_genomes(df, {"G1", "G3"})
    assert set(out["genome"]) == {"G1", "G3"}


def test_fill_taxonomy_defaults_adds_missing_columns_as_NA():
    # Real-world failure mode: user-provided taxonomy.txt with only some ranks
    # filled. The pure transform must NA-fill every missing rank rather than
    # crash on KeyError downstream.
    df = make_tax_df([{"genome": "G1", "species": "E. coli"}])
    out = pmt.fill_taxonomy_defaults(df)
    for col in ("domain", "kingdom", "phylum", "class", "order", "family", "genus"):
        assert col in out.columns
        assert (out[col] == "NA").all()
    assert out["species"].iloc[0] == "E. coli"


def test_fill_taxonomy_defaults_adds_system_columns():
    df = make_tax_df([{"genome": "G1"}])
    out = pmt.fill_taxonomy_defaults(df)
    assert (out["proteome_id"] == "NA").all()
    assert (out["proteome_type"] == "NA").all()
    assert (out["download_info"] == "user_provided").all()


def test_build_genome_taxonomy_end_to_end():
    df = make_tax_df([
        {"genome": "B_subtilis", "species": "Bacillus subtilis", "phylum": "Bacillota"},
        {"genome": "A_baumannii", "species": "Acinetobacter baumannii"},
        {"genome": "E_coli",     "species": "Escherichia coli"},
    ])
    out = pmt.build_genome_taxonomy(df)
    # Output schema must match the OUT_COLS contract exactly.
    assert list(out.columns) == pmt.OUT_COLS
    # Alphabetical organism_id assignment.
    assert list(out["genome"]) == ["A_baumannii", "B_subtilis", "E_coli"]
    assert list(out["organism_id"]) == [1, 2, 3]


def test_build_genome_taxonomy_filters_then_assigns_ids():
    # IDs must be assigned AFTER filtering — otherwise they enumerate rows
    # whose proteins are never in the database. Pinning this contract.
    df = make_tax_df([
        {"genome": "G1"},
        {"genome": "G2"},
        {"genome": "G3"},
    ])
    out = pmt.build_genome_taxonomy(df, selected={"G2", "G3"})
    assert list(out["genome"]) == ["G2", "G3"]
    assert list(out["organism_id"]) == [1, 2]


def test_build_genome_taxonomy_requires_genome_column():
    df = make_tax_df([{"species": "E. coli"}])
    with pytest.raises(ValueError, match="genome"):
        pmt.build_genome_taxonomy(df)
