"""Shared pytest fixtures + sys.path wiring so tests can import the snakemake
scripts as ordinary Python modules. The scripts gate their snakemake-bound
main() with `if __name__ == "__main__" or "snakemake" in globals():` so plain
imports run no I/O — only the pure functions get loaded."""

import os
import sys

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))

for rel in [
    "modules/search_space/_shared/scripts",
    "modules/search_space/genomes/scripts",
    "modules/search_space/unipept_hapid/scripts",
    "modules/_shared",
]:
    p = os.path.join(REPO_ROOT, rel)
    if p not in sys.path:
        sys.path.insert(0, p)
