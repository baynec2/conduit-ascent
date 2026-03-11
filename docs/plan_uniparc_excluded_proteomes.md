# Plan: Incorporating Uniparc-Only (Excluded) Proteomes in uniprot_proteome_ids Module

## Context

- **Current behaviour**: The module downloads proteomes from **UniProtKB** only (`conduitR::download_fasta_from_proteome_ids` → `rest.uniprot.org/uniprotkb/stream`).
- **Gap**: Some proteomes are **Excluded** from UniProtKB (e.g. redundant, low-quality, or only present in the archive). They appear in `taxonomy.txt` with `proteome_type == "Excluded"` and `download_info == "not_downloaded"`.
- **Source of Excluded IDs**: `experiments/<exp>/input/database_resources/taxonomy.txt`, column `proteome_id`, where `proteome_type == "Excluded"`.
- **Uniparc**: Excluded proteomes can be downloaded from **Uniparc** via the same proteome ID, e.g.  
  `https://rest.uniprot.org/uniparc/stream?query=proteome:<UPID>&format=fasta`  
  (verified for e.g. UP000014150).

## Options

### Option A: Fallback in existing download step (recommended if conduitR can be extended)

**Idea**: In the same step that downloads from UniProtKB, for each proteome ID try UniProtKB first; if the response is empty, call the Uniparc stream for that proteome and append sequences.

**Pros**  
- No workflow reordering; no new Snakemake rules.  
- Single place for “get FASTA for this proteome”.  
- No need for taxonomy.txt to exist first.

**Cons**  
- Requires logic (and possibly API) to detect “no UniProtKB sequences” per proteome.  
- Likely needs changes in **conduitR** (e.g. `download_fasta_from_proteome_ids` or a new helper that tries Uniparc on empty UniProtKB result).

**Implementation sketch**  
- In R (conduitR): for each `proteome_id`, call UniProtKB stream; if result is empty or zero sequences, call Uniparc stream for that proteome and append to the FASTA (and optionally record source as “uniparc” in a side table).

---

### Option B: New rule after `get_taxonomy` (recommended if Option A is not feasible)

**Idea**: Add a rule that runs after `get_taxonomy`, reads `taxonomy.txt`, filters `proteome_type == "Excluded"`, downloads those proteomes from Uniparc, then merges the result into the database FASTA (and optionally tracks which proteomes came from Uniparc).

**Pros**  
- Uses existing taxonomy and column `proteome_type`; no need to change conduitR download logic.  
- Clear separation: UniProtKB download → taxonomy → Uniparc backfill.

**Cons**  
- New rule and dependency; downstream rules must consume the merged FASTA (and possibly an updated proteome/ID list).  
- Uniparc FASTA headers are UPI-only (e.g. `>UPI00015329C3 status=active`), so downstream **organism_id** (and protein_info) must be handled (see below).

**Implementation sketch**

1. **New rule** (e.g. `get_fasta_from_uniparc_excluded`) in `uniprot_proteome_ids.smk`:
   - **Inputs**: `taxonomy.txt`, `database.fasta` (from `get_fasta_from_proteome_ids`).
   - **Outputs**: e.g. `database.fasta` (merged) or a dedicated `uniparc_excluded.fasta` that a small merge step combines with `database.fasta` so that a single path is used downstream.

2. **Script** (e.g. `scripts/get_fasta_from_uniparc_excluded.R`):
   - Read `taxonomy.txt`, filter `proteome_type == "Excluded"`, get unique `proteome_id` (and keep `organism_id` for header enrichment).
   - For each Excluded proteome ID, request  
     `https://rest.uniprot.org/uniparc/stream?query=proteome:<proteome_id>&format=fasta`.
   - Optionally rewrite headers to include `organism_id` (and `proteome_id`) so that `extract_fasta_info` (or a small adapter) can attach taxonomy in `get_protein_info_from_fasta` (e.g. `>UPI...|organism_id=123|proteome_id=UP0000...` or a format conduitR understands).
   - Append Uniparc sequences to the main FASTA (or write `uniparc_excluded.fasta` and merge in a separate rule/script).

3. **Downstream**:
   - `get_protein_info_from_fasta` (and any step that expects a single `database.fasta`) should use the **merged** FASTA (either produced by the new rule or by a small merge rule that combines `database.fasta` + `uniparc_excluded.fasta`).
   - Ensure `protein_info` can resolve `organism_id` for Uniparc-derived entries (via header enrichment or a join key from taxonomy).

---

## Downstream: FASTA headers and `organism_id`

- **UniProtKB** headers usually carry taxonomy (e.g. `OX=NCBI_TaxID=...`), which conduitR’s `extract_fasta_info` likely uses for the `left_join` with `taxonomy.txt` on `organism_id`.
- **Uniparc** stream returns headers like `>UPI00015329C3 status=active` (no organism_id).

So for Uniparc-sourced sequences you must either:

1. **Enrich headers when writing**  
   When saving Uniparc FASTA, add `organism_id` (and optionally `proteome_id`) into the header (e.g. from taxonomy), and ensure conduitR’s `extract_fasta_info` (or a wrapper) can parse this and join to taxonomy.

2. **Or extend conduitR**  
   e.g. a second path: “protein_info from Uniparc FASTA + mapping proteome_id → organism_id from taxonomy”, then merge with the main protein_info.

---

## Recommendation

- **Preferred**: **Option A** if you can extend conduitR so that the existing download step tries Uniparc when UniProtKB returns no sequences; minimal pipeline change and one place for “proteome → FASTA”.
- **Otherwise**: **Option B** (new rule after `get_taxonomy`, filter `proteome_type == "Excluded"`, download from Uniparc, merge FASTA and handle header/organism_id as above).

## Files to touch (for Option B)

| Component | File(s) |
|-----------|--------|
| New rule + script | `modules/search_space/uniprot_proteome_ids/uniprot_proteome_ids.smk`, `.../scripts/get_fasta_from_uniparc_excluded.R` |
| Merge / downstream | Either merge in the new script into one `database.fasta`, or add a merge rule and point `get_protein_info_from_fasta` (and any other consumers) at the merged FASTA |
| Header / protein_info | Enrich Uniparc headers with `organism_id` (from taxonomy) and/or adapt `get_protein_info_from_fasta` / conduitR so Uniparc entries get correct `organism_id` for the taxonomy join |

## Summary

- **Excluded proteome IDs**: from `taxonomy.txt` where `proteome_type == "Excluded"`.
- **Download**: Uniparc stream API `rest.uniprot.org/uniparc/stream?query=proteome:<UPID>&format=fasta`.
- **Integration**: Either fallback in current download (Option A) or new post-taxonomy rule + merge (Option B).
- **Critical**: Ensure Uniparc sequences have `organism_id` (and proteome_id if needed) for `protein_info` and taxonomy join (header enrichment or conduitR extension).
