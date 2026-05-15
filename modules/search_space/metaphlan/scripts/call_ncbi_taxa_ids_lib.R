# Pure-function helper for call_ncbi_taxa_ids.R. Kept in a separate file so the
# clade_name parsing logic is reachable from testthat without snakemake@ globals.

# Parse a merged MetaPhlAn profile table and dereplicate species-level NCBI
# taxon IDs whose relative abundance exceeds `threshold` in at least one sample.
#
# Expected input columns: clade_name, NCBI_tax_id, then one column per sample.
# MetaPhlAn emits clade_name as a pipe-delimited string of rank prefixes
# (k__, p__, c__, o__, f__, g__, s__, t__) and NCBI_tax_id as the matching
# pipe-delimited taxon IDs at each rank.
#
# Returns a one-column tibble: ncbi_taxonomy_id (character).
parse_metaphlan_profiles <- function(profiles_df, threshold) {
  sample_names <- names(profiles_df)[3:ncol(profiles_df)]

  profiles_df |>
    tidyr::separate(
      clade_name,
      into = c("kingdom", "phylum", "class", "order", "family", "genus", "species", "strain"),
      sep = "\\|",
      fill = "right"
    ) |>
    dplyr::mutate(
      kingdom = stringr::str_remove(kingdom, "^k__"),
      phylum  = stringr::str_remove(phylum, "^p__"),
      class   = stringr::str_remove(class, "^c__"),
      order   = stringr::str_remove(order, "^o__"),
      family  = stringr::str_remove(family, "^f__"),
      genus   = stringr::str_remove(genus, "^g__"),
      species = stringr::str_remove(species, "^s__"),
      strain  = stringr::str_remove(strain, "^t__")
    ) |>
    tidyr::separate(
      NCBI_tax_id,
      into = c("kingdom_ncbi", "phylum_ncbi", "class_ncbi", "order_ncbi",
               "family_ncbi", "genus_ncbi", "species_ncbi", "strain_ncbi"),
      sep = "\\|",
      fill = "right"
    ) |>
    dplyr::select(species, species_ncbi, dplyr::all_of(sample_names)) |>
    dplyr::filter(species_ncbi != "", !is.na(species_ncbi)) |>
    tidyr::pivot_longer(
      cols      = dplyr::all_of(sample_names),
      names_to  = "sample_name",
      values_to = "relative_abundance"
    ) |>
    dplyr::filter(relative_abundance > threshold) |>
    dplyr::pull(species_ncbi) |>
    unique() |>
    (\(x) tibble::tibble(ncbi_taxonomy_id = x))()
}
