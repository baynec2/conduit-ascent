# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages
start_time <- Sys.time()

conduitR::log_with_timestamp("Running call_ncbi_taxa_ids.R script")

# Defining Inputs
merged_profiles_fp  = snakemake@input[["merged_profiles"]]
# Defining Output
ncbi_taxonomy_id_fp = snakemake@output[["ncbi_taxa_ids"]]
# Defining Config Values
relative_abundance_threshold = snakemake@config$relative_abundance_threshold

conduitR::log_with_timestamp("Reading in merged metaphlan profiles")

merged_profiles = readr::read_delim(merged_profiles_fp)

conduitR::log_with_timestamp(
    paste0("Parsing merged metaphlan profiles ",
    "Filtering to include taxa detected in any sample at relative abundance > ",
    relative_abundance_threshold)
    )

# Extracting sample names, these are from the third column to the last.             
sample_names = names(merged_profiles)[3:ncol(merged_profiles)]
# Splitting taxa and ncbi ids into seperate columns
taxa_split <- merged_profiles |>
  # Split clade_name by pipe
  tidyr::separate(
    clade_name,
    into = c("kingdom","phylum","class","order","family","genus","species","strain"),
    sep = "\\|",
    fill = "right"
  )|>
    # Optional: remove the prefixes (k__, p__, etc.)
  dplyr::mutate(
    kingdom = stringr::str_remove(kingdom, "^k__"),
    phylum  = stringr::str_remove(phylum, "^p__"),
    class   = stringr::str_remove(class, "^c__"),
    order   = stringr::str_remove(order, "^o__"),
    family  = stringr::str_remove(family, "^f__"),
    genus   = stringr::str_remove(genus, "^g__"),
    species = stringr::str_remove(species, "^s__"),
    strain = stringr::str_remove(species, "^t__")
  )|>
    tidyr::separate(
    NCBI_tax_id,
    into = c("kingdom_ncbi","phylum_ncbi","class_ncbi","order_ncbi","family_ncbi","genus_ncbi","species_ncbi","strain_ncbi"),
    sep = "\\|",
    fill = "right"
  )|> 
  dplyr::select(species,species_ncbi,dplyr::all_of(sample_names))|>
  dplyr::filter(species_ncbi != "",
                !is.na(species_ncbi))|>
  tidyr::pivot_longer(cols = dplyr::all_of(sample_names),
                      names_to = "sample_name",
                      values_to = "relative_abundance") |>
                      dplyr::filter(relative_abundance > relative_abundance_threshold)

# Dereplicating organism ids that are beyond theshold in multiple samples.
organism_id = unique(dplyr::pull(taxa_split,species_ncbi))

# Constructing ncbi_taxonomy_id data frame to use as input to the next part of workflow
ncbi_taxonomy_ids = tibble::tibble(organism_id)|>
  dplyr::mutate(organism_type = dplyr::case_when(organism_id %in% c(9606,10090) ~ "host",
                                                 TRUE ~ "microbiome"))
  
conduitR::log_with_timestamp(
    paste0(
        "Writing ncbi taxon ids present at > ",
        relative_abundance_threshold,
        "to", ncbi_taxonomy_id_fp
        )
)

readr::write_delim(ncbi_taxonomy_ids,ncbi_taxonomy_id_fp)
# =============================================================================
# Cleanup and Logging
# =============================================================================
end_time <- Sys.time()
elapsed_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))

conduitR::log_with_timestamp(
  "Completed infer_species_presence.R script. Time taken: %.2f minutes", 
  elapsed_minutes
)

# Close log file connections
sink(type = "message")
sink()
close(zz)
