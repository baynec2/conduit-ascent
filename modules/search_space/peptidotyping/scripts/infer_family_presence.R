# =============================================================================
# Setup and Logging
# =============================================================================
# Open log file for both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

# Record start time
start_time <- Sys.time()
conduitR::log_with_timestamp("Starting infer_family_presence.R script")

# =============================================================================
#  Setting up Input and Output Files
# =============================================================================

# Test
#first_pass_diann_parquet <-  "experiments/defined_community/input/database_resources/proteotyping/report-family.parquet"
#presence_min_peptides <-2

# Input files
first_pass_diann_parquet <- snakemake@input[["first_pass_diann"]]

conduitR::log_with_timestamp("Input files: %s", first_pass_diann_parquet)

# Config values:
presence_min_peptides <- snakemake@config[["presence_min_peptides"]]


# Output files
ncbi_taxonomy_id_fp <- snakemake@output[["ncbi_taxonomy_id"]]

conduitR::log_with_timestamp("Reading in config values to specify theshold")

presence_min_peptides <- snakemake@config[["presence_min_peptides"]]

conduitR::log_with_timestamp("Output files: %s", ncbi_taxonomy_id_fp)

# =============================================================================
#  Generating First Pass Search Taxa Metrics
# =============================================================================
# Reading in precursor file with parquet
conduitR::log_with_timestamp("Reading in: %s", first_pass_diann_parquet)
precursors <- arrow::read_parquet(first_pass_diann_parquet)

# Summarising precursors to peptide
conduitR::log_with_timestamp("Summing precursors to peptides")

# Count unique peptides per detected_taxonomy vs true taxonomy
  pep <- precursors |> 
    dplyr::filter(Proteotypic == 1) |> 
    dplyr::mutate(detected_taxonomy = gsub(".*_", "", Protein.Names),
    ncbi_taxonomy_id = gsub(".*\\|", "", Protein.Group),
    ) |> 
    dplyr::group_by(Run, detected_taxonomy, ncbi_taxonomy_id, Stripped.Sequence) |> 
    dplyr::summarise(sum_intensity = sum(Precursor.Normalised), .groups = "drop")
  
  # Count number of unique peptides per taxon per sample
  sum_pep <- pep |> 
    dplyr::group_by(detected_taxonomy, ncbi_taxonomy_id) |> 
    dplyr::summarise(n_peptides = dplyr::n(), .groups = "drop") |>
    dplyr::filter(n_peptides >= presence_min_peptides)|>
    dplyr::select(ncbi_taxonomy_id,detected_taxonomy)|>
    dplyr::distinct()

  # Unique  
readr::write_delim(sum_pep, ncbi_taxonomy_id_fp)


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