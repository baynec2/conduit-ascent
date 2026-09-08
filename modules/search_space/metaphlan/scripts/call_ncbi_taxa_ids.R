# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages
start_time <- Sys.time()

conduitR::log_with_timestamp("Running call_ncbi_taxa_ids.R script")

# Pure-function helper (parse_metaphlan_profiles) lives in a sibling file so
# it's reachable from testthat without snakemake@ globals.
snakemake@source("call_ncbi_taxa_ids_lib.R")

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

ncbi_taxonomy_ids <- parse_metaphlan_profiles(merged_profiles, relative_abundance_threshold)
  
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
