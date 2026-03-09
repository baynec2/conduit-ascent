################################################################################
# Creating Taxonomy File from user provided MAG metadata file
################################################################################
## Opening Log File 
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

# Get input and output files from Snakemake workflow
mag_metadata_fp <- snakemake@input[["mag_metadata"]]
taxonomy_fp <- snakemake@output[["taxonomy"]]

start_time <- Sys.time()
conduitR::log_with_timestamp("Running get_mag_taxonomy.R script")
conduitR::log_with_timestamp(paste0("Input file: ", mag_metadata_fp))
conduitR::log_with_timestamp(paste0("Output file: ", taxonomy_fp))

conduitR::log_with_timestamp("Reading MAG taxonomy from input file.")
# Read MAG taxonomy file
mag_metadata <- readr::read_delim(mag_metadata_fp)

# Error if organism id not found
if(!("organism_id" %in% names(mag_metadata))) {
  stop("organism_id must be a provided column in the mag metadata")
}

# Extract organism IDs (NCBI taxonomy IDs)
organism_ids <- mag_metadata |>
  dplyr::pull(organism_id) |>
  as.character() |>
  unique()

conduitR::log_with_timestamp(paste0("Found ", length(organism_ids), " unique organism IDs in MAG metadata."))

conduitR::log_with_timestamp("Getting Full NCBI Taxonomy corresponding to organism_ids from NCBI API.")

# Pull all taxonomy information from NCBI API
taxonomy <- conduitR::get_ncbi_taxonomy(organism_ids)

conduitR::log_with_timestamp("Finished downloading Taxonomy Information from NCBI API.")

# Join with MAG taxonomy to preserve MAG information
taxonomy <- taxonomy |>
  # Set organism_type (MAGs are typically microbiome)
  dplyr::mutate(organism_type = dplyr::case_when(
    organism_id %in% c("9606", "10090") ~ "host",
    TRUE ~ "microbiome"
  ),
  download_info = "MAG_user_supplied_ncbi_taxonomy")

conduitR::log_with_timestamp("Writing Taxonomy Information to file.")
# Writing to file.
readr::write_delim(taxonomy, taxonomy_fp)

end_time <- Sys.time()

conduitR::log_with_timestamp("Completed get_mag_taxonomy.R script. Time taken: %.2f minutes", 
    as.numeric(difftime(end_time, start_time, units = "mins")))

# closing logfile connection
sink(type = "message")
sink()
close(zz)

