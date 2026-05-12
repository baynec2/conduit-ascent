################################################################################
# Getting uniprot proteome IDs from NCBI organism IDs
################################################################################
## Opening Log File 
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

start_time <- Sys.time()

# Now everything from print(), message(), warning() will go into the log file
conduitR::log_with_timestamp("Starting get_uniprot_proteome_ids.R script")
conduitR::log_with_timestamp("Input file: %s", snakemake@input[[1]])
conduitR::log_with_timestamp("Output file: %s", snakemake@output[[1]])

## Defining inputs and outputs from snakemake workflow
input_file <- snakemake@input[[1]]
proteome_ids_fp <- snakemake@output[[1]]

# Extracting additional ncbi taxa ids to append if specified
append_additional_ncbi_taxa_id <- snakemake@config$append_additional_ncbi_taxa_id

conduitR::log_with_timestamp("Making the database_resources_directory if it doesn't exist.")

conduitR::log_with_timestamp("Reading organism IDs from the input file.")
# Read organism IDs from the input file. Tolerates either a single-column
# format (just the numeric ID under any header, e.g. `ncbi_taxa_id`) or a
# multi-column TSV where the first column carries the IDs and additional
# columns carry user metadata (e.g. `organism_id\tsource`). The previous
# implementation used `readr::read_lines + as.integer`, which silently
# coerced lines like `820\tBacteroides_uniformis_ATCC_8492` to NA and then
# queried UniProt with NAs (HTTP 400).
organism_ids <- readr::read_tsv(input_file, show_col_types = FALSE) |>
  dplyr::pull(1) |>
  as.integer() |>
  (\(x) x[!is.na(x)])() |>
  unique()

# Append additional NCBI taxonomic IDs if specified by the user
if (!isFALSE(append_additional_ncbi_taxa_id)) {
  # Check if the user-specified ID is already in the list of organism IDs
  if (append_additional_ncbi_taxa_id %in% organism_ids) {
    conduitR::log_with_timestamp(
      paste0(
        "User-specified NCBI organism ID ", 
        append_additional_ncbi_taxa_id, 
        " is already present in the data."
      )
    )
  } else {
    # Log that we are appending the user-specified ID
    conduitR::log_with_timestamp(
      paste0(
        "Appending user-specified NCBI organism ID ", 
        append_additional_ncbi_taxa_id, 
        " to the NCBI taxonomy IDs to search."
      )
    )
    # Actually append the ID to the list
    organism_ids <- c(organism_ids, append_additional_ncbi_taxa_id)
  }
}

conduitR::log_with_timestamp("Finding the uniprot proteome ID cooresponding to each NCBI taxonomy ID")

# Using conduitR to get proteome ids corresponding to an organism id
proteome_ids_df <- conduitR::get_proteome_ids_from_organism_ids(organism_ids)

conduitR::log_with_timestamp("Searching for higher-quality UniProt proteome IDs from closely related taxa")

# UniProt API results may not always provide the optimal proteome for a given NCBI taxonomic ID.
# Returned entries can be redundant, absent from UniProtKB, or assigned to strains without curated proteomes.
# We therefore identify higher-quality proteomes from related taxa and select the best available match.
best_proteome_ids <- conduitR::get_better_proteome_ids(proteome_ids_df)

# Renaming to be compatible with uniprot_proteome_ids search space workflow. Essentially this will
# extract proteome_id from the dataframe column "proteome_id".
best_proteome_ids <- best_proteome_ids |>
dplyr::rename(initial_proteome_id = proteome_id,
proteome_id = selected_proteome_id)

# Save to file
conduitR::log_with_timestamp("Saving proteome ids to file.")

readr::write_tsv(best_proteome_ids,proteome_ids_fp)

end_time <- Sys.time()

conduitR::log_with_timestamp("Completed get_uniprot_proteome_ids.R script. Time taken: %.2f minutes", 
    as.numeric(difftime(end_time, start_time, units = "mins")))

# closing logfile connection
conduitR::log_with_timestamp("get_uniprot_proteome_ids.R script complete.")
# closing logfile connection
sink(type = "message")
sink()
close(zz)