################################################################################
# Getting Fasta Files from proteome IDS
################################################################################
## Opening Log File 
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

# Now everything from print(), message(), warning() will go into the log file
conduitR::log_with_timestamp("Starting get_fasta_from_proteome_ids.R script")
conduitR::log_with_timestamp("Input file: %s", snakemake@input[["proteome_id_input"]])
conduitR::log_with_timestamp("Output file: %s", snakemake@output[["proteome_id_output"]])
conduitR::log_with_timestamp("Output file: %s", snakemake@output[["fasta_output"]])

## Defining inputs and outputs from snakemake workflow
input_file <- snakemake@input[["proteome_id_input"]]
proteome_id_destination_fp <- snakemake@output[["proteome_id_output"]]
fasta_destination_fp <- snakemake@output[["fasta_output"]]

# Extacting additional uniprot proteome id from the config file.
append_additional_uniprot_proteome_id <- snakemake@config$append_additional_proteome_id

conduitR::log_with_timestamp(paste0(append_additional_uniprot_proteome_id, " specified as additional uniprot proteome id in config file"))

conduitR::log_with_timestamp("Making the database_resources_directory if it doesn't exist.")

# Make the database_resources_directory if it doesn't exist.
if (!dir.exists(dirname(fasta_destination_fp))) {
  dir.create(dirname(fasta_destination_fp))
}

conduitR::log_with_timestamp("Reading proteome IDs from the input file.")

# Read proteome IDs from the input file
proteome_ids <- readr::read_delim(input_file) |>
# If there is NA filter them out. Can occassionly happen when the taxa is present in uniprot, but there is no proteome. 
  dplyr::filter(!is.na(proteome_id))|>
  dplyr::pull(proteome_id)

# Append additional uniprot proteome IDs if specified by the user
if (!isFALSE(append_additional_uniprot_proteome_id)) {
  # Check if the user-specified ID is already in the list of organism IDs
  if (append_additional_uniprot_proteome_id %in% proteome_ids) {
    conduitR::log_with_timestamp(
      paste0(
        "User-specified Uniprot proteome ID ", 
        append_additional_uniprot_proteome_id, 
        " is already present in the data."
      )
    )
  } else {
    # Log that we are appending the user-specified ID
    conduitR::log_with_timestamp(
      paste0(
        "Appending user-specified Uniprot Proteome ID ", 
        append_additional_uniprot_proteome_id, 
        " to the Uniprot Proteome taxonomy IDs to search."
      )
    )
    # Actually append the ID to the list
    proteome_ids <- c(proteome_ids, append_additional_uniprot_proteome_id)
  }
}

conduitR::log_with_timestamp("Starting download of fasta files from the proteome ids.")

start_time <- Sys.time()

# Create a fasta file by hitting uniprot API via R.
conduitR::download_fasta_from_proteome_ids(proteome_ids,
                                        proteome_id_destination_fp = proteome_id_destination_fp,
                                        fasta_destination_fp = fasta_destination_fp)

end_time <- Sys.time()
conduitR::log_with_timestamp("Completed downloading fasta files. Time taken: %.2f minutes", 
    as.numeric(difftime(end_time, start_time, units = "mins")))

# closing logfile connection
conduitR::log_with_timestamp("get_fasta.R script complete.")
# closing logfile connection
sink(type = "message")
sink()
close(zz)