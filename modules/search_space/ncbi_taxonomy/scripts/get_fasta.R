################################################################################
# Getting Fasta Files from organism IDS
################################################################################
## Opening Log File 
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

# Now everything from print(), message(), warning() will go into the log file
conduitR::log_with_timestamp("Starting get_fasta.R script")
conduitR::log_with_timestamp("Input file: %s", snakemake@input[[1]])
conduitR::log_with_timestamp("Output file: %s", snakemake@output[[1]])
conduitR::log_with_timestamp("Output file: %s", snakemake@output[[2]])

## Defining inputs and outputs from snakemake workflow
input_file <- snakemake@input[[1]]
proteome_id_destination_fp <- snakemake@output[[1]]
fasta_destination_fp <- snakemake@output[[2]]


conduitR::log_with_timestamp("Making the database_resources_directory if it doesn't exist.")

# Make the database_resources_directory if it doesn't exist.
if (!dir.exists(dirname(fasta_destination_fp))) {
  dir.create(dirname(fasta_destination_fp))
}

conduitR::log_with_timestamp("Reading organism IDs from the input file.")
# Read organism IDs from the input file
organism_ids <- readr::read_delim(input_file) |>
  dplyr::pull(organism_id)

conduitR::log_with_timestamp("Starting download of fasta files from the organism IDs.")

start_time <- Sys.time()

# Create a fasta file by hitting uniprot API via R.
# Warnings are suppressed because the df containing reference proteomes has
# a number of rows that could change with future updates
suppressWarnings(conduitR::download_fasta_from_organism_ids(organism_ids,
                                        proteome_id_destination_fp = proteome_id_destination_fp,
                                        fasta_destination_fp = fasta_destination_fp))

end_time <- Sys.time()
conduitR::log_with_timestamp("Completed downloading fasta files. Time taken: %.2f minutes", 
    as.numeric(difftime(end_time, start_time, units = "mins")))

# closing logfile connection
conduitR::log_with_timestamp("get_fasta.R script complete.")
# closing logfile connection
sink(type = "message")
sink()
close(zz)