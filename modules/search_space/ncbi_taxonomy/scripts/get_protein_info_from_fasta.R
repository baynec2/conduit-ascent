################################################################################
# Extract Protein Info From Fasta
################################################################################
## Opening Log File 
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

# Get input and output files from Snakemake workflow
database_fasta <- snakemake@input[["database_fasta"]]
taxonomy_txt <- snakemake@input[["taxonomy_txt"]]
output_file <- snakemake@output[[1]]

start_time <- Sys.time()

# Now everything from print(), message(), warning() will go into the log file
conduitR::log_with_timestamp("Running get_protein_info_from_fasta.R script")
conduitR::log_with_timestamp(paste0("Input files: ", snakemake@input[[1]], " and ", snakemake@input[[2]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[[1]]))

# Loading organism annotation
conduitR::log_with_timestamp("Loading organism annotation from the taxonomy file.")
taxonomy_txt <- readr::read_delim(taxonomy_txt)

conduitR::log_with_timestamp("Extracting protein information from fasta database.")
# Extracting protein information from fasta database
protein_info <- conduitR::extract_fasta_info(database_fasta)

conduitR::log_with_timestamp("Adding organism annotation to protein info.")

# Adding organism annotation to protein info
protein_info <- dplyr::inner_join(protein_info, taxonomy_txt,
  by = "organism_id"
)

conduitR::log_with_timestamp("Writing protein info to file.")

readr::write_delim(protein_info, output_file, delim = "\t")

end_time <- Sys.time()
conduitR::log_with_timestamp("Completed get_protein_info_from_fasta.R script. Time taken: %.2f minutes", 
    as.numeric(difftime(end_time, start_time, units = "mins")))
# closing clogfile connection
sink(type = "message")
sink()
close(zz)
