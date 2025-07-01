################################################################################
# Extract Detected Proteins
################################################################################
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

# Get input and output files from Snakemake workflow
protein_info_df <- snakemake@input[["protein_info_df"]]
protein_info_fasta <- snakemake@input[["protein_info_fasta"]]
report_pg_matrix <- snakemake@input[["report_pg_matrix"]]
detected_protein_info_df <- snakemake@output[["detected_protein_info_df"]]
detected_protein_info_fasta<- snakemake@output[["detected_protein_info_fasta"]]

start_time <- Sys.time()

conduitR::log_with_timestamp("Running extract_detected_proteins.R script")
conduitR::log_with_timestamp(paste0("Input file: ", protein_info_df))
conduitR::log_with_timestamp(paste0("Input file: ", protein_info_fasta))
conduitR::log_with_timestamp(paste0("Input file: ", report_pg_matrix))
conduitR::log_with_timestamp(paste0("Output file: ", detected_protein_info_df))
conduitR::log_with_timestamp(paste0("Output file: ", detected_protein_info_fasta))

# Make the directory if it doesn't exist
if(!dir.exists(dirname(detected_protein_info_df))){
dir.create(dirname(detected_protein_info_df))
}
conduitR::log_with_timestamp(paste0("Reading in protein information from ", protein_info_df))
# Reading in protein information
protein_info = readr::read_tsv(protein_info_df)
conduitR::log_with_timestamp(paste0("Reading in report pg matrix from ", report_pg_matrix))
# Determining what the detected proteins from DIA-NN were detected
report_pg_matrix = readr::read_tsv(report_pg_matrix)
detected_uniprot_ids = report_pg_matrix |>
  dplyr::pull("Protein.Group") |>
  strsplit(";") |>
  unlist() |>
  unique()

conduitR::log_with_timestamp("Filtering detected protein information")
# Filtering detected protein information
detected_protein_info = protein_info |>
  dplyr::filter(protein_id %in% detected_uniprot_ids)

conduitR::log_with_timestamp(paste0("Writing detected protein information to ", detected_protein_info_df))
# Writing to file - data frame
readr::write_tsv(detected_protein_info,detected_protein_info_df)

### Filtering Fasta file to only include detected proteins ###
conduitR::log_with_timestamp("Filtering fasta file to only include detected proteins")
# Reading in fasta file
sequences <- Biostrings::readAAStringSet(protein_info_fasta)
# Extract UniProt IDs from FASTA headers
uniprot_ids <- names(sequences) |> stringr::str_extract("(?<=\\|)[A-Z0-9]+(?=\\|)")
# Create a named vector mapping UniProt IDs to full headers
header_map <- setNames(names(sequences), uniprot_ids)

# Filtering the header map to include only the ids detected in study.
detected_headers = header_map[names(header_map) %in% detected_protein_info$protein_id]

# Producing a filtered_fasta
filtered_fasta <- sequences[names(sequences) %in% detected_headers]

conduitR::log_with_timestamp(paste0("Writing filtered fasta file to ", detected_protein_info_fasta))

# Writing Filtered Fasta file
Biostrings::writeXStringSet(filtered_fasta,detected_protein_info_fasta)

end_time <- Sys.time()

conduitR::log_with_timestamp("Completed extract_detected_proteins.R script. Time taken: %.2f minutes", 
    as.numeric(difftime(end_time, start_time, units = "mins")))

# closing clogfile connection
sink(type = "message")
sink()
close(zz)
