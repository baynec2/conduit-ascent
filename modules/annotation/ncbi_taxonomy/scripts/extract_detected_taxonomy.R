################################################################################
# Extracting Detected Taxonomy
################################################################################
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

start_time <- Sys.time()

conduitR::log_with_timestamp("Running extract_detected_taxonomy.R script")
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[1]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[[1]]))

# Loading files
uniprot_annotated_protein_info = snakemake@input[["uniprot_annotated_protein_info"]]
detected_taxonomy =snakemake@output[["detected_taxonomy"]]

conduitR::log_with_timestamp(paste0("Reading in annotated protein info from ", uniprot_annotated_protein_info))
# Select only relevant columns
uniprot_annotated_protein_info = readr::read_delim(uniprot_annotated_protein_info) |>
  dplyr::select(protein_id,
                organism_type,
                domain,
                kingdom,
                phylum,
                class,
                order,
                family,
                genus,
                species)

conduitR::log_with_timestamp(paste0("Writing detected taxonomy to ", detected_taxonomy))
readr::write_delim(uniprot_annotated_protein_info,detected_taxonomy)

end_time <- Sys.time()
conduitR::log_with_timestamp("Completed extract_detected_taxonomy.R script. Time taken: %.2f minutes", 
    as.numeric(difftime(end_time, start_time, units = "mins")))

# closing clogfile connection
sink(type = "message")
sink()
close(zz)