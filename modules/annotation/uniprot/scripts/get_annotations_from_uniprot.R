################################################################################
# Get Detected Proteins Annotation From Uniprot
################################################################################
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

start_time <- Sys.time()

conduitR::log_with_timestamp("Running get_annotations_from_uniprot.R script")
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[1]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[[1]]))

# Get input and output files from Snakemake workflow
detected_protein_info <- snakemake@input[["detected_protein_info"]]
uniprot_annotated_protein_info <- snakemake@output[["uniprot_annotated_protein_info"]]

conduitR::log_with_timestamp(paste0("Reading in detected protein information from ", detected_protein_info))
detected_protein_info = readr::read_delim(detected_protein_info) 

detected_protein_ids = detected_protein_info |>
  dplyr::pull("protein_id")

conduitR::log_with_timestamp("Getting annotations from Uniprot")
# Annotating with the relevant information.
uniprot_annotations <- conduitR::get_annotations_from_uniprot(detected_protein_ids)

conduitR::log_with_timestamp("Appending uniprot annotations to protein info")
# Appending uniprot annotations to protein info
protein_info <- dplyr::inner_join(detected_protein_info, uniprot_annotations,
                                 by = c("protein_id" = "accession"))

conduitR::log_with_timestamp(paste0("Writing to file ", uniprot_annotated_protein_info))
# Writing to file
readr::write_delim(protein_info,uniprot_annotated_protein_info)

end_time <- Sys.time()
conduitR::log_with_timestamp("Completed get_annotations_from_uniprot.R script. Time taken: %.2f minutes", 
    as.numeric(difftime(end_time, start_time, units = "mins")))

# closing clogfile connection
sink(type = "message")
sink()
close(zz)
