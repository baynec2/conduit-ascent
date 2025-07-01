################################################################################
#Extracting GO Information from Uniprot Annotations
################################################################################
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

start_time <- Sys.time()

conduitR::log_with_timestamp("Running extract_go_info.R script")
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[1]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[[1]]))

uniprot_annotated_protein_info = snakemake@input[["uniprot_annotated_protein_info"]]
go_annotations_out = snakemake@output[["go_annotations"]]

conduitR::log_with_timestamp(paste0("Reading in annotated protein info from ", uniprot_annotated_protein_info))
# Loading annotated protein info
uniprot_annotated_protein_info = readr::read_delim(uniprot_annotated_protein_info)
# Defining columns in output
shared_columns = c("protein_id","organism_type","domain","kingdom",
                   "phylum","class","order","family","genus","species")

conduitR::log_with_timestamp("Extracting GO information from Uniprot annotations")
# Saving GO data from the annotations
go_annotations = uniprot_annotated_protein_info |>
  dplyr::select("go",dplyr::all_of(shared_columns)) |>
  # Split by semicolon and unnest into separate rows
  dplyr::mutate(go = stringr::str_split(go, "; ")) |>
  tidyr::unnest(go) |>
  # Extract description and GO term
  dplyr::mutate(
    description = stringr::str_extract(go, "^[^\\[]+"), # Extract the part before the bracket
    go_term = stringr::str_extract(go, "\\[GO:[^\\]]+\\]") # Extract the GO term inside the brackets
  ) |>
  # Clean up the columns (remove extra spaces and brackets)
  dplyr::mutate(
    description = stringr::str_trim(description),
    go_term = stringr::str_replace_all(go_term, "\\[|\\]", "")# Remove square brackets from GO term
  ) |>
  # Select the relevant columns
  dplyr::select(description, go_term,protein_id,dplyr::all_of(shared_columns))

conduitR::log_with_timestamp(paste0("Writing GO annotations to ", go_annotations_out))
# Write to file.
readr::write_delim(go_annotations,go_annotations_out,delim = "\t")

end_time <- Sys.time()
conduitR::log_with_timestamp("Completed extract_go_info.R script. Time taken: %.2f minutes", 
    as.numeric(difftime(end_time, start_time, units = "mins")))

# closing clogfile connection
sink(type = "message")
sink()
close(zz)