################################################################################
# Preparing GO Annotations 
################################################################################
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

start_time <- Sys.time()

conduitR::log_with_timestamp("Running get_go_info.R script")

conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[["uniprot_annotated_protein_info"]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[["go_info"]]))

# Defining Input Files
uniprot_annotated_protein_info_fp = snakemake@input[["uniprot_annotated_protein_info"]]

# Defining Output Files
go_info_fp = snakemake@output[["go_info"]]

uniprot_annotated_protein_info <- readr::read_delim(uniprot_annotated_protein_info_fp)

go_long <- uniprot_annotated_protein_info |>
 dplyr::select(protein_id,go)|>
 # Extract description and GO ID together
  dplyr::mutate(matches = stringr::str_match_all(go, "\\s*([^;\\[]+?)\\s*\\[(GO:\\d{7})\\]")) |>
  tidyr::unnest(cols = c(matches)) |>
  # matches matrix: [,2] = description, [,3] = GO ID
  dplyr::mutate(
    description = stringr::str_trim(matches[,2]),  # remove any leading/trailing whitespace
    go_id = matches[,3],
    annotation_type ="go"
  ) |>
  dplyr::select(protein_id,annotation_type,term = go_id, description)|>
  dplyr::distinct()

conduitR::log_with_timestamp("Writing go annotation info to file")
readr::write_delim(go_long,go_info_fp)

end_time <- Sys.time()
conduitR::log_with_timestamp("Completed get_go_info.R script. Time taken: %.2f minutes", 
    as.numeric(difftime(end_time, start_time, units = "mins")))
# closing clogfile connection
sink(type = "message")
sink()
close(zz)