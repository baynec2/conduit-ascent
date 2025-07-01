################################################################################
#Extracting Subcellular Location Annotations from Uniprot Data
################################################################################
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

conduitR::log_with_timestamp("Running extract_subcellular_locations.R script")
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[1]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[[1]]))

uniprot_annotated_protein_info = snakemake@input[["uniprot_annotated_protein_info"]]
subcellular_locations_out = snakemake@output[["subcellular_locations"]]

# Getting Cellular Location Information
# Defining columns in output
shared_columns = c("protein_id","organism_type","domain","kingdom",
                   "phylum","class","order","family","genus","species")

conduitR::log_with_timestamp(paste0("Reading in annotated protein info from ", uniprot_annotated_protein_info))

uniprot_annotated_protein_info =readr::read_delim(uniprot_annotated_protein_info)

conduitR::log_with_timestamp("Extracting subcellular location information from Uniprot annotations")
subcellular_locations = uniprot_annotated_protein_info |>
  dplyr::select("cc_subcellular_location",dplyr::all_of(shared_columns)) |>
  dplyr::mutate(cc_subcellular_location = stringr::str_split(cc_subcellular_location, "; ")) |>  # Split into list
  tidyr::unnest(cc_subcellular_location) |>  # Expand to long format
  dplyr::mutate(
    subcellular_location = stringr::str_extract(cc_subcellular_location, "^[^{]+"),  # Extract location (before {)
    evidence = stringr::str_extract(cc_subcellular_location, "\\{([^}]+)\\}")  # Extract text inside curly braces
  ) |>
  dplyr::mutate(
    evidence = stringr::str_remove_all(evidence, "[{}]"),  # Clean up curly braces
    cc_subcellular_location = gsub("SUBCELLULAR LOCATION: ", "",
                                   cc_subcellular_location),
    cc_subcellular_location = gsub("\\{.*","",cc_subcellular_location),
    cc_subcellular_location = gsub("\\.","",cc_subcellular_location),
    note = stringr::str_extract(cc_subcellular_location, "Note.*"),
    cc_subcellular_location = gsub("Note.*","",cc_subcellular_location)) |>
  dplyr::select(cc_subcellular_location, evidence,note, protein_id,
                dplyr::all_of(shared_columns))  # Keep relevant columns

conduitR::log_with_timestamp(paste0("Writing subcellular locations to ", subcellular_locations_out))
# Writing to file.
readr::write_delim(subcellular_locations,
                   subcellular_locations_out,
                   delim = "\t")
conduitR::log_with_timestamp("extract_subcellular_locations.R script completed")
# closing clogfile connection
sink(type = "message")
sink()
close(zz)
