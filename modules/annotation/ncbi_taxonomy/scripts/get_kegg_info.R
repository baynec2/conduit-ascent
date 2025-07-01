################################################################################
# Getting KEGG Information for each KEGG ID
################################################################################
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

start_time <- Sys.time()
conduitR::log_with_timestamp("Running get_kegg_info.R script")
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[1]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[[1]]))


# Loading files
uniprot_annotated_protein_info = snakemake@input[["uniprot_annotated_protein_info"]]
kegg_annotations_out =snakemake@output[[1]]

conduitR::log_with_timestamp(paste0("Reading in detected protein info from ", uniprot_annotated_protein_info))
# Reading in detected protein info
uniprot_annotated_protein_info = readr::read_delim(uniprot_annotated_protein_info)

# Defining columns that all outputs will share.
shared_columns = c("protein_id","organism_type","domain","kingdom",
                   "phylum","class","order","family","genus","species")

# Getting KEGG Information for each KEGG ID

# Pulling out the KEGG ids
kegg_ids = uniprot_annotated_protein_info |>
  dplyr::filter(!is.na("xref_kegg")) |>
  dplyr::pull("xref_kegg")

conduitR::log_with_timestamp(paste0("Retrieving KEGG information for ", length(kegg_ids)," KEGG IDs"))
# Getting the terms for each annotation
kegg_annotations = conduitR::get_kegg_in_batches(kegg_ids)

# Joining back to original annotation file
kegg_annotations_final = uniprot_annotated_protein_info |>
  tidyr::separate_rows(sep =";",xref_kegg) |>
  dplyr::left_join(kegg_annotations,by = c("xref_kegg" = "kegg_id")) |>
  dplyr::select(kegg_pathway,kegg_id = xref_kegg,ko,ko_description,
                dplyr::all_of(shared_columns)) |>
  dplyr::filter(kegg_id != "")

conduitR::log_with_timestamp(paste0("Writing KEGG annotations to ", kegg_annotations_out))
# Writing to File
readr::write_delim(kegg_annotations_final,kegg_annotations_out,delim = "\t")

end_time <- Sys.time()
conduitR::log_with_timestamp("Completed get_kegg_info.R script. Time taken: %.2f minutes", 
    as.numeric(difftime(end_time, start_time, units = "mins")))

# closing clogfile connection
sink(type = "message")
sink()
close(zz)
