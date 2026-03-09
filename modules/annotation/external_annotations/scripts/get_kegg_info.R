################################################################################
# Getting Kegg Annotations for each Detected Protein ID
################################################################################
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

start_time <- Sys.time()

conduitR::log_with_timestamp("Running get_kegg_info.R script")
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[["uniprot_annotated_protein_info"]]))
conduitR::log_with_timestamp(paste0("Output file: ",snakemake@output[["kegg_pathway_info"]]))
conduitR::log_with_timestamp(paste0("Output file: ",snakemake@output[["kegg_map_pathway_info"]]))
conduitR::log_with_timestamp(paste0("Output file: ",snakemake@output[["kegg_orthology_info"]]))
#Defining files
# Inputs
uniprot_annotated_protein_info_fp = snakemake@input[["uniprot_annotated_protein_info"]]

# Outputs
kegg_pathway_info_fp = snakemake@output[["kegg_pathway_info"]]
kegg_map_pathway_info_fp = snakemake@output[["kegg_map_pathway_info"]]
kegg_orthology_info_fp = snakemake@output[["kegg_orthology_info"]]


# Reading in files
conduitR::log_with_timestamp("Reading in files")
uniprot_annotated_protein_info = readr::read_delim(uniprot_annotated_protein_info_fp)

# Converting Annotations to long format 
conduitR::log_with_timestamp("Converting uniprot annotated protein info to long format.")  


# Converting KEGG Terms to Long Format
conduitR::log_with_timestamp("Processing KEGG xrefs")  

kegg_long <- uniprot_annotated_protein_info |>
    dplyr::select(protein_id,
    xref_kegg) |>
  # remove trailing semicolon if present
  dplyr::mutate(xref_kegg = stringr::str_remove(xref_kegg, ";$")) |>
  # separate into long format, keeping IDs and descriptions aligned
  tidyr::separate_rows(xref_kegg, sep = ";")|>
  dplyr::filter(!is.na(xref_kegg))

# Getting KEGG pathway information from IDs.
conduitR::log_with_timestamp("Getting data from Kegg database")  
unique_kegg = unique(kegg_long$xref_kegg)
kegg_database = conduitR::get_kegg_in_batches(unique_kegg)

kegg_combined = kegg_long |>
dplyr::left_join(kegg_database,by = c("xref_kegg" = "kegg_id"))

# Kegg Pathway 
conduitR::log_with_timestamp("Extracting Kegg Pathway Info")  

kegg_pathway_long = kegg_combined |>
dplyr::select(protein_id,kegg_pathway_id,kegg_pathway)|>
dplyr::mutate(annotation_type = "kegg_pathway")|>
dplyr::select(protein_id,annotation_type,term = kegg_pathway_id, description = kegg_pathway)|>
dplyr::filter(!is.na(term))

# Generic Kegg Pathways
conduitR::log_with_timestamp("Converting Species Specific Kegg Pathways to General Kegg Pathways")  

kegg_map_pathway_long = kegg_combined |>
dplyr::select(protein_id,kegg_pathway_id,kegg_pathway)|>
dplyr::mutate(kegg_map_pathway_id = gsub("^[a-z]{3}","map",kegg_pathway_id))|>
dplyr::mutate(annotation_type = "kegg_map_pathway")|>
dplyr::select(protein_id,annotation_type,term = kegg_map_pathway_id, description = kegg_pathway)|>
dplyr::filter(!is.na(term))

# Kegg Orthology
conduitR::log_with_timestamp("Extracting Kegg orthology info")  

ko_long = kegg_combined |>
dplyr::mutate(annotation_type = "kegg_orthology")|>
dplyr::select(protein_id,annotation_type,term = ko,description =ko_description)|>
dplyr::distinct()|>
dplyr::filter(!is.na(term))

conduitR::log_with_timestamp("Saving Kegg Pathway Annotations")
readr::write_delim(kegg_pathway_long,kegg_pathway_info_fp)

conduitR::log_with_timestamp("Saving Kegg Pathway Map Annotations")
readr::write_delim(kegg_map_pathway_long,kegg_map_pathway_info_fp)

conduitR::log_with_timestamp("Saving Kegg Orthology Annotations")
readr::write_delim(ko_long,kegg_orthology_info_fp)

end_time <- Sys.time()
conduitR::log_with_timestamp("Completed get_kegg_info.R script. Time taken: %.2f minutes", 
    as.numeric(difftime(end_time, start_time, units = "mins")))
# closing clogfile connection
sink(type = "message")
sink()
close(zz)

