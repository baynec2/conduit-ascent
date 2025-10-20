################################################################################
# Preparing Annotations for Conduit Class
################################################################################
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

start_time <- Sys.time()

conduitR::log_with_timestamp("Running prepare_annotations.R script")

conduitR::log_with_timestamp(paste0("Input file: ", snakemake@output[["annotated_qf"]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[["uniprot_annotated_protein_info"]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[["conduit_annotations"]]))

#Defining files
# Inputs
annotated_qf_fp = snakemake@input[["annotated_qf"]]
uniprot_annotated_protein_info_fp = snakemake@input[["uniprot_annotated_protein_info"]]

# Outputs
conduit_annotations_fp = snakemake@output[["conduit_annotations"]]

# Reading in files
conduitR::log_with_timestamp("Reading in files")
uniprot_annotated_protein_info = readr::read_delim(uniprot_annotated_protein_info_fp)
annotated_qf = readRDS(annotated_qf_fp)

# Converting Annotations to long format 
conduitR::log_with_timestamp("Converting uniprot annotated protein info to long format.")  

sel = uniprot_annotated_protein_info |>
        dplyr::select(protein_id,
                      go,xref_kegg)

# Converting GO Terms to Long Format
conduitR::log_with_timestamp("Processing GO terms")  

go_long <- sel |>
  dplyr::select(protein_id, go) |>
  # Extract description and GO ID together
  dplyr::mutate(matches = stringr::str_match_all(go, "\\s*([^;\\[]+?)\\s*\\[(GO:\\d{7})\\]")) |>
  tidyr::unnest(cols = c(matches)) |>
  # matches matrix: [,2] = description, [,3] = GO ID
  dplyr::mutate(
    description = stringr::str_trim(matches[,2]),  # remove any leading/trailing whitespace
    go_id = matches[,3],
    annotation_type ="go"
  ) |>
  dplyr::select(protein_id,annotation_type,term = go_id, description)

# Converting KEGG Terms to Long Format
conduitR::log_with_timestamp("Processing KEGG terms")  

kegg_long <- sel |>
  dplyr::select(protein_id, xref_kegg) |>
  # remove trailing semicolon if present
  dplyr::mutate(xref_kegg = stringr::str_remove(xref_kegg, ";$")) |>
  # separate into long format, keeping IDs and descriptions aligned
  tidyr::separate_rows(xref_kegg, sep = ";")|>
  dplyr::select(protein_id,xref_kegg) |>
  dplyr::filter(!is.na(xref_kegg))

# Getting KEGG pathway information from IDs.
conduitR::log_with_timestamp("Getting data from Kegg database")  
unique_kegg = unique(kegg_long$xref_kegg)
kegg_database = conduitR::get_kegg_in_batches(unique_kegg)

kegg_combined = kegg_long |>
dplyr::left_join(kegg_database,by = c("xref_kegg" = "kegg_id"))

kegg_pathway_long = sel |>
dplyr::select(-xref_kegg)|>
dplyr::left_join(kegg_combined, by = "protein_id") |>
dplyr::mutate(annotation_type = "kegg_pathway")|>
dplyr::select(protein_id,annotation_type,term = kegg_pathway_id, description = kegg_pathway)

conduitR::log_with_timestamp("Combining annotations")
combined_annotations = dplyr::bind_rows(go_long,kegg_pathway_long)

conduitR::log_with_timestamp("Converting protein IDs to protein groups")
annotated_qf = readRDS(annotated_qf_fp)
taxonomy_columns = c("domain","kingdom","phylum","class","order","family","genus","species","lca")
protein_groups <- SummarizedExperiment::rowData(annotated_qf[["protein_groups"]])[, c("Protein.Group",
taxonomy_columns)] |>
  tibble::as_tibble() |>
  dplyr::rename(Protein.Group.Temp = Protein.Group) |>
  dplyr::mutate(Protein.Group = Protein.Group.Temp)|>
  tidyr::separate_rows(Protein.Group.Temp, sep = ";") |>
  dplyr::select(Protein.Group,protein_id = Protein.Group.Temp,dplyr::everything())

conduit_annotations = protein_groups |>
  dplyr::left_join(combined_annotations,protein_groups, by = "protein_id")|>
  dplyr::select(Protein.Group,protein_id,species,lca,annotation_type,term,description)|>
  dplyr::distinct()|>
  dplyr::filter(!is.na(term))

conduitR::log_with_timestamp("Writing conduit annotations to file")
readr::write_tsv(conduit_annotations,conduit_annotations_fp)

end_time <- Sys.time()
conduitR::log_with_timestamp("Completed prepare_annotations.R script. Time taken: %.2f minutes", 
    as.numeric(difftime(end_time, start_time, units = "mins")))
# closing clogfile connection
sink(type = "message")
sink()
close(zz)
