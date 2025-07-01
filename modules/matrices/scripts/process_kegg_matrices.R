# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

start_time <- Sys.time()
conduitR::log_with_timestamp("Running process_kegg_matrices.R script")

conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[["kegg_annotations"]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[["pg_matrix"]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[["kegg_pathway_matrix"]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[["kegg_ko_matrix"]]))

################################################################################
# Defining Inputs and Outputs
################################################################################
# Inputs
input_kegg_annotations = readr::read_delim(snakemake@input[["kegg_annotations"]])
input_pg_matrix = readr::read_tsv(snakemake@input[["pg_matrix"]])

# Outputs
output_kegg_pathway_matrix = snakemake@output[["kegg_pathway_matrix"]]
output_kegg_ko_matrix = snakemake@output[["kegg_ko_matrix"]]

################################################################################
# Normalizing protein group matrix
################################################################################
conduitR::log_with_timestamp("Normalizing protein group matrix")
# Adjusting the matrix to be positive and median centered  (deals with sample loading fluctuations)
# Finding first and last column containing data
first_col = which(names(input_pg_matrix) == "First.Protein.Description") + 1
last_col = length(input_pg_matrix)

selected_columns = names(input_pg_matrix)[first_col:last_col]

pg_matrix_mod <- input_pg_matrix |>
  dplyr::mutate(
    # Median centering to adjust for sample loading.
    dplyr::across(dplyr::all_of(selected_columns),
                  ~ .x - median(.x, na.rm = TRUE)),
    # Adjusting negative values to be positive. Offsetting by smallest value.
    dplyr::across(dplyr::all_of(selected_columns),
                  ~ .x + abs(min(.x, na.rm = TRUE)) + 1)  # Adjusting to make all values positive
  ) |>
  dplyr::mutate(master_protein_id = gsub(";.*","",Protein.Group),.before = "Protein.Group")

# Combining KEGG annotations with protein group matrix
combined = dplyr::inner_join(input_kegg_annotations,pg_matrix_mod,by = c("protein_id" = "master_protein_id"))
################################################################################
# Processing Pathway Matrix
################################################################################
conduitR::log_with_timestamp("Processing KEGG Pathway Matrix")
kegg_pathway_matrix = combined |>
  dplyr::mutate(kegg_pathway_id = paste0(organism_type,"_",kegg_pathway)) |>
  dplyr::group_by(kegg_pathway_id,organism_type,kegg_pathway) |>
  dplyr::summarise(dplyr::across(dplyr::all_of(selected_columns),
                   ~ sum(.x, na.rm = TRUE)),
                   .groups = "drop")

################################################################################
# Processing KO Matrix
################################################################################
conduitR::log_with_timestamp("Processing KEGG KO Matrix")

ko_matrix = combined |>
  dplyr::mutate(ko_id = paste0(organism_type,"_",ko)) |>
  dplyr::group_by(ko_id,organism_type,ko,ko_description) |>
  dplyr::summarise(dplyr::across(dplyr::all_of(selected_columns),
                          ~ sum(.x, na.rm = TRUE)),
                          .groups = "drop")
                          
################################################################################
# Saving Outputs
################################################################################
readr::write_tsv(kegg_pathway_matrix,output_kegg_pathway_matrix)
readr::write_tsv(ko_matrix,output_kegg_ko_matrix)

end_time <- Sys.time()
conduitR::log_with_timestamp("Completed process_kegg_matrices.R script. Time taken: %.2f minutes", 
    as.numeric(difftime(end_time, start_time, units = "mins")))

# closing clogfile connection
sink(type = "message")
sink()
close(zz)