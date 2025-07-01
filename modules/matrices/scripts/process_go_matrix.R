################################################################################
# Processing GO Matrices
################################################################################
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

start_time <- Sys.time()

conduitR::log_with_timestamp("Running process_go_matrix.R script")
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[1]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[2]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[[1]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[[2]]))

# Getting inputs and outputs from snakemake.
report_pg_matrix = snakemake@input[["report_pg_matrix"]]
go_annotations = snakemake@input[["go_annotations"]]
go_annotations_matrix_fp = snakemake@output[["go_annotations_matrix"]]
go_annotations_matrix_taxa_fp = snakemake@output[["go_annotations_taxa_matrix"]]

conduitR::log_with_timestamp(paste0("Reading in report_pg_matrix from ", report_pg_matrix))
report_pg_matrix = readr::read_tsv(report_pg_matrix)

conduitR::log_with_timestamp(paste0("Reading in go_annotations from ", go_annotations))
go_annotations = readr::read_delim(go_annotations)

# Finding first and last column containing data
first_col = which(names(report_pg_matrix) == "First.Protein.Description") + 1
last_col = length(report_pg_matrix)

selected_columns = names(report_pg_matrix)[first_col:last_col]

# Adjusting the matrix to be positive and median centered  (deals with sample loading fluctuations)
report_pg_matrix_mod <- report_pg_matrix |>
  dplyr::mutate(
    # Median centering to adjust for sample loading.
    dplyr::across(first_col:last_col,
                  ~ .x - median(.x, na.rm = TRUE)),
    # Adjusting negative values to be positive. Offsetting by smallest value.
    dplyr::across(first_col:last_col,
                  ~ .x + abs(min(.x, na.rm = TRUE)) + 1)  # Adjusting to make all values positive
  ) |>
  dplyr::mutate(master_protein_id = gsub(";.*","",Protein.Group),.before = "Protein.Group")

combined = dplyr::inner_join(go_annotations,
                             report_pg_matrix_mod,
                             by = c("protein_id" = "master_protein_id"))

go_annotations_matrix = combined |>
  dplyr::mutate(go_id = paste0(organism_type,"_",description)) |>
  dplyr::group_by(go_id,organism_type,description) |>
  dplyr::summarise(dplyr::across(dplyr::all_of(selected_columns),
                   ~ sum(.x, na.rm = TRUE))
  )

conduitR::log_with_timestamp(paste0("Writing go_annotations_matrix to ", go_annotations_matrix_fp))
# Writing to file
readr::write_tsv(go_annotations_matrix,go_annotations_matrix_fp)

conduitR::log_with_timestamp("Processing taxonomic go_annotations_matrix")
go_annotations_matrix_taxa = combined |>
  dplyr::mutate(go_id = paste0(organism_type,"_", domain,"_",kingdom,"_",
                               phylum,"_",class,"_",order,"_",family,"_",genus,
                               "_",species,"_",description),.before = 1

  ) |>
  dplyr::group_by(go_id,organism_type,description,domain,kingdom,phylum,
                  class,order,family,genus,species) |>
  dplyr::summarise(dplyr::across(dplyr::all_of(selected_columns),
                                 ~ sum(.x, na.rm = TRUE)))

conduitR::log_with_timestamp(paste0("Writing taxonomic go_annotations_matrix to ", go_annotations_matrix_taxa_fp))
# Writing to file
readr::write_tsv(go_annotations_matrix_taxa,go_annotations_matrix_taxa_fp)

end_time <- Sys.time()
conduitR::log_with_timestamp("Completed process_go_matrix.R script. Time taken: %.2f minutes", 
    as.numeric(difftime(end_time, start_time, units = "mins")))

# closing clogfile connection
sink(type = "message")
sink()
close(zz)