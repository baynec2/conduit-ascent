################################################################################
# Processing Subcellular Locations Matrices
################################################################################
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

start_time <- Sys.time()

conduitR::log_with_timestamp("Running 02_process_subcellular_locations_matrix.R script")
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[1]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[2]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[[1]]))

# Loading from Snakemake
report_pg_matrix=snakemake@input[["report_pg_matrix"]]
subcellular_locations = snakemake@input[["subcellular_locations"]]
subcellular_locations_matrix = snakemake@output[["subcellular_locations_matrix"]]

# Reading files
conduitR::log_with_timestamp(paste0("Reading in report_pg_matrix from ", report_pg_matrix))
report_pg_matrix = readr::read_tsv(report_pg_matrix)
conduitR::log_with_timestamp(paste0("Reading in subcellular_locations from ", subcellular_locations))
subcellular_locations = readr::read_delim(subcellular_locations)

# Finding first and last column containing data
first_col = which(names(report_pg_matrix) == "First.Protein.Description") + 1
last_col = length(report_pg_matrix)

conduitR::log_with_timestamp("Processing subcellular locations matrix")

# Dividing each intensity by the number of proteins it maps to
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

combined = dplyr::inner_join(subcellular_locations,report_pg_matrix_mod,by = c("protein_id" = "master_protein_id")) |>
  dplyr::mutate(cc_subcellular_location = gsub("\\s+$","",cc_subcellular_location))

# Finding first and last column containing data
first_col = which(names(combined) == "First.Protein.Description") + 1
last_col = ncol(combined)

selected_columns <- names(combined)[first_col:last_col]

# Summing values across the selected columns for each group
sum = combined |>
  dplyr::mutate(subcellular_location_id = paste0(organism_type,"_",cc_subcellular_location),.before = 1) |>
  dplyr::group_by(subcellular_location_id,cc_subcellular_location,organism_type) |>
  dplyr::summarise(
    dplyr::across(dplyr::all_of(selected_columns),
                  ~ sum(.x, na.rm = TRUE))
  )

conduitR::log_with_timestamp(paste0("Writing subcellular locations matrix to ", subcellular_locations_matrix))
# Writing to file.
readr::write_tsv(sum,subcellular_locations_matrix)

end_time <- Sys.time()
conduitR::log_with_timestamp("Completed 02_process_subcellular_locations_matrix.R script. Time taken: %.2f minutes", 
    as.numeric(difftime(end_time, start_time, units = "mins")))

# closing clogfile connection
sink(type = "message")
sink()