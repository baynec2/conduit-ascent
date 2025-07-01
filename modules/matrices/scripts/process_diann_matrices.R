################################################################################
# Processing Diann Matrices
################################################################################
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

start_time <- Sys.time()

conduitR::log_with_timestamp("Running process_diann_matrices.R script")
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[1]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[2]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[[1]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[[2]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[[3]]))

# Defining files
# Inputs
report_pg_matrix = snakemake@input[["report_pg_matrix"]]
report_pr_matrix = snakemake@input[["report_pr_matrix"]]
# Outputs
protein_group_matrix = snakemake@output[["protein_group_matrix"]]
precursor_matrix = snakemake@output[["precursor_matrix"]]
peptide_matrix = snakemake@output[["peptide_matrix"]]

conduitR::log_with_timestamp(paste0("Reading in report_pr_matrix from ", report_pr_matrix))
conduitR::log_with_timestamp("Processing report_pr_matrix")

# Converting precursor format - ID needs to be first.
precursor_matrix_out = readr::read_tsv(report_pr_matrix) |>
  dplyr::select(Precursor.Id,dplyr::everything()) |>
  # Renaming the samples to not be the full path
  dplyr::rename_with(
    .fn = ~ ifelse(
      stringr::str_ends(., ".raw"),
      stringr::str_remove(stringr::str_extract(., "[^/]+(?=\\.raw$)"), "\\.raw$"),
      .
    )
  )

conduitR::log_with_timestamp(paste0("Writing precursor matrix to ", precursor_matrix))
readr::write_tsv(precursor_matrix_out,precursor_matrix)

# Creating peptide format - summing peptides with same sequence (including mods)
first_col = names(precursor_matrix_out)[which(names(precursor_matrix_out) == "Precursor.Charge") + 1]
last_col = names(precursor_matrix_out[ncol(precursor_matrix_out)])

conduitR::log_with_timestamp("Processing precursor matrix to peptide matrix")
peptide_matrix_out = precursor_matrix_out |>
  dplyr::group_by(Modified.Sequence,Protein.Group,Protein.Ids,Protein.Names,
                  First.Protein.Description,Stripped.Sequence) |>
  dplyr::summarise(dplyr::across(first_col:last_col,
                   ~ sum(.x, na.rm = TRUE))
  ) |>
  # Renaming the samples to not be the full path
  dplyr::rename_with(
    .fn = ~ ifelse(
      stringr::str_ends(., ".raw"),
      stringr::str_remove(stringr::str_extract(., "[^/]+(?=\\.raw$)"), "\\.raw$"),
      .
    )
  )
conduitR::log_with_timestamp(paste0("Writing peptide matrix to ", peptide_matrix))
readr::write_tsv(peptide_matrix_out,peptide_matrix)

# Protein Group format stays the same
protein_group_matrix_out = readr::read_tsv(report_pg_matrix) |>
  # Renaming the samples to not be the full path
  dplyr::rename_with(
    .fn = ~ ifelse(
      stringr::str_ends(., ".raw"),
      stringr::str_remove(stringr::str_extract(., "[^/]+(?=\\.raw$)"), "\\.raw$"),
      .
    )
  )

conduitR::log_with_timestamp(paste0("saving protein groupmatrix to ", protein_group_matrix))

readr::write_tsv(protein_group_matrix_out,protein_group_matrix)

end_time <- Sys.time()
conduitR::log_with_timestamp("Completed process_diann_matrices.R script. Time taken: %.2f minutes", 
    as.numeric(difftime(end_time, start_time, units = "mins")))

# closing clogfile connection
sink(type = "message")
sink()
close(zz)
