# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages
start_time <- Sys.time()

conduitR::log_with_timestamp("Running merge_profiles.R script")

# Defining inputs
metaphlan_profiles_fp = snakemake@input[["metaphlan_profiles"]]

# Defining outputs
merged_profiles_fp = snakemake@output[["merged_profiles"]]

conduitR::log_with_timestamp("Reading in metaphlan files")

# 2. Function to read a single MetaPhlAn profile
read_metaphlan <- function(file) {
  sample_name <- base::gsub("_profile\\.txt$", "", base::basename(file))
  readr::read_tsv(file,
     comment = "#",
     col_names = c("clade_name", "NCBI_tax_id","relative_abundance","additional_species"), 
     col_types = readr::cols(clade_name = readr::col_character(),
                             NCBI_tax_id = readr::col_character(),
                             relative_abundance = readr::col_double(),
                             additional_species = readr::col_character())) |>
    dplyr::mutate(sample = sample_name)
}

# 3. Read all files and combine
all_profiles <- purrr::map_dfr(metaphlan_profiles_fp, read_metaphlan)

# 4. Separate NCBI_tax_id if there are multiple IDs (like 2|976)
conduitR::log_with_timestamp("Converting to wide format")

# 5. Pivot so that each sample is a column (wide format)
merged_profiles <- tidyr::pivot_wider(
  all_profiles |>
    dplyr::select(clade_name, NCBI_tax_id, relative_abundance, sample),
  names_from = sample,
  values_from = relative_abundance,
  values_fill = 0
)
conduitR::log_with_timestamp("Writing merged results to file")

readr::write_delim(merged_profiles,merged_profiles_fp)

end_time <- Sys.time()
conduitR::log_with_timestamp("merge_profiles.R script completed. Time taken: %.2f minutes", 
                              as.numeric(difftime(end_time, start_time, units = "mins")))
# closing clogfile connection
sink(type = "message")
sink()
close(zz)