################################################################################
# Extracting Metrics
################################################################################
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages
start_time <- Sys.time()

conduitR::log_with_timestamp("Running extract_metrics.R script")
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[1]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[2]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[[1]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[[2]]))

# Inputs
protein_info = snakemake@input[["protein_info"]]
detected_protein_info = snakemake@input[["detected_protein_info"]]

# Defining outputs
database_taxonomy=snakemake@output[["database_taxonomy"]]
database_metrics=snakemake@output[["database_metrics"]]
detected_protein_taxonomy = snakemake@output[["detected_protein_taxonomy"]]
detected_protein_metrics = snakemake@output[["detected_protein_metrics"]]
combined_metrics = snakemake@output[["combined_metrics"]]

# Database information
conduitR::log_with_timestamp(paste0("Reading in protein info from ", protein_info))
protein_info = readr::read_tsv(protein_info)

conduitR::log_with_timestamp("Generating database taxonomy")
# Generating database_taxonomy
database_taxonomy_out = protein_info |>
  dplyr::select(protein_id,organism_type,domain,kingdom,phylum,class,
                order,family,genus,species)

conduitR::log_with_timestamp(paste0("Writing database taxonomy to ", database_taxonomy))
# Writing to file
readr::write_tsv(database_taxonomy_out,database_taxonomy)

conduitR::log_with_timestamp("Generating database metrics")
# Generating metrics
database_metrics_out= protein_info |>
  dplyr::select(-organism_name,-organism_id,-sequence) |>
  tidyr::pivot_longer(1:10,names_to = "metric") |>
  dplyr::group_by(metric) |>
  dplyr::summarise(n_in_db = dplyr::n_distinct(value)) |>
  dplyr::ungroup() |>
  dplyr::arrange(factor(metric, levels = c("organism_type","domain","kingdom","phylum","class","order","family","genus","species","protein_id")))

conduitR::log_with_timestamp(paste0("Writing database metrics to ", database_metrics))
# Writing to file
readr::write_tsv(database_metrics_out,database_metrics)

# Detected protein Information
conduitR::log_with_timestamp(paste0("Reading in detected protein info from ", detected_protein_info))
detected_protein_info = readr::read_tsv(detected_protein_info)

conduitR::log_with_timestamp("Generating detected protein taxonomy")
# Generating database_taxonomy
detected_protein_taxonomy_out = detected_protein_info |>
  dplyr::select(protein_id,organism_type,domain,kingdom,phylum,class,
                order,family,genus,species)

conduitR::log_with_timestamp(paste0("Writing detected protein taxonomy to ", detected_protein_taxonomy))
# Writing to file
readr::write_tsv(detected_protein_taxonomy_out,detected_protein_taxonomy)

# Generating metrics
detected_protein_metrics_out= detected_protein_info |>
  dplyr::select(-organism_name,-organism_id,-sequence) |>
  tidyr::pivot_longer(1:10,names_to = "metric") |>
  dplyr::group_by(metric) |>
  dplyr::summarise(n_detected = dplyr::n_distinct(value)) |>
  dplyr::ungroup() |>
  dplyr::arrange(factor(metric, levels = c("organism_type","domain","kingdom","phylum","class","order","family","genus","species","protein_id")))

# Writing to file
readr::write_tsv(detected_protein_metrics_out,detected_protein_metrics)

conduitR::log_with_timestamp("Combining metrics into one tsv file")
# Combining metrics into one file.
combined_metrics_out = dplyr::inner_join(database_metrics_out,detected_protein_metrics_out,
                            by = "metric") |>
  dplyr::mutate(per_detected = round((n_detected/n_in_db) * 100)) |>
  dplyr::ungroup() |>
  dplyr::arrange(factor(metric, levels = c("organism_type","domain","kingdom","phylum","class","order","family","genus","species","protein_id")))

conduitR::log_with_timestamp(paste0("Writing combined metrics to ", combined_metrics) )
readr::write_tsv(combined_metrics_out,combined_metrics)
end_time <- Sys.time()  
conduitR::log_with_timestamp("extract_metrics.R script completed. Time taken: %.2f minutes", 
                              as.numeric(difftime(end_time, start_time, units = "mins")))
# closing clogfile connection
sink(type = "message")
sink()
close(zz)
