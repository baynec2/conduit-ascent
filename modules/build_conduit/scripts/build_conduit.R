################################################################################
# Building a Conduit Object 
################################################################################
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

start_time <- Sys.time()
conduitR::log_with_timestamp("Running build_conduit.R script")

conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[["diann_stats"]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[["qfeatures"]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[["database"]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[["annotations"]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[["taxonomy"]]))

conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[["conduit"]]))

#Defining files
# Inputs
diann_stats_fp = snakemake@input[["diann_stats"]]
QFeatures_fp = snakemake@input[["qfeatures"]]
database_fp = snakemake@input[["database"]]
annotations_fp = snakemake@input[["annotations"]]
taxonomy_fp = snakemake@input[["taxonomy"]]
# Params
workflow_version  <- snakemake@params[["workflow_version"]]
snakemake_version <- snakemake@params[["snakemake_version"]]
# Outputs
conduit_fp = snakemake@output[["conduit"]]
# Reading files
conduitR::log_with_timestamp("Reading in input files")
  # Reading in QFeatures Object
QFeatures <- readRDS(QFeatures_fp)
  
  # Reading in Metrics
diann_stats <- readr::read_tsv(diann_stats_fp) |>
    dplyr::mutate(File.Name = tools::file_path_sans_ext(basename(File.Name)))
  
metrics <- list(diann_stats = diann_stats)
  
  # Reading in Database (in tabular format)
database <- readr::read_tsv(database_fp) |>
    # Only keeping protein ids and corresponding organism ids
    dplyr::select(protein_id,organism_id) |>
    # Saving as factors to reduce memory footprint
    dplyr::mutate(dplyr::across(where(is.character), as.factor))

  # Reading in Annotation
annotations <- readr::read_delim(annotations_fp) |>
    dplyr::select(Protein.Group, annotation_type,
    term,description) |>
    # Saving as factors to reduce memory footprint
    dplyr::mutate(dplyr::across(where(is.character), as.factor))
  
taxonomy <- readr::read_delim(taxonomy_fp) |>
    # Saving as factors to reduce memory footprint
    dplyr::mutate(dplyr::across(where(is.character), as.factor))

conduitR::log_with_timestamp("Building provenance metadata")

parse_diann_cfg <- function(path) {
  lines <- readLines(path)
  lines <- lines[nzchar(trimws(lines)) & !startsWith(trimws(lines), "#")]
  tibble::tibble(parameter = paste0("line_", seq_along(lines)), value = lines)
}

config_list <- list(
  snakemake_yaml = tibble::tibble(
    parameter = names(snakemake@config),
    value     = vapply(snakemake@config, function(x) {
      if (is.null(x) || length(x) == 0) return("")
      if (is.list(x)) paste(names(x), unlist(x), sep = "=", collapse = ", ")
      else as.character(x)
    }, character(1))
  ),
  diann_spectral_library_cfg = parse_diann_cfg(
    snakemake@input[["diann_spectral_lib_config"]]
  ),
  diann_run_cfg = parse_diann_cfg(
    snakemake@input[["diann_run_config"]]
  ),
  runtime = tibble::tibble(
    parameter = "snakemake_version",
    value     = snakemake_version
  )
)

provenance <- conduitR::create_provenance(
  workflow_version = workflow_version,
  config           = config_list
)

conduitR::log_with_timestamp("Constructing Conduit object from snakemake workflow files")

  # Create the  conduit object
conduit <- new("conduit",
QFeatures = QFeatures,
metrics = metrics,
database = database,
annotations = annotations,
taxonomy = taxonomy,
provenance = provenance
)

conduitR::log_with_timestamp("Calculating protein coverage per taxonomy, adding to Conduit metric slot")

conduit <- conduitR::add_protein_coverage_taxa_metrics(conduit)

conduitR::log_with_timestamp(paste0("Writing ",conduit_fp,"to file"))
saveRDS(conduit,conduit_fp)
end_time <- Sys.time()
conduitR::log_with_timestamp("Completed build_conduit.R script. Time taken: %.2f minutes", 
    as.numeric(difftime(end_time, start_time, units = "mins")))

# closing clogfile connection
sink(type = "message")
sink()
close(zz)