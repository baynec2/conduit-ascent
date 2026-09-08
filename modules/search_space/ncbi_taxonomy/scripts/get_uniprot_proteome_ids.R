################################################################################
# Getting uniprot proteome IDs from NCBI organism IDs
################################################################################
## Opening Log File 
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

## Cap conduitR's parallel worker pool to this rule's Snakemake allocation.
# conduitR::get_proteome_ids_from_organism_ids() sizes its future/furrr pool
# from future::availableCores() - 1 (= parallelly::availableCores()), which
# otherwise reports every physical core on the node and oversubscribes when
# Snakemake scheduled this rule with fewer threads. Setting the `custom`
# availableCores() method makes it return this rule's thread count (parallelly
# takes the min across methods, so it never exceeds the node's real cores).
n_threads <- as.integer(snakemake@threads[[1]])
options(parallelly.availableCores.custom = function() n_threads)

start_time <- Sys.time()

# Now everything from print(), message(), warning() will go into the log file
conduitR::log_with_timestamp("Starting get_uniprot_proteome_ids.R script")
conduitR::log_with_timestamp("Input file: %s", snakemake@input[[1]])
conduitR::log_with_timestamp("Output file: %s", snakemake@output[[1]])

# Pure-function helper (parse_organism_ids) lives in a sibling file so it's
# reachable from testthat without snakemake@ globals.
snakemake@source("get_uniprot_proteome_ids_lib.R")

## Defining inputs and outputs from snakemake workflow
input_file <- snakemake@input[[1]]
proteome_ids_fp <- snakemake@output[[1]]

# Extracting additional ncbi taxa ids to append if specified
append_additional_ncbi_taxa_id <- snakemake@config$append_additional_ncbi_taxa_id

conduitR::log_with_timestamp("Making the database_resources_directory if it doesn't exist.")

conduitR::log_with_timestamp("Reading organism IDs from the input file.")
raw_df <- readr::read_tsv(input_file, show_col_types = FALSE)

# Empty detection (e.g. a first pass that identified nothing): no organism IDs
# to resolve. Emit an empty, schema-correct proteome table so the run resolves
# to an empty (no-detection) conduit downstream instead of querying UniProt with
# an empty set.
if (nrow(raw_df) == 0L && isFALSE(snakemake@config$append_additional_ncbi_taxa_id)) {
  conduitR::log_with_timestamp("No organism IDs — writing empty proteome_ids table.")
  empty_cols <- c("initial_proteome_id", "organism_id", "organism", "protein_count",
                  "proteome_type", "redundant_to", "genome_assembly_id",
                  "genome_assembly_level", "annotation_score", "parent_id",
                  "child_rank", "proteome_id")
  empty_df <- stats::setNames(
    lapply(empty_cols, function(x) character(0)), empty_cols
  ) |> tibble::as_tibble()
  readr::write_tsv(empty_df, proteome_ids_fp)
  sink(type = "message"); sink(); close(zz)
  quit(save = "no", status = 0)
}

# Compute the pre-append id set so we can log "already present" vs "appended"
# accurately. The parse_organism_ids helper handles both cases internally and
# returns the final deduplicated vector.
pre_append_ids <- parse_organism_ids(raw_df, append_id = FALSE)
organism_ids   <- parse_organism_ids(raw_df, append_additional_ncbi_taxa_id)

if (!isFALSE(append_additional_ncbi_taxa_id)) {
  if (as.integer(append_additional_ncbi_taxa_id) %in% pre_append_ids) {
    conduitR::log_with_timestamp(
      "User-specified NCBI organism ID %s is already present in the data.",
      append_additional_ncbi_taxa_id
    )
  } else {
    conduitR::log_with_timestamp(
      "Appending user-specified NCBI organism ID %s to the NCBI taxonomy IDs to search.",
      append_additional_ncbi_taxa_id
    )
  }
}

conduitR::log_with_timestamp("Finding the uniprot proteome ID cooresponding to each NCBI taxonomy ID")

# Using conduitR to get proteome ids corresponding to an organism id
proteome_ids_df <- conduitR::get_proteome_ids_from_organism_ids(organism_ids)

conduitR::log_with_timestamp("Searching for higher-quality UniProt proteome IDs from closely related taxa")

# UniProt API results may not always provide the optimal proteome for a given NCBI taxonomic ID.
# Returned entries can be redundant, absent from UniProtKB, or assigned to strains without curated proteomes.
# We therefore identify higher-quality proteomes from related taxa and select the best available match.
best_proteome_ids <- conduitR::get_better_proteome_ids(proteome_ids_df)

# Renaming to be compatible with uniprot_proteome_ids search space workflow. Essentially this will
# extract proteome_id from the dataframe column "proteome_id".
best_proteome_ids <- best_proteome_ids |>
dplyr::rename(initial_proteome_id = proteome_id,
proteome_id = selected_proteome_id)

# Save to file
conduitR::log_with_timestamp("Saving proteome ids to file.")

readr::write_tsv(best_proteome_ids,proteome_ids_fp)

end_time <- Sys.time()

conduitR::log_with_timestamp("Completed get_uniprot_proteome_ids.R script. Time taken: %.2f minutes", 
    as.numeric(difftime(end_time, start_time, units = "mins")))

# closing logfile connection
conduitR::log_with_timestamp("get_uniprot_proteome_ids.R script complete.")
# closing logfile connection
sink(type = "message")
sink()
close(zz)