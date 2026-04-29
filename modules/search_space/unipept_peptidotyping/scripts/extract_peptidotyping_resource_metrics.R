# =============================================================================
# Extract Peptidotyping Resources Metrics (Streaming, Memory-Safe)
# =============================================================================
# This script extracts metrics from the peptidotyping resources.
# It counts the number of peptides that belong to each taxa in the sequence index,
# and appends the full taxonomy creating a dataframe that is easy to work with.
# =============================================================================

# =============================================================================
# Setup and Logging
# =============================================================================
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz, append = TRUE)
sink(zz, type = "message")

start_time <- Sys.time()
conduitR::log_with_timestamp(
  "Starting extract_peptidotyping_resource_metrics.R script"
)

# =============================================================================
# Input / Output
# =============================================================================
sequences_file <- snakemake@input[["sequences"]]
taxons_file <- snakemake@input[["taxons"]]

peptidotyping_resource_metrics_fp <- snakemake@output[["peptidotyping_resource_metrics"]]

conduitR::log_with_timestamp("Input files: %s, %s", sequences_file, taxons_file)
conduitR::log_with_timestamp("Output file: %s", peptidotyping_resource_metrics_fp)

# =============================================================================
# Load Taxonomy (small file, safe to load in memory)
# =============================================================================
conduitR::log_with_timestamp("Reading taxon file")

taxons <- readr::read_tsv(
  pipe(paste("lz4 -d -c", shQuote(taxons_file))),
  col_names = c("id", "name", "rank", "parent_id"),
  col_types = "cccc"
)

conduitR::log_with_timestamp("Loaded %d taxons", nrow(taxons))

# =============================================================================
# Process Sequences File in Chunks (Memory-Safe)
# =============================================================================
# We only need the first four columns of the sequence file. 
# The first four columns of the sequence file are the id, peptide sequence, and LCA.
# 1. id: Internal identifier of this sequence. This identifier is used by the peptides.tsv.gz to refer to sequences in this file.
# 2. sequence: String-representation of the sequence of amino acids that this peptide consists of.
# 3. lca: Lowest common ancestor of the taxa associated with all proteins that contain this peptide sequence (in the case that I and L are not considered equal).
# 4. lca_il: Lowest common ancestor of the taxa associated with all proteins that contain this peptide sequence (in the case that I and L are considered equal).

conduitR::log_with_timestamp("Counting peptides per taxa using shell tools (most efficient for large files)")

# For massive files (292GB+), use shell tools to count directly
# This is much faster and more memory-efficient than reading into R first
# Pipeline: decompress -> extract column 4 (lca_il) -> sort -> count with uniq -c
conduitR::log_with_timestamp("Running: lz4 -d -c | cut -f4 | sort | uniq -c")

counts_con <- pipe(paste(
  "lz4 -d -c", shQuote(sequences_file),
  "| cut -f4",
  "| sort",
  "| uniq -c",
  "| awk '{print $2\"\\t\"$1}'"
))

peptide_counts <- read.table(
  counts_con,
  col.names = c("lca_il", "n_peptides"),
  stringsAsFactors = FALSE
)

conduitR::log_with_timestamp("Found %d unique taxa", nrow(peptide_counts))

# =============================================================================
# Join with Taxonomy and Create Final Metrics
# =============================================================================
conduitR::log_with_timestamp("Joining peptide counts with taxonomy")

peptidotyping_resource_metrics <- peptide_counts |>
  dplyr::left_join(taxons, by = c("lca_il" = "id")) |>
  dplyr::select(lca_il, name, rank, parent_id, n_peptides) |>
  dplyr::arrange(dplyr::desc(n_peptides))

conduitR::log_with_timestamp("Created metrics for %d taxa", nrow(peptidotyping_resource_metrics))

# =============================================================================
# Write Output
# =============================================================================
conduitR::log_with_timestamp("Writing metrics to file")

readr::write_tsv(peptidotyping_resource_metrics, peptidotyping_resource_metrics_fp)

# =============================================================================
# Cleanup and Logging
# =============================================================================
end_time <- Sys.time()
elapsed_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))

conduitR::log_with_timestamp(
  "Completed extract_peptidotyping_resource_metrics.R script. Time taken: %.2f minutes", 
  elapsed_minutes
)

# Close log file connections
sink(type = "message")
sink()
close(zz)