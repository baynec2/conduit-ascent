################################################################################
# Append Additional Organisms or Proteomes to MAG FASTA & Taxonomy
################################################################################

## Snakemake logging setup -----------------------------------------------------
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz, append = TRUE)       # stdout
sink(zz, type = "message")  # stderr/messages

start_time <- Sys.time()

## Helpers --------------------------------------------------------------------
is_missing <- function(x) {
  is.null(x) || isFALSE(x) || length(x) == 0 || identical(x, "")
}

append_proteomes <- function(fasta_set, taxonomy, proteome_ids, fasta_dir) {
  # Download UniProt FASTA(s)
  conduitR::get_fasta_file(proteome_ids, fasta_dir = fasta_dir)

  fasta_files <- list.files(fasta_dir, full.names = TRUE)
  if (length(fasta_files) == 0) {
    stop("No UniProt FASTA files found in fasta_dir")
  }

  uniprot_fasta <- Biostrings::readAAStringSet(fasta_files)
  additional_taxonomy <- conduitR::get_taxonomy_from_proteome_ids(proteome_ids)|> dplyr::pull(organism_id)
  
  additional_taxonomy_df <- conduitR::get_ncbi_taxonomy(additional_taxonomy)|>
   dplyr::mutate(organism_id = as.numeric(organism_id))

  list(
    fasta = c(fasta_set, uniprot_fasta),
    taxonomy = dplyr::bind_rows(taxonomy, additional_taxonomy_df)
  )
}

## Inputs ---------------------------------------------------------------------
mag_fasta_fp <- snakemake@input[["mag_fasta"]]
mag_taxonomy_fp <- snakemake@input[["mag_taxonomy"]]

## Outputs --------------------------------------------------------------------
uniprot_fasta_dir <- snakemake@output[["uniprot_fasta_dir"]]
fasta_fp <- snakemake@output[["fasta"]]
taxonomy_fp <- snakemake@output[["taxonomy"]]

## Config ---------------------------------------------------------------------
additional_proteome_id <- snakemake@config[["append_additional_proteome_id"]]
additional_ncbi_taxa_id <- snakemake@config[["append_additional_ncbi_taxa_id"]]

conduitR::log_with_timestamp("Running append_additional_organisms_or_proteomes.R")
conduitR::log_with_timestamp(
  paste0("Input files: ", mag_fasta_fp, " ", mag_taxonomy_fp)
)
conduitR::log_with_timestamp(
  paste0("Output files: ", fasta_fp, " ", taxonomy_fp)
)

## Read inputs ----------------------------------------------------------------
# FASTA via Biostrings
fasta <- Biostrings::readAAStringSet(mag_fasta_fp)

# Taxonomy table
taxonomy <- readr::read_delim(mag_taxonomy_fp, col_types = readr::cols())
conduitR::log_with_timestamp("Creating directory for uniprot fasta file")
dir.create(uniprot_fasta_dir)
## Control flow ---------------------------------------------------------------
if (is_missing(additional_proteome_id) && is_missing(additional_ncbi_taxa_id)) {
  conduitR::log_with_timestamp(
    "No additional proteomes or taxa specified; passing MAGs through unchanged"
  )

} else if (!is_missing(additional_proteome_id) && is_missing(additional_ncbi_taxa_id)) {
  conduitR::log_with_timestamp(
    paste0("Appending UniProt proteome ID(s): ", additional_proteome_id)
  )
  res <- append_proteomes(
    fasta_set = fasta,
    taxonomy = taxonomy,
    proteome_ids = additional_proteome_id,
    fasta_dir = uniprot_fasta_dir
  )

  fasta <- res$fasta
  taxonomy <- res$taxonomy

} else if (!is_missing(additional_ncbi_taxa_id) && is_missing(additional_proteome_id)) {
  conduitR::log_with_timestamp(
    paste0("Resolving NCBI taxa ID(s): ", additional_ncbi_taxa_id)
  )

  proteome_ids <- conduitR::get_proteome_ids_from_organism_ids(
    additional_ncbi_taxa_id
  ) |> 
  dplyr::pull(proteome_id)

  res <- append_proteomes(
    fasta_set = fasta,
    taxonomy = taxonomy,
    proteome_ids = proteome_ids,
    fasta_dir = uniprot_fasta_dir
  )

  fasta <- res$fasta
  taxonomy <- res$taxonomy

} else {
  stop(
    "Only one of append_additional_proteome_id or append_additional_ncbi_taxa_id may be provided"
  )
}

## Write outputs --------------------------------------------------------------
conduitR::log_with_timestamp("Writing updated FASTA and taxonomy files, unless unchaged")

Biostrings::writeXStringSet(fasta, fasta_fp)
readr::write_delim(taxonomy, taxonomy_fp)

end_time <- Sys.time()

conduitR::log_with_timestamp(
  sprintf(
    "Completed append_additional_organisms_or_proteomes.R (%.2f minutes)",
    as.numeric(difftime(end_time, start_time, units = "mins"))
  )
)

## Cleanup --------------------------------------------------------------------
sink(type = "message")
sink()
close(zz)
