# =============================================================================
# Setup and Logging
# =============================================================================
# Open log file for both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

# Record start time
start_time <- Sys.time()
conduitR::log_with_timestamp("Starting infer_species_presence.R script")

# =============================================================================
#  Setting up Input and Output Files
# =============================================================================
# Input files
first_pass_diann_parquet <- snakemake@input[["first_pass_diann"]]
taxon_specific_peptide_db <- snakemake@input[["taxon_specific_peptide_db"]]

# Test 
#first_pass_diann_parquet <-  "experiments/defined_community/input/database_resources/proteotyping/first_pass_diann.parquet"
#taxon_specific_peptide_db <- "resources/databases/lca_filtered_taxa.tsv"

conduitR::log_with_timestamp("Input files: %s, %s", first_pass_diann_parquet, taxon_specific_peptide_db)

# Output files
first_pass_search_taxa_metrics_fp <- snakemake@output[["first_pass_search_taxa_metrics"]]
called_taxa_per_run_fp <- snakemake@output[["called_taxa_per_run"]]
ncbi_taxonomy_id_fp <- snakemake@output[["ncbi_taxonomy_id"]]

conduitR::log_with_timestamp("Output files: %s, %s, %s", 
first_pass_search_taxa_metrics_fp,
called_taxa_per_run_fp,
 ncbi_taxonomy_id_fp)

# =============================================================================
#  Generating First Pass Search Taxa Metrics
# =============================================================================

# Reading in precursor file with parquet
conduitR::log_with_timestamp("Reading in: %s", first_pass_diann_parquet)
precursors <- arrow::read_parquet(first_pass_diann_parquet)

# Summarising precursors to peptide
conduitR::log_with_timestamp("Summing precursors to peptides")

peptides <- precursors |> 
  dplyr::group_by(Run,Protein.Names,Stripped.Sequence,Protein.Group) |> 
  dplyr::summarise(n_precursors = dplyr::n(),
            PG.MaxLFQ = mean(PG.MaxLFQ, na.rm=TRUE))

conduitR::log_with_timestamp("Counting the number of peptides detected per species or strain")

# Loading first pass database
conduitR::log_with_timestamp("Loading first pass database composition")

first_pass_db <- readr::read_tsv(taxon_specific_peptide_db)|>  
dplyr::select(id,sequence,parent_id,organism_id = lca_il,rank, name)

conduitR::log_with_timestamp("Determining peptide metrics")

detected_peptide_count <- peptides |> 
  dplyr::left_join(first_pass_db,by = c("Stripped.Sequence" = "sequence"))|> 
  dplyr::group_by(Run,parent_id, organism_id, rank, taxa = Protein.Names) |> 
  dplyr::summarise(n_peptides_detected = dplyr::n())


database_peptide_count <- first_pass_db |>  
  dplyr::group_by(organism_id) |>  
  dplyr::summarise(n_peptides_in_db = dplyr::n())

# total number of proteotypic peptides across all organisms.
n_total_peptides <- nrow(first_pass_db)

# Total peptides observed (after FDR filtering)
n_total_detected_peptides <- nrow(peptides)

conduitR::log_with_timestamp("Combining detected peptide information with database info")
conduitR::log_with_timestamp("Calculating P values")

# Combine everything
first_pass_metrics <- dplyr::left_join(detected_peptide_count,
                                      database_peptide_count,
                                      by ="organism_id") |>  
  tidyr::separate(taxa,sep = "_",into = c("rank","name")) |>
  dplyr::mutate(genus = stringr::word(name, 1, sep = "-"),
                species = stringr::word(name, 2, sep = "-"),
                species = paste0(genus," ",species),
                .after = name) |>  
  dplyr::mutate(n_total_peptides = n_total_peptides,
         n_total_detected_peptides = n_total_detected_peptides,
         # This general approach is very simplistic and flawed. Doesn't work 
         # Super well. Need to try some sort of Bayesian Approach.
         # Priors need to be set by some sort of reasonable metric.
         p0 = n_peptides_in_db/n_total_peptides,
         expected = n_total_detected_peptides * 0.01 * p0,
         coverage = round(n_peptides_detected/n_peptides_in_db * 100,1),
         p_value = pbinom(n_peptides_detected - 1,
                     size = n_total_detected_peptides,
                     prob = 0.01 * p0,
                     lower.tail = FALSE),
            # Add multiple testing correction
        p_value_fdr = p.adjust(p_value, method = "fdr")
         )

conduitR::log_with_timestamp("Writing first pass metrics to file")

readr::write_tsv(first_pass_metrics, first_pass_search_taxa_metrics_fp)

# =============================================================================
#  Reducing Strain level to species level
# =============================================================================
conduitR::log_with_timestamp("Aggregating strain level taxonomy at species level")

species_level_first_pass_metrics <- first_pass_metrics |>
  dplyr::ungroup()|>
  dplyr::mutate(
    organism_id = dplyr::case_when(
      rank == "strain" ~ parent_id,
      TRUE ~ organism_id
    )
  ) |> dplyr::select(-parent_id,-p_value,-p_value_fdr)|>
  dplyr::group_by(Run,organism_id,species)|>
  dplyr::summarise(n_peptides_detected = sum(n_peptides_detected),
  n_peptides_in_db = sum(n_peptides_in_db),
  coverage = n_peptides_detected/n_peptides_in_db * 100)

# =============================================================================
#  Applying Thresholds to Call Species
# =============================================================================
# Applying thresholds to call what species are present.
# Note these are currently somewhat arbitrary and need to be improved
# Worked well with E.coli test.
conduitR::log_with_timestamp("Applying thresholds to call species")

called_taxa_per_run <- species_level_first_pass_metrics |>
# Thresholds are applied here
  dplyr::filter(n_peptides_detected > 1,
         coverage > 2,
         n_peptides_in_db > 20) 


conduitR::log_with_timestamp("Writing file with taxa called per run (MS File)")

# Write to file, keeping all of the run info. 
readr::write_tsv(called_taxa_per_run,called_taxa_per_run_fp)

conduitR::log_with_timestamp("Writing identified taxa to ncbi_taxa_ids.tsv file")

ncbi_taxa_id <- called_taxa_per_run |>  
  dplyr::ungroup()|> 
  dplyr::select("organism_id") |>  
  dplyr::mutate(organism_type = dplyr::case_when(organism_id %in% c(9606,10090) ~ "host",
                                                 TRUE ~ "microbiome"))

readr::write_tsv(ncbi_taxa_id,ncbi_taxonomy_id_fp)

# =============================================================================
# Cleanup and Logging
# =============================================================================
end_time <- Sys.time()
elapsed_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))

conduitR::log_with_timestamp(
  "Completed infer_species_presence.R script. Time taken: %.2f minutes", 
  elapsed_minutes
)

# Close log file connections
sink(type = "message")
sink()
close(zz)