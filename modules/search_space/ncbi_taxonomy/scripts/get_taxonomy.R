################################################################################
# Getting Taxonomy From Organism IDs
################################################################################
## Opening Log File 
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

# Get input and output files from Snakemake workflow
proteome_ids_fp <- snakemake@input[["proteome_ids"]]
output_file <- snakemake@output[["taxonomy"]]

start_time <- Sys.time()
conduitR::log_with_timestamp("Running get_taxonomy.R script")
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[1]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[[1]]))

conduitR::log_with_timestamp("Reading selected proteome ids from the input file.")

# Read organism IDs from the input file
proteome_ids <- readr::read_delim(proteome_ids_fp,
                                  col_types = "cc")|>
                                  dplyr::pull(selected_proteome_id) |>
  unique()

conduitR::log_with_timestamp("Getting NCBI Taxonomy Ids corresponding to selected proteome from uniprot.")

organism_ids = conduitR::get_taxonomy_from_proteome_ids(proteome_ids)|>
  dplyr::pull(organsim_id)|>
  unique()

conduitR::log_with_timestamp("Getting Full NCBI Taxonomy corresponding to NCBI ID from NCBI API. ")
# Pull all taxonomy information
taxonomy = conduitR::get_ncbi_taxonomy(organism_ids)

conduitR::log_with_timestamp("Finished downloading Taxonomy Information from NCBI API.")

taxonomy = taxonomy |>
  dplyr::left_join(proteome_ids,by = c("organism_id"= "organism_id"))|>
  # This is probably not the best approach, but I can't think of a better way to do it for now
  dplyr::mutate(organism_type = dplyr::case_when(organism_id %in% c(9606,10090) ~ "host",
  TRUE ~ "microbiome"))

conduitR::log_with_timestamp("Writing Taxonomy Information to file.")
# Writing to file.
readr::write_delim(taxonomy,output_file)

end_time <- Sys.time()

conduitR::log_with_timestamp("Completed get_taxonomy.R script. Time taken: %.2f minutes", 
    as.numeric(difftime(end_time, start_time, units = "mins")))

# closing clogfile connection
sink(type = "message")
sink()
close(zz)