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

conduitR::log_with_timestamp("Reading proteome ids from the input file.")

# Read organism IDs from the input file
proteome_id_df <- readr::read_delim(proteome_ids_fp,
                                  col_types = "cc")
                                  
                                  
proteome_ids <- proteome_id_df |>
  dplyr::pull(proteome_id) |>
  unique()

conduitR::log_with_timestamp("Getting NCBI Taxonomy Ids corresponding to selected proteome from uniprot.")

organism_ids = conduitR::get_taxonomy_from_proteome_ids(proteome_ids)|>
  dplyr::pull(organism_id)|>
  unique()

conduitR::log_with_timestamp("Getting Full NCBI Taxonomy corresponding to NCBI ID from NCBI API. ")

# Pull all taxonomy information from NCBI API
taxonomy = conduitR::get_ncbi_taxonomy(organism_ids)

conduitR::log_with_timestamp("Finished downloading Taxonomy Information from NCBI API.")

taxonomy = taxonomy |>
  dplyr::left_join(proteome_id_df, by = c("organism_id" = "organism_id"))

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