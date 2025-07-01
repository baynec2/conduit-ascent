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
input_file <- snakemake@input[[1]]
proteome_ids_file <- snakemake@input[[2]]
output_file <- snakemake@output[[1]]

start_time <- Sys.time()
conduitR::log_with_timestamp("Running get_taxonomy.R script")
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[1]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[[1]]))

conduitR::log_with_timestamp("Reading organism IDs from the input file.")

# Read organism IDs from the input file
organism_txt <- readr::read_delim(input_file,
                                  col_types = "cc")

organism_ids = organism_txt |>
  dplyr::pull(organism_id)

conduitR::log_with_timestamp("Downloaing Taxonomy Information from NCBI API.")

# Pull all taxonomy information
taxonomy = conduitR::get_ncbi_taxonomy(organism_ids) |>
  dplyr::mutate(organism_id = as.character(organism_id))|>
  dplyr::inner_join(organism_txt,by = "organism_id") |>
  dplyr::select("organism_type","organism_id","domain","kingdom","phylum",
                "class","order","family","genus","species"
    )

conduitR::log_with_timestamp("Finished downloading Taxonomy Information from NCBI API.")

# Modify taxonomy file to include proteome ids
proteome_ids = readr::read_delim(proteome_ids_file) |>
dplyr::select("Proteome Id","organism_id" = "Organism Id","reference","downloaded_by_conduit","download_info")|>
dplyr::mutate(organism_id = as.character(organism_id))


taxonomy = taxonomy |>
  dplyr::left_join(proteome_ids,by = c("organism_id"= "organism_id"))

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