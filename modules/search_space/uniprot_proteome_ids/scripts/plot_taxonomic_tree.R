################################################################################
# Plot Taxonomic Tree
################################################################################
## Opening Log File 
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

conduitR::log_with_timestamp("Running plot_taxonomic_tree.R script")
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[1]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[[1]]))

start_time <- Sys.time()

# Get input and output files from Snakemake workflow
taxonomy_data_fp <- snakemake@input[[1]]
output_file_fp <- snakemake@output[[1]]

conduitR::log_with_timestamp(paste0("Loading taxonomy data from ", taxonomy_data_fp))
taxonomy_data_fp = readr::read_delim(taxonomy_data_fp,
                                     col_types = "cc") |>
  conduitR::plot_taxa_tree(node_color = "download_info")

conduitR::log_with_timestamp(paste0("Saving the plot to ", output_file_fp))
ggplot2::ggsave(output_file_fp,
                taxonomy_data_fp,
                height = 12,
                width = 12,
                units = "in")

end_time <- Sys.time()
conduitR::log_with_timestamp("Completed plot_taxonomic_tree.R script. Time taken: %.2f minutes", 
    as.numeric(difftime(end_time, start_time, units = "mins")))


# closing clogfile connection
sink(type = "message")
sink()
close(zz)
