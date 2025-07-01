################################################################################
# Creating Conduit Object
################################################################################
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages
start_time <- Sys.time()

conduitR::log_with_timestamp("Running create_conduit_obj.R script")
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[["qf"]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[["combined_metrics"]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[["database_taxonomy"]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[["detected_protein_taxonomy"]]))

conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[["conduit_obj"]]))

# Defining inputs
qf = snakemake@input[["qf"]]
combined_metrics = snakemake@input[["combined_metrics"]]
database_taxonomy = snakemake@input[["database_taxonomy"]]
detected_protein_taxonomy = snakemake@input[["detected_protein_taxonomy"]]

# Defining outputs
conduit_obj = snakemake@output[["conduit_obj"]]

conduitR::log_with_timestamp("Creating Conduit object")
conduit = conduitR::create_conduit_obj(qf,
                                       combined_metrics,
                                       database_taxonomy,
                                       detected_protein_taxonomy)
conduitR::log_with_timestamp(paste0("Saving Conduit object to ", conduit_obj))
# Saving to File
saveRDS(conduit,conduit_obj)

end_time <- Sys.time()
conduitR::log_with_timestamp("create_conduit_obj.R script completed. Time taken: %.2f minutes", 
                              as.numeric(difftime(end_time, start_time, units = "mins")))
# closing clogfile connection
sink(type = "message")
sink()
close(zz)