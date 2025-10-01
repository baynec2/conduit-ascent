################################################################################
# Processing Diann Parquet file to make a QFeatures object
################################################################################
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

start_time <- Sys.time()

conduitR::log_with_timestamp("Running generate_qfeatures_from_diann_parquet.R script")
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[["diann_parquet"]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[["qf"]]))

# Defining files
# Inputs
diann_parquet_fp = snakemake@input[["diann_parquet"]]
# Outputs
qf_fp = snakemake@output[["qf"]]
conduitR::log_with_timestamp(paste0("Reading in diann parquet file from ", diann_parquet_fp))
conduitR::log_with_timestamp("Processing diann parquet file to qfeatures object")

qf <- conduitR::diann_to_qfeatures(diann_parquet_fp)
conduitR::log_with_timestamp(paste0("Writing Qfeatures object to ", qf_fp))
saveRDS(qf,qf_fp)
end_time <- Sys.time()
conduitR::log_with_timestamp("Completed generate_qfeatures_from_diann_parquet script. Time taken: %.2f minutes", 
    as.numeric(difftime(end_time, start_time, units = "mins")))

# closing clogfile connection
sink(type = "message")
sink()
close(zz)
