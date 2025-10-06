################################################################################
# Building a Conduit Object 
################################################################################
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

start_time <- Sys.time()

conduitR::log_with_timestamp("Running build_conduit.R script")

conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[["diann_stats"]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[["qfeatures"]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[["database"]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[["annotations"]]))

conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[["conduit"]]))

#Defining files
# Inputs
diann_stats_fp = snakemake@input[["diann_stats"]]
QFeatures_fp = snakemake@input[["qfeatures"]]
database_fp = snakemake@input[["database"]]
annotations_fp = snakemake@input[["annotations"]]
# Outputs
conduit_fp = snakemake@output[["conduit"]]

conduitR::log_with_timestamp("Constructing Conduit object from snakemake workflow files")

# Constructing Conduit object with the files that were produced
conduit <- conduitR::create_conduit_obj(QFeatures_fp,
                                       diann_stats_fp,
                                       database_fp,
                                       annotations_fp)

conduitR::log_with_timestamp("Calculating protein coverage per taxonomy, adding to Conduit metric slot")

conduit <- conduitR::add_protein_coverage_taxa_metrics(conduit)

conduitR::log_with_timestamp(paste0("Writing ",conduit_fp,"to file"))
saveRDS(conduit,conduit_fp)
end_time <- Sys.time()
conduitR::log_with_timestamp("Completed build_conduit.R script. Time taken: %.2f minutes", 
    as.numeric(difftime(end_time, start_time, units = "mins")))

# closing clogfile connection
sink(type = "message")
sink()
close(zz)