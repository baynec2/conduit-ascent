################################################################################
# Adding Protein Group Annotations to QFeatures object
################################################################################
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

start_time <- Sys.time()

conduitR::log_with_timestamp("Running add_annotations_to_qfeatures.R script")

conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[["qf"]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[["uniprot_annotated_protein_info"]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[["annotated_qf"]]))

#Defining files
# Inputs
qf_fp = snakemake@input[["qf"]]
uniprot_annotated_protein_info_fp = snakemake@input[["uniprot_annotated_protein_info"]]

# Outputs
annotated_qf_fp = snakemake@output[["annotated_qf"]]

# Reading in files
conduitR::log_with_timestamp("Reading in files")
uniprot_annotated_protein_info = readr::read_delim(uniprot_annotated_protein_info_fp)
qf = readRDS(qf_fp)
# Adding taxonomy annotations
conduitR::log_with_timestamp("Adding taxonomy information to QFeatures, handling assay links, and summarizing")
qf = conduitR::add_taxonomy_to_qf(qf,uniprot_annotated_protein_info)
# Adding go annotations
conduitR::log_with_timestamp("Adding GO information to QFeatures, handling assay links, and summarizing")
qf = conduitR::add_annotation_to_qf(qf,uniprot_annotated_protein_info)
# Adding kegg annotations
conduitR::log_with_timestamp("Adding KEGG information to QFeatures, handling assay links, and summarizing")
qf = conduitR::add_annotation_to_qf(qf,uniprot_annotated_protein_info,xref_kegg,"[^;]+(?=;)")

conduitR::log_with_timestamp(paste0("Writing Qfeatures object to ", annotated_qf_fp))

saveRDS(qf,annotated_qf_fp)
end_time <- Sys.time()
conduitR::log_with_timestamp("Completed add_annotations_to_qfeatures script. Time taken: %.2f minutes", 
    as.numeric(difftime(end_time, start_time, units = "mins")))

# closing clogfile connection
sink(type = "message")
sink()
close(zz)