################################################################################
# Get Supplementary Annotations From Uniprot
# Here we are taking the uniref matches provided by bakta, and using those 
# to annotate the data. This works conviently with the rest of the conduit 
# infastructure. 
################################################################################
## Opening Log File 
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

# Get input and output files from Snakemake workflow
bakta_annotated_protein_info_fp <- snakemake@input[["bakta_annotated_protein_info"]]

# Output
uniprot_annotated_protein_info_fp <- snakemake@output[["uniprot_annotated_protein_info"]]


start_time <- Sys.time()

# Logging Inputs and outputs
conduitR::log_with_timestamp("Running get_supplementary_annotations_from_uniprot.R script")
conduitR::log_with_timestamp(paste0("Input file: ",bakta_annotated_protein_info_fp))
conduitR::log_with_timestamp(paste0("Output file: ", bakta_annotated_protein_info_fp))

bakta_annotated_protein_info <- readr::read_delim(bakta_annotated_protein_info_fp)

uniprot_ids <- bakta_annotated_protein_info$xref_uniprot[!is.na(bakta_annotated_protein_info$xref_uniprot)]

# Dealing with annotations derived from bakta vs uniref protein ids 
# If a Bakta annotation is there, we will use that. If not, we will use the associated uniprot annotation. 
# Which annotation is being used will not explicitly be recorded, but I think that is an okay approach.
conduitR::log_with_timestamp("Getting Uniprot Annotations Via Uniprot API")
uniprot_annotations <- conduitR::get_annotations_from_uniprot(uniprot_ids)|>
dplyr::distinct()

conduitR::log_with_timestamp("Merging Uniprot Annotations with Bakta Annotations")
uniprot_annotated_protein_info <- bakta_annotated_protein_info |>
# Join the dataframes - columns that exist in both will get suffix to designate the one from bakta
# Uniprot ones will have no suffix consistant with the rest of the workflow.
  dplyr::left_join(uniprot_annotations, by = c("xref_uniprot" = "accession"), suffix = c("_bakta", ""))


# Writing to file 
conduitR::log_with_timestamp(paste0("Writing annotations to ", uniprot_annotated_protein_info_fp))
readr::write_delim(uniprot_annotated_protein_info,
                   uniprot_annotated_protein_info_fp)

end_time <- Sys.time()

conduitR::log_with_timestamp("Completed get_supplementary_annotations_from_uniprot.R script. Time taken: %.2f minutes", 
    as.numeric(difftime(end_time, start_time, units = "mins")))

# closing logfile connection
sink(type = "message")
sink()
close(zz)