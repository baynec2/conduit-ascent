################################################################################
# Getting Pfam Annotations for each Detected Protein ID
################################################################################

# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

start_time <- Sys.time()

conduitR::log_with_timestamp("Running get_pfam_info.R script")
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[["uniprot_annotated_protein_info"]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[["pfam_db"]]))

conduitR::log_with_timestamp(paste0("Output file: ",snakemake@output[["pfam_info"]]))

# Defining input files
uniprot_annotated_protein_info_fp = snakemake@input[["uniprot_annotated_protein_info"]]
pfam_db_fp = snakemake@input[["pfam_db"]]

# Defining output files
pfam_info_fp = snakemake@output[["pfam_info"]]

conduitR::log_with_timestamp("Loading pfam database resource")

# See https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/userman.txt
pfam_db = readr::read_delim(pfam_db_fp,col_select = c(1,4),col_names = FALSE)
colnames(pfam_db) <- c("term","description")

conduitR::log_with_timestamp("Loading detected pfams")

detected_pfams = readr::read_delim(uniprot_annotated_protein_info_fp)|>
dplyr::select(protein_id,term = xref_pfam)|>
tidyr::separate_longer_delim(cols = "term",delim = ";")|>
dplyr::filter(term != "")|>
dplyr::filter(!is.na(term))

conduitR::log_with_timestamp("Mapping detected pfams to their identities")

pfam_info = dplyr::inner_join(detected_pfams,pfam_db, by = "term") |>
dplyr::mutate(annotation_type = "pfam") |>
dplyr::select(protein_id,annotation_type,term,description)|>
dplyr::distinct()

conduitR::log_with_timestamp("Saving annotations to pfam_info.txt")
readr::write_delim(pfam_info,pfam_info_fp)

end_time <- Sys.time()
conduitR::log_with_timestamp("Completed get_pfam_info.R script. Time taken: %.2f minutes", 
    as.numeric(difftime(end_time, start_time, units = "mins")))
# closing clogfile connection
sink(type = "message")
sink()
close(zz)
