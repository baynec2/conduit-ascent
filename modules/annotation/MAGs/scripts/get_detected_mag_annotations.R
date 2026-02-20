################################################################################
# Get Detected MAG Annotations
################################################################################
## Opening Log File 
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

# Get input and output files from Snakemake workflow
detected_protein_info_fp <- snakemake@input[["detected_protein_info"]]
mag_annotations_fp <- snakemake@input[["mag_annotations"]]

# Output
bakta_annotated_protein_info_fp <- snakemake@output[["bakta_annotated_protein_info"]]

# Testing
#detected_protein_info_fp = "/home/nanopore-catalyst/conduit/experiments/example/input/database_resources/protein_info.txt"
#mag_annotations_fp = "/home/nanopore-catalyst/conduit/experiments/example/input/database_resources/bakta/mag_annotations.txt"
#bakta_annotated_protein_info_fp = "/home/nanopore-catalyst/conduit/experiments/example/input/database_resources/bakta_annotated_protein_info.txt"

start_time <- Sys.time()
# Logging Inputs and outputs
conduitR::log_with_timestamp("Running get_detected_mag_annotations.R script")
conduitR::log_with_timestamp(paste0("Input files: ",detected_protein_info_fp," and ",mag_annotations_fp))
conduitR::log_with_timestamp(paste0("Output file: ", bakta_annotated_protein_info_fp))


# Reading in files
detected_protein_info <- readr::read_delim(detected_protein_info_fp)

mag_annotations <- readr::read_delim(mag_annotations_fp) |> 
dplyr::select(`Locus Tag`,dplyr::starts_with(("xref_")),
protein_name = Product)

# Merging files
bakta_annotated_protein_info <- dplyr::left_join(detected_protein_info,mag_annotations, by = c("protein_id" = "Locus Tag")) 

# Get the best uniprot annotation from available columns (priority: uniref100 > uniref90 > uniref50)
# Check which columns exist, then use coalesce with only available columns
available_uniref_cols <- names(bakta_annotated_protein_info)[grepl("xref_uniref.*",names(bakta_annotated_protein_info))]
sorted_uniref_cols <- available_uniref_cols[order(as.numeric(sub("xref_uniref", "", available_uniref_cols)),decreasing = TRUE)]

bakta_annotated_protein_info <- bakta_annotated_protein_info |>
  dplyr::mutate(
    xref_uniprot = dplyr::coalesce(!!!dplyr::syms(sorted_uniref_cols))
  )

mag_ids <- mag_annotations |> dplyr::pull("Locus Tag")

# If a uniprot proteome was appended, it won't have the uniprot_xref column. Adding that here. 
bakta_annotated_protein_info <- bakta_annotated_protein_info |>
  dplyr::mutate(
    xref_uniprot = dplyr::if_else(
      !is.na(protein_id) & !protein_id %in% mag_ids,
      protein_id,
      xref_uniprot
    )
  )
# Writing to file 
readr::write_delim(bakta_annotated_protein_info,bakta_annotated_protein_info_fp)

end_time <- Sys.time()

conduitR::log_with_timestamp("Completed get_detected_mag_annotations.R script. Time taken: %.2f minutes", 
    as.numeric(difftime(end_time, start_time, units = "mins")))

# closing logfile connection
sink(type = "message")
sink()
close(zz)
