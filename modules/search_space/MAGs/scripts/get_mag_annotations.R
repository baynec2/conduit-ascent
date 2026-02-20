################################################################################
# Get MAG Annotations From Bakta Output
################################################################################
## Opening Log File 
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

# Get input and output files from Snakemake workflow
bakta_dirs <- snakemake@input[["bakta_dirs"]]
# Output
mag_annotations_fp <- snakemake@output[["mag_annotations"]]

start_time <- Sys.time()
mag_names <- basename(bakta_dirs)
cds_files <- paste0(bakta_dirs,"/",mag_names,".tsv")

conduitR::log_with_timestamp("Running get_mag_annotations.R script")
conduitR::log_with_timestamp(paste0("Input files: ", cds_files))
conduitR::log_with_timestamp(paste0("Output file: ", mag_annotations_fp))

conduitR::log_with_timestamp("Consolidating MAG annotations into one file")

# Function to read in annotation files and give them magname
read_annotation_files = function(filepath){
mag_name <- gsub(".tsv","",basename(filepath))
    data <- readr::read_tsv(filepath,skip=5) |>
    dplyr::mutate(mag = mag_name)
    return(data)
}

combined_annotations <- purrr::map_df(cds_files,read_annotation_files)|>
dplyr::select(mag,dplyr::everything())

# Formatting 
conduitR::log_with_timestamp("Formatting MAG annotations")

long_annotations = combined_annotations |>
    tidyr::separate_rows(DbXrefs, sep = ",\\s*") |> # each DbXref in its own row
    dplyr::mutate(xref_name = dplyr::case_when(grepl("^SO.*",DbXrefs) ~ "xref_SO",
                                           grepl("^UniRef:UniRef50_.*",DbXrefs) ~ "xref_uniref50",
                                           grepl("^UniRef:UniRef90_.*",DbXrefs) ~ "xref_uniref90",
                                           grepl("^UniRef:UniRef100_.*",DbXrefs) ~ "xref_uniref100",
                                           grepl("^PFAM.*",DbXrefs) ~ "xref_pfam",
                                           grepl("^GO.*",DbXrefs) ~ "xref_go",
                                           grepl("^COG.*",DbXrefs) ~ "xref_cog",
                                           grepl("^EC.*",DbXrefs) ~ "xref_ec",
                                           TRUE ~ "unsuported_annotation"
                                           ),
                    DbXrefs = gsub("^UniRef:UniRef[0-9]*_|PFAM:","",DbXrefs))
# Step 2: Separate db prefix from value
mag_annotations <- long_annotations  |> 
tidyr::pivot_wider(names_from = xref_name,values_from = DbXrefs, values_fn = list(DbXrefs = ~paste(., collapse = ";")))

conduitR::log_with_timestamp("Writing MAG annotations to file.")

readr::write_delim(mag_annotations,mag_annotations_fp)

end_time <- Sys.time()

conduitR::log_with_timestamp("Completed get_mag_annotations.R script. Time taken: %.2f minutes", 
    as.numeric(difftime(end_time, start_time, units = "mins")))

# closing logfile connection
sink(type = "message")
sink()
close(zz)
