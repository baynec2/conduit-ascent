################################################################################
# Get Genome Annotations From Bakta Output
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
genome_annotations_fp <- snakemake@output[["genome_annotations"]]

start_time <- Sys.time()
genome_names <- basename(bakta_dirs)
cds_files <- paste0(bakta_dirs,"/",genome_names,".tsv")

conduitR::log_with_timestamp("Running get_genome_annotations.R script")
conduitR::log_with_timestamp(paste0("Input files: ", cds_files))
conduitR::log_with_timestamp(paste0("Output file: ", genome_annotations_fp))

conduitR::log_with_timestamp("Consolidating genome annotations into one file")

# Function to read in annotation files and give them genome name
read_annotation_files = function(filepath){
genome_name <- gsub(".tsv","",basename(filepath))
    data <- readr::read_tsv(filepath,skip=5) |>
    dplyr::mutate(genome = genome_name)
    return(data)
}

combined_annotations <- purrr::map_df(cds_files,read_annotation_files)|>
dplyr::select(genome,dplyr::everything())

# Formatting
conduitR::log_with_timestamp("Formatting genome annotations")

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
                                           TRUE ~ "unsupported_annotation"
                                           ),
                    DbXrefs = gsub("^UniRef:UniRef[0-9]*_|PFAM:","",DbXrefs))

# Warn if any DbXrefs could not be categorised
unsupported <- long_annotations |>
  dplyr::filter(xref_name == "unsupported_annotation") |>
  dplyr::pull(DbXrefs) |>
  unique()
if (length(unsupported) > 0) {
  warning(paste0(length(unsupported), " unsupported DbXref annotation(s) were found and will be dropped: ",
                 paste(head(unsupported, 10), collapse = ", ")))
}

# Step 2: Separate db prefix from value
genome_annotations <- long_annotations |>
  dplyr::filter(xref_name != "unsupported_annotation") |>
tidyr::pivot_wider(names_from = xref_name,values_from = DbXrefs, values_fn = list(DbXrefs = ~paste(., collapse = ";")))

conduitR::log_with_timestamp("Writing genome annotations to file.")

readr::write_delim(genome_annotations,genome_annotations_fp)

end_time <- Sys.time()

conduitR::log_with_timestamp("Completed get_genome_annotations.R script. Time taken: %.2f minutes",
    as.numeric(difftime(end_time, start_time, units = "mins")))

# closing logfile connection
sink(type = "message")
sink()
close(zz)
