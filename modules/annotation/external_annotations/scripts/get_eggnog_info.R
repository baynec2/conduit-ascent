################################################################################
# Getting Eggnog Annotations for each Detected Protein ID
################################################################################
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

start_time <- Sys.time()

conduitR::log_with_timestamp("Running get_eggnog_info.R script")
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[["uniprot_annotated_protein_info"]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[["eggnog_resource"]]))
conduitR::log_with_timestamp(paste0("Output file: ",snakemake@output[["eggnog_info"]]))

#Defining files
# Inputs
uniprot_annotated_protein_info_fp = snakemake@input[["uniprot_annotated_protein_info"]]
eggnog_resource_fp = snakemake@input[["eggnog_resource"]]

# Outputs
eggnog_info_fp = snakemake@output[["eggnog_info"]]
eggnog_code_info_fp = snakemake@output[["eggnog_code_info"]]

# Reading in files
conduitR::log_with_timestamp("Reading in files")

uniprot_annotated_protein_info = readr::read_delim(uniprot_annotated_protein_info_fp)

# Empty input (UniParc / non-reference proteome -> no UniProtKB annotations):
# emit empty, correctly-typed tables. NB these are the UniProt-derived eggNOG
# xrefs; eggNOG-mapper's own annotations are produced on a separate path.
.empty_annotation <- tibble::tibble(
  protein_id = character(), annotation_type = character(),
  term = character(), description = character()
)
if (nrow(uniprot_annotated_protein_info) == 0L) {
  conduitR::log_with_timestamp("No UniProtKB annotations; writing empty eggNOG tables")
  readr::write_delim(.empty_annotation, eggnog_info_fp)
  readr::write_delim(.empty_annotation, eggnog_code_info_fp)
  sink(type = "message"); sink(); close(zz)
  quit(save = "no", status = 0)
}

detected_eggnog = uniprot_annotated_protein_info |>
dplyr::select(protein_id,xref_eggnog)|>
tidyr::separate_longer_delim(xref_eggnog,delim = ";")|>
dplyr::filter(xref_eggnog != "")|>
dplyr::distinct()

eggnog_resource = readr::read_tsv(eggnog_resource_fp,col_names = c("row","xref_eggnog","code","description"))

eggnog_combined = dplyr::inner_join(detected_eggnog,eggnog_resource,by = "xref_eggnog")

# OG lookup table
eggnog_code_lookup <- tibble::tibble(
  code = c("A","B","C","D","E","F","G","H","I","J",
           "K","L","M","N","O","P","Q","R","S","T",
           "U","V","W","Y","Z"),
  code_description = c(
    "RNA processing and modification",
    "Chromatin structure and dynamics",
    "Energy production and conversion",
    "Cell cycle control, cell division",
    "Amino acid transport and metabolism",
    "Nucleotide transport and metabolism",
    "Carbohydrate transport and metabolism",
    "Coenzyme transport and metabolism",
    "Lipid transport and metabolism",
    "Translation, ribosomal structure and biogenesis",
    "Transcription",
    "Replication, recombination and repair",
    "Cell wall/membrane/envelope biogenesis",
    "Cell motility",
    "Posttranslational modification, protein turnover, chaperones",
    "Inorganic ion transport and metabolism",
    "Secondary metabolites biosynthesis, transport and catabolism",
    "General function prediction only",
    "Function unknown",
    "Signal transduction mechanisms",
    "Intracellular trafficking, secretion, and vesicular transport",
    "Defense mechanisms",
    "Extracellular structures",
    "Nuclear structure",
    "Cytoskeleton"
  )
)

# Eggnog code
conduitR::log_with_timestamp("Extracting Eggnog codes")

eggnog_codes = eggnog_combined |>
  dplyr::mutate(code = stringr::str_split(code, "")) |>  # split into list of single letters
  tidyr::unnest(code) |> # make one row per letter
  dplyr::inner_join(eggnog_code_lookup,by = "code") |>
  dplyr::select(protein_id,term = code,description = code_description)|>
  dplyr::mutate(annotation_type = "uniprot_eggnog_code")|>
  dplyr::select(protein_id,annotation_type,term,description)|>
  dplyr::distinct()

# Eggnog description
conduitR::log_with_timestamp("Extracting Eggnog descriptions")

eggnog <- eggnog_combined |>
dplyr::select(protein_id,term = xref_eggnog,description)|>
dplyr::mutate(annotation_type = "uniprot_eggnog") |>
dplyr::select(protein_id,annotation_type,term,description)|>
dplyr::distinct()

conduitR::log_with_timestamp("Writing eggnog code annotaions to file")
readr::write_delim(eggnog_codes,eggnog_code_info_fp)

conduitR::log_with_timestamp("Writing eggnog annotations to file")
readr::write_delim(eggnog,eggnog_info_fp)

end_time <- Sys.time()
conduitR::log_with_timestamp("Completed get_eggnog_info.R script. Time taken: %.2f minutes", 
    as.numeric(difftime(end_time, start_time, units = "mins")))
# closing clogfile connection
sink(type = "message")
sink()
close(zz)