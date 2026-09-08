################################################################################
# Consolidating annotations
################################################################################
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

start_time <- Sys.time()

# Defining inputs
# Annotations
go_info_fp = snakemake@input[["go_info"]]
kegg_pathway_info_fp = snakemake@input[["kegg_pathway_info"]]
kegg_map_pathway_info_fp = snakemake@input[["kegg_map_pathway_info"]]
kegg_orthology_info_fp = snakemake@input[["kegg_orthology_info"]]
pfam_info_fp = snakemake@input[["pfam_info"]]
cazy_class_info_fp = snakemake@input[["cazy_class_info"]]
cazy_family_info_fp = snakemake@input[["cazy_family_info"]]
eggnog_info_fp = snakemake@input[["eggnog_info"]]
eggnog_code_info_fp = snakemake@input[["eggnog_code_info"]]
emapper_annotations_fp = snakemake@input[["emapper_annotations"]]

# Authoritative term-name dictionaries (for description backfill)
go_obo_fp       = snakemake@input[["go_obo"]]
kegg_ko_fp      = snakemake@input[["kegg_ko"]]
kegg_pathway_fp = snakemake@input[["kegg_pathway"]]
kegg_module_fp  = snakemake@input[["kegg_module"]]
kegg_brite_fp   = snakemake@input[["kegg_brite"]]
enzyme_dat_fp   = snakemake@input[["enzyme_dat"]]
pfam_clans_fp   = snakemake@input[["pfam_clans"]]

# qf
qf_fp = snakemake@input[["qf"]]


# Defining Output 
conduit_annotations_fp = snakemake@output[["conduit_annotations"]]

conduitR::log_with_timestamp("Combining all annotations together")

# Establishing vector of filenames
input_annotations <- c(go_info_fp,kegg_pathway_info_fp,
kegg_map_pathway_info_fp,kegg_orthology_info_fp,pfam_info_fp,
cazy_class_info_fp,cazy_family_info_fp,eggnog_info_fp,
eggnog_code_info_fp,emapper_annotations_fp)

# Read each annotation file and combine
combined_annotations <- purrr::map_dfr(
  input_annotations,
  ~ readr::read_delim(.x)  # adjust delim if needed
)
conduitR::log_with_timestamp("Converting protein IDs of annotations to the protein groups that were detected")

qf = readRDS(qf_fp)

# Getting the protein group associated with each protein annotation
protein_groups <- SummarizedExperiment::rowData(qf[["protein_groups"]]) |>
  tibble::as_tibble() |>
  dplyr::rename(Protein.Group.Temp = Protein.Group) |>
  dplyr::mutate(Protein.Group = Protein.Group.Temp)|>
  tidyr::separate_rows(Protein.Group.Temp, sep = ";") |>
  dplyr::select(Protein.Group,protein_id = Protein.Group.Temp,dplyr::everything())

#
conduit_annotations = protein_groups |>
  dplyr::left_join(combined_annotations, by = "protein_id")|>
  dplyr::select(Protein.Group,annotation_type,term,description)|>
  # If the proteinids in a protein group have the same content, they will only be counted once.
  dplyr::distinct()|>
  dplyr::filter(!is.na(term))

# Fill missing descriptions for the eggNOG-mapper-derived accession types from
# authoritative external dictionaries. Each vocabulary is tagged with the exact
# annotation_type it describes; add_term_descriptions() only fills blanks and
# never borrows a description from another source, so provenance stays explicit.
conduitR::log_with_timestamp("Filling descriptions from authoritative dictionaries")

term_dictionary <- dplyr::bind_rows(
  conduitR::parse_go_obo(go_obo_fp)        |> dplyr::mutate(annotation_type = "go"),
  conduitR::read_kegg_list(kegg_ko_fp)     |> dplyr::mutate(annotation_type = "kegg_orthology"),
  conduitR::read_kegg_list(kegg_pathway_fp)|> dplyr::mutate(annotation_type = "kegg_map_pathway"),
  conduitR::read_kegg_list(kegg_module_fp) |> dplyr::mutate(annotation_type = "kegg_module"),
  conduitR::read_kegg_list(kegg_brite_fp)  |> dplyr::mutate(annotation_type = "brite"),
  conduitR::parse_enzyme_dat(enzyme_dat_fp)|> dplyr::mutate(annotation_type = "ec_number"),
  conduitR::parse_pfam_clans(pfam_clans_fp)|> dplyr::mutate(annotation_type = "pfam")
) |>
  dplyr::select(annotation_type, term, description)

n_blank_before <- sum(is.na(conduit_annotations$description) |
                        conduit_annotations$description == "", na.rm = TRUE)

conduit_annotations <- conduitR::add_term_descriptions(conduit_annotations, term_dictionary)

n_blank_after <- sum(is.na(conduit_annotations$description) |
                       conduit_annotations$description == "", na.rm = TRUE)
conduitR::log_with_timestamp(sprintf(
  "Descriptions filled: blank rows %d -> %d of %d total",
  n_blank_before, n_blank_after, nrow(conduit_annotations)))

conduitR::log_with_timestamp("Writing conduit annotations to file")

readr::write_delim(conduit_annotations,conduit_annotations_fp)

end_time <- Sys.time()
conduitR::log_with_timestamp("Completed consolidate_annotations.R script. Time taken: %.2f minutes", 
    as.numeric(difftime(end_time, start_time, units = "mins")))
# closing clogfile connection
sink(type = "message")
sink()
close(zz)