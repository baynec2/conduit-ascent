################################################################################
# Adding Protein Group Annotations to QFeatures object
################################################################################
# These annotations include taxonomy and functional annotations that were 
# obtained from uniprot or NCBI

# Some annotations rely on external references. These are contained in the 
# annotaions conduit slot (since multiple annotations can go to a single protein.)

# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

library(SummarizedExperiment)
library(QFeatures)

start_time <- Sys.time()
conduitR::log_with_timestamp("Running add_annotations_to_qfeatures.R script")

conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[["qf"]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[["uniprot_annotated_protein_info"]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[["annotated_qf"]]))

#Defining files
# Inputs
qf_fp = snakemake@input[["qf"]]
uniprot_annotated_protein_info_fp = snakemake@input[["uniprot_annotated_protein_info"]]
conduit_annotations_fp = snakemake@input[["conduit_annotations"]]

# Outputs
annotated_qf_fp = snakemake@output[["annotated_qf"]]

# Reading in files
conduitR::log_with_timestamp("Reading in files")
uniprot_annotated_protein_info = readr::read_delim(uniprot_annotated_protein_info_fp)
qf = readRDS(qf_fp)
conduitR::log_with_timestamp("Working with annotations contained in the Uniprot_annotated_protein_info_file")

# Adding taxonomy annotations
conduitR::log_with_timestamp("Adding taxonomy information to QFeatures, handling assay links, and summarizing")
qf = conduitR::add_taxonomy_to_qf(qf,uniprot_annotated_protein_info)

# Adding all of the annotations that we have extracted. 
# First, we need to pivot these to wide format. 
conduitR::log_with_timestamp("Transforming annotations contained in conduit_annotations.txt file into wide format")

# Pivoting to the proper format.
conduit_annotations_wide = readr::read_delim(conduit_annotations_fp) |>
  dplyr::select(Protein.Group, annotation_type, term) |>  # keep columns of interest
  # pivot so each annotation_type becomes a column
  tidyr::pivot_wider(
    names_from = annotation_type,
    values_from = term,
    values_fn = \(x) paste(unique(x), collapse = ";")  # collapse multiple terms per protein
  )

  # Helper: only add annotation if the column is present in conduit_annotations_wide.
  # Annotation types are absent when no proteins in the dataset have that annotation.
  maybe_add_annotation <- function(qf, column_name, ...) {
    col_str <- rlang::as_string(rlang::ensym(column_name))
    if (col_str %in% colnames(conduit_annotations_wide)) {
      conduitR::add_annotation_to_qf(qf, ..., column_name = !!rlang::sym(col_str))
    } else {
      conduitR::log_with_timestamp("Skipping %s — no annotations present for this dataset", col_str)
      qf
    }
  }

  # --- UniProt-derived annotations ---
conduitR::log_with_timestamp("Adding uniprot_go annotations to QFeatures")
qf = maybe_add_annotation(qf,
                                    id_column = Protein.Group,
                                    conduit_annotations = conduit_annotations_wide,
                                    column_name = uniprot_go)

conduitR::log_with_timestamp("Adding uniprot_pfam annotations to QFeatures")
qf = maybe_add_annotation(qf,
                                    id_column = Protein.Group,
                                    conduit_annotations = conduit_annotations_wide,
                                    column_name = uniprot_pfam)

conduitR::log_with_timestamp("Adding uniprot_eggnog annotations to QFeatures")
qf = maybe_add_annotation(qf,
                                    id_column = Protein.Group,
                                    conduit_annotations = conduit_annotations_wide,
                                    column_name = uniprot_eggnog)

conduitR::log_with_timestamp("Adding uniprot_eggnog_code annotations to QFeatures")
qf = maybe_add_annotation(qf,
                                    id_column = Protein.Group,
                                    conduit_annotations = conduit_annotations_wide,
                                    column_name = uniprot_eggnog_code)

conduitR::log_with_timestamp("Adding uniprot_kegg_pathway annotations to QFeatures")
qf = maybe_add_annotation(qf,
                                    id_column = Protein.Group,
                                    conduit_annotations = conduit_annotations_wide,
                                    column_name = uniprot_kegg_pathway)

conduitR::log_with_timestamp("Adding uniprot_kegg_map_pathway annotations to QFeatures")
qf = maybe_add_annotation(qf,
                                    id_column = Protein.Group,
                                    conduit_annotations = conduit_annotations_wide,
                                    column_name = uniprot_kegg_map_pathway)

conduitR::log_with_timestamp("Adding uniprot_kegg_orthology annotations to QFeatures")
qf = maybe_add_annotation(qf,
                                    id_column = Protein.Group,
                                    conduit_annotations = conduit_annotations_wide,
                                    column_name = uniprot_kegg_orthology)

conduitR::log_with_timestamp("Adding uniprot_cazy_class annotations to QFeatures")
qf = maybe_add_annotation(qf,
                                    id_column = Protein.Group,
                                    conduit_annotations = conduit_annotations_wide,
                                    column_name = uniprot_cazy_class)

conduitR::log_with_timestamp("Adding uniprot_cazy_family annotations to QFeatures")
qf = maybe_add_annotation(qf,
                                    id_column = Protein.Group,
                                    conduit_annotations = conduit_annotations_wide,
                                    column_name = uniprot_cazy_family)

# --- eggNOG-mapper-derived annotations ---
conduitR::log_with_timestamp("Adding go annotations to QFeatures")
qf = maybe_add_annotation(qf,
                                    id_column = Protein.Group,
                                    conduit_annotations = conduit_annotations_wide,
                                    column_name = go)

conduitR::log_with_timestamp("Adding pfam annotations to QFeatures")
qf = maybe_add_annotation(qf,
                                    id_column = Protein.Group,
                                    conduit_annotations = conduit_annotations_wide,
                                    column_name = pfam)

conduitR::log_with_timestamp("Adding eggnog annotations to QFeatures")
qf = maybe_add_annotation(qf,
                                    id_column = Protein.Group,
                                    conduit_annotations = conduit_annotations_wide,
                                    column_name = eggnog)

conduitR::log_with_timestamp("Adding eggnog_code annotations to QFeatures")
qf = maybe_add_annotation(qf,
                                    id_column = Protein.Group,
                                    conduit_annotations = conduit_annotations_wide,
                                    column_name = eggnog_code)

conduitR::log_with_timestamp("Adding kegg_map_pathway annotations to QFeatures")
qf = maybe_add_annotation(qf,
                                    id_column = Protein.Group,
                                    conduit_annotations = conduit_annotations_wide,
                                    column_name = kegg_map_pathway)

conduitR::log_with_timestamp("Adding kegg_orthology annotations to QFeatures")
qf = maybe_add_annotation(qf,
                                    id_column = Protein.Group,
                                    conduit_annotations = conduit_annotations_wide,
                                    column_name = kegg_orthology)

conduitR::log_with_timestamp("Adding cazy_class annotations to QFeatures")
qf = maybe_add_annotation(qf,
                                    id_column = Protein.Group,
                                    conduit_annotations = conduit_annotations_wide,
                                    column_name = cazy_class)

conduitR::log_with_timestamp("Adding cazy_family annotations to QFeatures")
qf = maybe_add_annotation(qf,
                                    id_column = Protein.Group,
                                    conduit_annotations = conduit_annotations_wide,
                                    column_name = cazy_family)

conduitR::log_with_timestamp(paste0("Annotations sucessfully added. Writing Qfeatures object to ", annotated_qf_fp))

saveRDS(qf,annotated_qf_fp)
end_time <- Sys.time()
conduitR::log_with_timestamp("Completed add_annotations_to_qfeatures script. Time taken: %.2f minutes", 
    as.numeric(difftime(end_time, start_time, units = "mins")))

# closing clogfile connection
sink(type = "message")
sink()
close(zz)