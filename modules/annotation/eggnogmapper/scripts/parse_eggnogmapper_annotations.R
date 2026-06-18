################################################################################
# Parse eggNOG-mapper Annotations
################################################################################
# Parses emapper.annotations into the standard long-format annotation table:
#   (protein_id, annotation_type, term, description)
#
# Covers the following annotation types:
#   go               <- GOs column
#   kegg_orthology   <- KEGG_ko column (strips "ko:" prefix)
#   kegg_map_pathway <- KEGG_Pathway column (already "map" prefixed)
#   eggnog           <- eggNOG_OGs column (OG ID before "@")
#   eggnog_code      <- COG_category column (individual letter codes)
#   pfam             <- PFAMs column
#   cazy_class       <- CAZy column (class prefix, e.g. "GH" from "GH1")
#   cazy_family      <- CAZy column (full family ID, e.g. "GH1")
#   gene_symbol      <- Preferred_name column (symbol is its own description)
#   ec_number        <- EC column (descriptions filled from ENZYME dictionary)
#   kegg_module      <- KEGG_Module column (descriptions filled from KEGG list)
#   brite            <- BRITE column (descriptions filled from KEGG list)
#
# Descriptions for go/kegg_*/pfam/ec_number/brite are intentionally left NA here
# and filled in consolidate_annotations.R from authoritative offline dictionaries
# (see conduitR::add_term_descriptions). KEGG_Reaction, KEGG_rclass, KEGG_TC and
# BiGG_Reaction columns are intentionally not parsed (low fill / niche).
#
# eggNOG-mapper v2 output format notes:
#   - Tab-separated; ## lines are metadata comments; column header starts with #
#   - Missing values represented as "-"
#   - Multiple values within a column are comma-separated
#   - eggNOG_OGs format: "OG@tax_id|tax_name,..." (e.g. "COG0001@2|Bacteria")
#   - KEGG_ko format: "ko:K00001,ko:K00002"
#   - KEGG_Pathway format: "map00010,map00020"
#   - COG_category: one or more concatenated letters (e.g. "C", "CE", "GEK")
################################################################################

# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz, append = TRUE)
sink(zz, type = "message")

start_time <- Sys.time()

conduitR::log_with_timestamp("Running parse_eggnogmapper_annotations.R script")
conduitR::log_with_timestamp(paste0("Input file:  ", snakemake@input[["annotations"]]))
conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[["emapper_annotations"]]))

annotations_fp        <- snakemake@input[["annotations"]]
emapper_annotations_fp <- snakemake@output[["emapper_annotations"]]

################################################################################
# Read emapper output
# "##" lines are metadata comments; the column-header line begins with "#"
# so readr will name the first column "#query_name" — rename to protein_id.
################################################################################
conduitR::log_with_timestamp("Reading emapper annotations file")

raw <- readr::read_tsv(
  annotations_fp,
  comment   = "##",
  col_names = TRUE,
  name_repair = "minimal",
  show_col_types = FALSE
)

colnames(raw)[1] <- "protein_id"

# Normalize protein IDs: strip leading "sp|" / "tr|" and trailing "|NAME_SPECIES"
# e.g. "tr|A0A1D5NW85|A0A1D5NW85_CHICK" -> "A0A1D5NW85"
raw <- raw |>
  dplyr::mutate(protein_id = stringr::str_remove(protein_id, "^[a-z]+\\|") |>
                               stringr::str_remove("\\|.*$"))

# Replace the "-" sentinel with NA throughout
raw <- raw |>
  dplyr::mutate(dplyr::across(where(is.character), ~ dplyr::na_if(.x, "-")))

conduitR::log_with_timestamp(paste0("Read ", nrow(raw), " protein entries"))

################################################################################
# Helper: expand a comma-delimited column into long format
################################################################################
expand_col <- function(df, col) {
  df |>
    dplyr::select(protein_id, value = {{ col }}) |>
    dplyr::filter(!is.na(value)) |>
    tidyr::separate_longer_delim(value, delim = ",") |>
    dplyr::mutate(value = stringr::str_trim(value)) |>
    dplyr::filter(value != "", !is.na(value)) |>
    dplyr::distinct()
}

################################################################################
# COG category lookup table
################################################################################
eggnog_code_lookup <- tibble::tibble(
  code = c("A","B","C","D","E","F","G","H","I","J",
           "K","L","M","N","O","P","Q","R","S","T",
           "U","V","W","Y","Z"),
  description = c(
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

################################################################################
# CAZy class lookup table
################################################################################
cazy_class_lookup <- tibble::tibble(
  class       = c("GH",  "GT",                 "PL",                  "CE",                    "AA",               "CBM"),
  description = c("glycoside_hydrolase", "glycosyl_transferase", "polysaccharide_lyase",
                  "carbohydrate_esterase", "auxiliary_activity", "carbohydrate_binding_module")
)

################################################################################
# 1. GO terms
################################################################################
conduitR::log_with_timestamp("Extracting GO annotations")

go_annotations <- expand_col(raw, GOs) |>
  dplyr::filter(stringr::str_starts(value, "GO:")) |>
  dplyr::mutate(annotation_type = "go", description = NA_character_) |>
  dplyr::select(protein_id, annotation_type, term = value, description)

################################################################################
# 2. KEGG orthology
# Format: "ko:K00001,ko:K00002" — strip the "ko:" prefix
################################################################################
conduitR::log_with_timestamp("Extracting KEGG orthology annotations")

kegg_orthology_annotations <- expand_col(raw, KEGG_ko) |>
  dplyr::mutate(value = stringr::str_remove(value, "^ko:")) |>
  dplyr::filter(stringr::str_starts(value, "K")) |>
  dplyr::mutate(annotation_type = "kegg_orthology", description = NA_character_) |>
  dplyr::select(protein_id, annotation_type, term = value, description) |>
  dplyr::distinct()

################################################################################
# 3. KEGG map pathways
# Format: "map00010,map00020" — already standardised to "map" prefix
################################################################################
conduitR::log_with_timestamp("Extracting KEGG pathway annotations")

kegg_pathway_annotations <- expand_col(raw, KEGG_Pathway) |>
  dplyr::filter(stringr::str_detect(value, "^map\\d{5}$")) |>
  dplyr::mutate(annotation_type = "kegg_map_pathway", description = NA_character_) |>
  dplyr::select(protein_id, annotation_type, term = value, description) |>
  dplyr::distinct()

################################################################################
# 4. eggNOG orthologous groups
# Format: "COG0001@1|root,COG0001@2|Bacteria,..."
# Extract the OG identifier (everything before "@").
# Use the emapper Description column as description.
################################################################################
conduitR::log_with_timestamp("Extracting eggNOG OG annotations")

protein_descriptions <- raw |>
  dplyr::select(protein_id, emapper_description = Description) |>
  dplyr::filter(!is.na(emapper_description))

eggnog_annotations <- expand_col(raw, eggNOG_OGs) |>
  dplyr::mutate(value = stringr::str_extract(value, "^[^@]+")) |>
  dplyr::filter(!is.na(value)) |>
  dplyr::left_join(protein_descriptions, by = "protein_id") |>
  dplyr::mutate(annotation_type = "eggnog") |>
  dplyr::select(protein_id, annotation_type, term = value, description = emapper_description) |>
  dplyr::distinct()

################################################################################
# 5. COG / eggNOG functional categories
# COG_category can be one or more concatenated letters (e.g. "C", "CE", "GEK")
################################################################################
conduitR::log_with_timestamp("Extracting eggNOG code annotations")

eggnog_code_annotations <- raw |>
  dplyr::select(protein_id, COG_category) |>
  dplyr::filter(!is.na(COG_category)) |>
  dplyr::mutate(code = stringr::str_split(COG_category, "")) |>
  tidyr::unnest(code) |>
  dplyr::inner_join(eggnog_code_lookup, by = "code") |>
  dplyr::mutate(annotation_type = "eggnog_code") |>
  dplyr::select(protein_id, annotation_type, term = code, description) |>
  dplyr::distinct()

################################################################################
# 6. Pfam domains
################################################################################
conduitR::log_with_timestamp("Extracting Pfam annotations")

pfam_annotations <- expand_col(raw, PFAMs) |>
  dplyr::filter(stringr::str_starts(value, "PF")) |>
  dplyr::mutate(annotation_type = "pfam", description = NA_character_) |>
  dplyr::select(protein_id, annotation_type, term = value, description) |>
  dplyr::distinct()

################################################################################
# 7. CAZy class and family
# Standard CAZy families are prefixed by class code (GH, GT, PL, CE, AA, CBM)
# followed by a number (e.g. "GH1", "GT2").
################################################################################
conduitR::log_with_timestamp("Extracting CAZy annotations")

cazy_expanded <- expand_col(raw, CAZy) |>
  # Extract leading alpha class prefix
  dplyr::mutate(class = stringr::str_extract(value, "^[A-Za-z]+"))

cazy_class_annotations <- cazy_expanded |>
  dplyr::inner_join(cazy_class_lookup, by = "class") |>
  dplyr::mutate(annotation_type = "cazy_class") |>
  dplyr::select(protein_id, annotation_type, term = class, description) |>
  dplyr::distinct()

cazy_family_annotations <- cazy_expanded |>
  # inner_join restricts to known classes and attaches the class-level
  # description (e.g. "GH1" -> "glycoside_hydrolase"); this is a static,
  # authoritative lookup, so it is consistent with the dictionary-only policy.
  dplyr::inner_join(cazy_class_lookup, by = "class") |>
  dplyr::mutate(annotation_type = "cazy_family") |>
  dplyr::select(protein_id, annotation_type, term = value, description) |>
  dplyr::distinct()

################################################################################
# 8. Gene symbol (Preferred_name)
# Self-describing: the symbol is its own description (uniform readout), so the
# term is mirrored into the description column.
################################################################################
conduitR::log_with_timestamp("Extracting gene symbol annotations")

gene_symbol_annotations <- expand_col(raw, Preferred_name) |>
  dplyr::mutate(annotation_type = "gene_symbol", description = value) |>
  dplyr::select(protein_id, annotation_type, term = value, description) |>
  dplyr::distinct()

################################################################################
# 9. EC numbers
# Format: "1.1.1.1,2.7.7.7" (partial ECs like "1.3.-.-" are kept).
# Descriptions are filled downstream from the ENZYME dictionary.
################################################################################
conduitR::log_with_timestamp("Extracting EC number annotations")

ec_annotations <- expand_col(raw, EC) |>
  dplyr::filter(stringr::str_detect(value, "^[0-9]")) |>
  dplyr::mutate(annotation_type = "ec_number", description = NA_character_) |>
  dplyr::select(protein_id, annotation_type, term = value, description) |>
  dplyr::distinct()

################################################################################
# 10. KEGG modules
# Format: "M00087,M00131". Descriptions filled downstream from KEGG list/module.
################################################################################
conduitR::log_with_timestamp("Extracting KEGG module annotations")

kegg_module_annotations <- expand_col(raw, KEGG_Module) |>
  dplyr::filter(stringr::str_starts(value, "M")) |>
  dplyr::mutate(annotation_type = "kegg_module", description = NA_character_) |>
  dplyr::select(protein_id, annotation_type, term = value, description) |>
  dplyr::distinct()

################################################################################
# 11. BRITE functional hierarchies
# Format: "ko00000,ko04147,br08901". Descriptions filled downstream from
# KEGG list/brite.
################################################################################
conduitR::log_with_timestamp("Extracting BRITE annotations")

brite_annotations <- expand_col(raw, BRITE) |>
  dplyr::filter(stringr::str_detect(value, "^(ko|br)")) |>
  dplyr::mutate(annotation_type = "brite", description = NA_character_) |>
  dplyr::select(protein_id, annotation_type, term = value, description) |>
  dplyr::distinct()

# NOTE: KEGG_Reaction, KEGG_rclass, KEGG_TC and BiGG_Reaction columns are
# intentionally not parsed (low fill rate / niche use); add them here if needed.

################################################################################
# Combine all annotation types and write output
################################################################################
conduitR::log_with_timestamp("Combining all annotation types")

emapper_annotations <- dplyr::bind_rows(
  go_annotations,
  kegg_orthology_annotations,
  kegg_pathway_annotations,
  eggnog_annotations,
  eggnog_code_annotations,
  pfam_annotations,
  cazy_class_annotations,
  cazy_family_annotations,
  gene_symbol_annotations,
  ec_annotations,
  kegg_module_annotations,
  brite_annotations
) |>
  dplyr::filter(!is.na(term)) |>
  dplyr::distinct()

conduitR::log_with_timestamp(paste0(
  "Total annotations: ", nrow(emapper_annotations),
  " across ", dplyr::n_distinct(emapper_annotations$annotation_type), " types"
))

conduitR::log_with_timestamp(paste0("Writing to ", emapper_annotations_fp))
readr::write_delim(emapper_annotations, emapper_annotations_fp)

end_time <- Sys.time()
conduitR::log_with_timestamp(
  "Completed parse_eggnogmapper_annotations.R script. Time taken: %.2f minutes",
  as.numeric(difftime(end_time, start_time, units = "mins"))
)

sink(type = "message")
sink()
close(zz)
