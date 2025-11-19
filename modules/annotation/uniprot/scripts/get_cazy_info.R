################################################################################
# Getting Cazyme Annotations for each Detected Protein ID
################################################################################
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

start_time <- Sys.time()

conduitR::log_with_timestamp("Running get_eggnog_info.R script")
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[["uniprot_annotated_protein_info"]]))
conduitR::log_with_timestamp(paste0("Output file: ",snakemake@output[["cazy_class_info"]]))
conduitR::log_with_timestamp(paste0("Output file: ",snakemake@output[["cazy_family_info"]]))

# Defining Input Files
uniprot_annotated_protein_info_fp = snakemake@input[["uniprot_annotated_protein_info"]]
cazy_resource_fp = snakemake@input[["cazy_resource"]]

# Defining Output Files
cazy_class_info_fp = snakemake@output[["cazy_class_info"]]
cazy_family_info_fp = snakemake@output[["cazy_family_info"]]

conduitR::log_with_timestamp("Reading in detected cazyme annotations")
  
detected_cazy = readr::read_delim(uniprot_annotated_protein_info_fp) |>
  dplyr::select(protein_id,xref_cazy)|>
  dplyr::filter(!is.na(xref_cazy))|>
  tidyr::separate_longer_delim(xref_cazy,delim = ";")|>
  dplyr::filter(xref_cazy != "")|>
  dplyr::mutate(class = gsub("[0-9]*","",xref_cazy))

# cazy class info
conduitR::log_with_timestamp("Exracting cazy class annotations")

cazy_class_lookup = tibble::tibble(class = c("GH","GT","PL","CE","AA","CBM"),
description = c("glycoside_hydrolase","glycosyl_transferase", "polysaccharide_lyase", 
"carbohydrate_esterase", "auxiliary_activity","carbohydrate_binding_module")
)

cazy_class_info = dplyr::inner_join(detected_cazy,cazy_class_lookup, by = "class") |>
dplyr::mutate(annotation_type = "cazy_class") |>
dplyr::select(protein_id,annotation_type,term = class,description)

conduitR::log_with_timestamp("Exracting cazy family annotations")

cazy_resource = readr::read_delim(cazy_resource_fp,skip = 2,col_names = c("xref_cazy","description"))

cazy_family_info = dplyr::inner_join(detected_cazy,cazy_resource, by = "xref_cazy") |>
dplyr::mutate(annotation_type = "cazy_family")|>
dplyr::select(protein_id,annotation_type,term = xref_cazy,description)

conduitR::log_with_timestamp("Saving cazy class annotations to file")
readr::write_delim(cazy_class_info, cazy_class_info_fp)

conduitR::log_with_timestamp("Saving cazy family annotations to file")
readr::write_delim(cazy_family_info,cazy_family_info_fp)

end_time <- Sys.time()
conduitR::log_with_timestamp("Completed get_cazy_info.R script. Time taken: %.2f minutes", 
    as.numeric(difftime(end_time, start_time, units = "mins")))
# closing clogfile connection
sink(type = "message")
sink()
close(zz)