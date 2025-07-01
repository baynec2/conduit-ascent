################################################################################
# Preparing QFeatures object from matrcies
################################################################################
# Open the log file to write both stdout and stderr
logfile <- snakemake@log[[1]]
zz <- file(logfile, open = "a")
sink(zz,append = TRUE)       # redirect stdout
sink(zz, type = "message")  # redirect stderr/messages

start_time <- Sys.time()

conduitR::log_with_timestamp("Running prepare_qf.R script")
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[1]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[2]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[3]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[4]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[5]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[6]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[7]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[8]])) 
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[9]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[10]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[11]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[12]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[14]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[15]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[16]]))
conduitR::log_with_timestamp(paste0("Input file: ", snakemake@input[[17]]))

conduitR::log_with_timestamp(paste0("Output file: ", snakemake@output[[1]]))

# Inputs
annotation = snakemake@input[["annotation"]]
precursor_matrix = snakemake@input[["precursor_matrix"]]
peptide_matrix = snakemake@input[["peptide_matrix"]]
protein_group_matrix = snakemake@input[["protein_group_matrix"]]
domain_matrix = snakemake@input[["domain_matrix"]]
kingdom_matrix= snakemake@input[["kingdom_matrix"]]
phylum_matrix= snakemake@input[["phylum_matrix"]]
class_matrix= snakemake@input[["class_matrix"]]
order_matrix= snakemake@input[["order_matrix"]]
family_matrix= snakemake@input[["family_matrix"]]
genus_matrix= snakemake@input[["genus_matrix"]]
species_matrix= snakemake@input[["species_matrix"]]
go_matrix= snakemake@input[["go_matrix"]]
go_taxa_matrix = snakemake@input[["go_taxa_matrix"]]
subcellular_locations_matrix = snakemake@input[["subcellular_locations_matrix"]]
kegg_pathway_matrix = snakemake@input[["kegg_pathway_matrix"]]
kegg_ko_matrix = snakemake@input[["kegg_ko_matrix"]]
# Output
qf = snakemake@output[["qf"]]

vector_of_matrix_fps = c(precursor_matrix,peptide_matrix,protein_group_matrix,
                         domain_matrix,kingdom_matrix,phylum_matrix,
                         class_matrix,order_matrix,family_matrix,genus_matrix,
                         species_matrix,go_matrix,go_taxa_matrix,subcellular_locations_matrix,
                         kegg_pathway_matrix,kegg_ko_matrix)


conduitR::log_with_timestamp("Preparing QFeatures object")
# Loading QFeatures object
QF = conduitR::prepare_qfeature(sample_annotation_fp = annotation,
                                         vector_of_matrix_fps)

conduitR::log_with_timestamp(paste0("Saving QFeatures object to ", qf))        

saveRDS(QF,qf)
end_time <- Sys.time()
conduitR::log_with_timestamp("prepare_qf.R script completed. Time taken: %.2f minutes", 
                              as.numeric(difftime(end_time, start_time, units = "mins")))
# closing clogfile connection
sink(type = "message")
sink()
close(zz)
