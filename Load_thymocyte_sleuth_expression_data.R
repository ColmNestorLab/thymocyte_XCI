source("plot_parameters.R")

#Import metadata
library(data.table)
message("loading metadata")
meta_thymo_transition <- fread("96_metadata/thymocytes_RNAsequencing_meta.csv")
meta_thymo_transition <- meta_thymo_transition[order(meta_thymo_transition$Sample),]
meta_thymo_transition <- meta_thymo_transition[meta_thymo_transition$Patient_ID != "T1"]

#Load sleuth object for differential expression
library(sleuth)
message("loading sleuth object")
so_thymo_transcripts <- sleuth_load("01_processed_data/Bulk_RNAseq_quantification/sleuth/thymocyte_transcript_mode_default_filter_no_Turner_no_ribosomal_proteins_RNAs")
so_thymo_genes <- sleuth_load("01_processed_data/Bulk_RNAseq_quantification/sleuth/thymocyte_gene_mode_default_filter_no_Turner_no_ribosomal_proteins_RNAs")

# Collapse to gene
refseq <- unique(fread("97_indexes_annotations/hg38_annotation_for_sleuth_sans_chrY_PAR_and_HGNC_ribosomal_proteins_and_RNAs.tsv",header=F))
txi_thymo_transition <- summarizeSleuthToGene(so_thymo_transcripts,refseq)

# Add metadata to sleuth object
so_thymo_transcripts$sample_to_covariates <- cbind(so_thymo_transcripts$sample_to_covariates,meta_thymo_transition[,-1])
so_thymo_genes$sample_to_covariates <- cbind(so_thymo_genes$sample_to_covariates,meta_thymo_transition[,-1])