library(data.table)
library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)
library(ggbeeswarm)
library(cowplot)
library(ggsci)
library(rstatix)
options(scipen = 10000)


# source
source("plot_parameters.R")

# read in metadata
metadata <- fread("96_metadata/metadata_RNAseq_full_for_tissues.tsv")

###################################################################################################################################
############################################################ LOAD DATA ############################################################
###################################################################################################################################

##############################################################
######## Pseudoalignment (Salmon -> Wasabi -> Sleuth) ########
##############################################################
df_pseduo <- data.table(read_excel("02_tidy_data/Suppl_tables/Supplementary Tables.xlsx", sheet = 9, skip = 1, col_types = c("text", "numeric", "numeric", "text", "numeric", "numeric", "numeric", "text")))

# remove NAs and headers left while merging
df_pseduo <- df_pseduo[!is.na(df_pseduo$log2FC),]
df_pseduo <- subset(df_pseduo, target_id != "target_id")

# make short for merging (i.e. remove unwanted columns)
df_pseduo_short <- dplyr::select(df_pseduo, target_id, log2FC, tissue_id)

#####################################################################
######## reference-alignment (STAR -> featureCount -> edgeR) ########
#####################################################################
df_reference <- data.table(read_excel("02_tidy_data/Suppl_tables/Supplementary Tables.xlsx", sheet = 8, skip = 1, col_types = c("text", "numeric", "numeric", "text", "numeric", "numeric", "text")))

# make short for merging (i.e. remove unwanted columns) and melt.
df_reference_short <- melt(dplyr::select(df_reference, -start, -pval, -end, -seqnames))
df_reference_short <- df_reference_short[,c("target_id", "tissue_id", "value")]

# remove genes with several ensembl_ids, i.e. count above than 1 from this code.
excludez <- df_reference_short[,c("target_id", "tissue_id", "value")] %>% dplyr::group_by(target_id) %>% dplyr::count()
excludez <- excludez[excludez$n > 10,]

df_reference_short_filt <- df_reference_short[!df_reference_short$target_id %in% excludez$target_id,]

# change colnames to match pseudo df
colnames(df_reference_short_filt) <- c("target_id", "tissue_id", "log2FC")

# declare vector of tissues to keep
keep <- c(gsub(sort(unique(df_reference_short_filt$tissue_id)), pattern = "logFC.", replacement = ""), "Whole_Blood")

################################################################
####################### Bind together dfs ######################
################################################################

df_pseduo_short$study <- "Pseudo-alignment"
df_reference_short_filt$study <- "Reference-alignment"

df_smash <- rbind(df_pseduo_short, df_reference_short_filt)
df_smash <- df_smash[!is.na(df_smash$target_id),]

df_smash$tissue_id <- gsub(df_smash$tissue_id, pattern = "logFC.", replacement = "")

################################################################
###################### Read in annotations #####################
################################################################

# Read in the Tukiainen log2FC data, transform it to long format.
tuki_df <- fread("97_indexes_annotations/Suppl.Table.2.csv")
tuki_df[tuki_df$`target_id name` == "6-Sep"]$`target_id name` <- "SEPT6"
tuki_log2fc_df <- tuki_df %>% dplyr:: select('Gene name', ends_with("logFC"))
melt_tuki_log2fc_df <- melt(tuki_log2fc_df)
melt_tuki_log2fc_df$variable <- gsub(melt_tuki_log2fc_df$variable, pattern = "_logFC", replacement = "")
names(melt_tuki_log2fc_df)[names(melt_tuki_log2fc_df) == 'Gene name'] <- 'target_id'

#read in annotations (PAR and the Tukiainen annotation)
PAR <- rbind(fread("97_indexes_annotations/group-715_PAR1.csv"), fread("97_indexes_annotations/group-716_PAR2.csv"))
esc <- fread("97_indexes_annotations/landscape.Suppl.Table.13.csv")
esc[esc$Gene_name == "6-Sep",]$Gene_name <- "SEPT6"
esc <- esc[,-1]
esc <- esc[!(esc$Gene_name == "IDS" & esc$Reported_XCI_status == "Unknown"),]

# Read in GCF target_id annotation
annotation_gene_table <- subset(fread("97_indexes_annotations/GCF_000001405.38_GRCh38.p12_genomic_with_correct_contig_names.tsv"), type == "gene_id" & is.na(pseudo))

# merge annotations
chrx_annotation <- merge(esc, PAR, by.x = "Gene_name" , by.y = "Approved symbol" , all = T)[-1,]
chrx_annotation$category <- chrx_annotation$Reported_XCI_status
chrx_annotation[chrx_annotation$PAR == "PAR1" | chrx_annotation$PAR == "PAR2",]$category <- "PAR"
chrx_annotation[chrx_annotation$Gene_name == "XG",]$category <- "PAR"

chrx_annotation[chrx_annotation$category == "Unknown",]$category <- "Potential"

################################################################
###################### Add annotation to df ####################
################################################################

df_smash_anno_chrx <- merge(df_smash, chrx_annotation, by.x = "target_id", by.y = "Gene_name", all.x = F)
df_smash_anno_all_chrs <- merge(df_smash, chrx_annotation, by.x = "target_id", by.y = "Gene_name", all.x = T)

# small fixes
df_smash_anno_chrx$log2FC <- as.numeric(df_smash_anno_chrx$log2FC)
df_smash_anno_chrx[df_smash_anno_chrx$tissue_id == "Whole_blood",]$tissue_id <- "Whole_Blood"


################################################################
######################### Add TPMs to df #######################
################################################################

# filter metadata
TPMz <- melt(fread("02_tidy_data/tissues_TPM_matrix.tsv"))
metadata_tpm <- dplyr::select(metadata, -2)

# change "." to "-".
metadata_tpm$V2 <- gsub(metadata_tpm$`entity:sample_id`, pattern = "\\.", replacement = "-")

# add meta
tpmz_mit_meta <- merge(TPMz, metadata_tpm, by.x = "variable", by.y = "entity:sample_id")

# calculate tpms per gene and tissue.
tpmz_per_tissue <- tpmz_mit_meta %>% dplyr::group_by(V1, tissue_id) %>% dplyr::summarise(meanz = mean(value))

# merge
df_smash_anno_chrx_tpmz <- merge(df_smash_anno_chrx, tpmz_per_tissue, by.x = c("target_id", "tissue_id"), by.y = c("V1", "tissue_id"))

###################################################################################################################################
############################################################ PLOT DATA ############################################################
###################################################################################################################################

# filter metadata
metadata <- metadata[metadata$`entity:sample_id` %in% tpmz_mit_meta$variable,]

# make stats df, including only genes with a TPM > 1
df_smash_anno_chrx_stats <- df_smash_anno_chrx_tpmz[df_smash_anno_chrx_tpmz$meanz > 1,] %>% dplyr::group_by(tissue_id, study, category) %>% rstatix::get_summary_stats(log2FC, type = "common")

# quite a variance in normality, going with median and boxplots
shappurooo <- df_smash_anno_chrx %>% dplyr::group_by(study, category, tissue_id) %>% shapiro_test(log2FC)

# remove genes with a mean TPM < & cap log2FC calues for cleaner plotting.
fin_filt <- df_smash_anno_chrx_tpmz[df_smash_anno_chrx_tpmz$meanz > 1,]
fin_filt[fin_filt$log2FC > 1,]$log2FC <- 1.1
fin_filt[fin_filt$log2FC < -1,]$log2FC <- -1.1

#plot
allez <- 
  ggplot(fin_filt[fin_filt$study != "Affymetrix",], 
         aes(x=factor(category, levels = c("PAR", "Escape", "Variable", "Inactive", "Potential")),y=log2FC,col=study)) + 
  geom_quasirandom(dodge.width = .75) +  
  geom_boxplot(outlier.shape = NA) + 
  geom_hline(yintercept = 0) + 
  theme_AL_box() +
  labs(x="", y="FC (log2)", title = "All tissues")+
  theme(legend.title = element_blank(), legend.position = "top")+
  coord_cartesian(ylim=c(-1.1,1.1))+
  scale_y_continuous(breaks = c(-1,0,1))+
  geom_hline(yintercept = c(-1,1),lty=2)+
  scale_color_manual(values=c("Pseudo-alignment"="black","Reference-alignment"="#5F7D7B"))


ggsave2("04_plots/Suppl_Figure_3b.pdf", height = 8, width = 8,
        plot_grid(allez))

# export source data
write.table(fin_filt[fin_filt$study != "Affymetrix",], file = "06_source_data/Suppl_fig_3b.tsv", quote = F, row.names = F, sep = "\t")



