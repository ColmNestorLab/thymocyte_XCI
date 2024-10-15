library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ggbeeswarm)
library(readxl)

# source
source("plot_parameters.R")

# Read in 
df_reference <- data.table(read_excel("Supplementary Data.xlsx", sheet = 6, skip = 1, col_types = c("text", "numeric", "numeric", "text", "numeric", "numeric", "text")))

# make short for merging (i.e. remove unwanted columns) and melt.
df_reference_short <- melt(dplyr::select(df_reference, -start, -end, -seqnames, -pval))[,c("gene", "tissue_id", "value")]

# remove genes with several ensembl_ids, i.e. count above than 1 from this code.
excludez <- df_reference_short %>% dplyr::group_by(gene) %>% dplyr::count()
excludez <- excludez[excludez$n > 10,]

df_reference_short_filt <- df_reference_short[!df_reference_short$gene %in% excludez$gene,]

# change colnames
colnames(df_reference_short_filt) <- c("target_id", "tissue_id", "b")

# declare vector of tissues to keep
keep <- c(gsub(sort(unique(df_reference_short_filt$tissue_id)), pattern = "logFC.", replacement = ""), "Whole_Blood")

# read in metadata
metadata <- fread("96_metadata/METADATA_RNAseq_affymetrix_matched.tsv")

################################################################
######## reference-alignment (affymetrix, preprocessed) ########
################################################################
df_affymetrix <- fread("02_tidy_data/GTEX/method_comparison/affymetrix/GSE45878_MA_data.txt", header = T)

# melt and change colnames
melt_df_affymetrix <- melt(df_affymetrix)
colnames(melt_df_affymetrix) <- c("id_ref", "sample", "val")

# read in affymetrix metadata
meta_affy <- dplyr::select(fread("96_metadata/metadata_GSE45878.tsv", header = F), V1, V4, V5, V8)

# change tissue names
meta_affy[meta_affy$V8 == "Adipose - Subcutaneous"]$V8 <- "Adipose_Subcutaneous"
meta_affy[meta_affy$V8 == "Artery - Tibial"]$V8 <- "Artery_Tibial"
meta_affy[meta_affy$V8 == "Brain - Cortex"]$V8 <- "Brain_Cortex"
meta_affy[meta_affy$V8 == "Heart - Left Ventricle"]$V8 <- "Heart_Left_Ventricle"
meta_affy[meta_affy$V8 == "Lung"]$V8 <- "Lung"
meta_affy[meta_affy$V8 == "Muscle - Skeletal"]$V8 <- "Muscle_Skeletal"
meta_affy[meta_affy$V8 == "Nerve - Tibial"]$V8 <- "Nerve_Tibial"
meta_affy[meta_affy$V8 == "Skin - Sun Exposed (Lower leg)"]$V8 <- "Skin_Sun_Exposed_Lower_leg"
meta_affy[meta_affy$V8 == "Whole Blood"]$V8 <- "Whole_Blood"

# merge metadata and data, change colnames
melt_df_affymetrix <- merge(melt_df_affymetrix, meta_affy, by.x = c("sample"), by.y = "V1", all = T)
colnames(melt_df_affymetrix) <- c("sample", "id_ref", "val", "sex", "age", "tissue_id")

# split
melt_df_affymetrix <- melt_df_affymetrix %>% separate(id_ref, into = c("genes", NA), sep = "_", remove = F)

# convert ensembl IDs to gene symbols.
ensembl_id_to_symbol <- fread("97_indexes_annotations/ensembl_to_symbol.tsv", header = F)
colnames(ensembl_id_to_symbol) <- c("enseml_id", "gene")
ensembl_id_to_symbol_fix <- ensembl_id_to_symbol %>% separate(enseml_id, into = c("genes", NA), sep = "\\.", remove = F)
melt_df_affymetrix <- merge(melt_df_affymetrix, ensembl_id_to_symbol_fix, by.x =c("genes"), by.y = "genes", all.x =T)

# change colnames
colnames(melt_df_affymetrix) <- c("ensembl_short", "sample", "id_ref", "val", "sex", "age", "tissue_id","ensembl_id", "target_id")

# remove excess columns
melt_df_affymetrix_short <- dplyr::select(melt_df_affymetrix, target_id, sample, val, sex, tissue_id)

# get means
transformation_before_log2FC_calc <- melt_df_affymetrix_short %>% dplyr::group_by(target_id, sex, tissue_id) %>% dplyr::summarise(meanz = mean(val))

# dcast
sex_split <- dcast(transformation_before_log2FC_calc[transformation_before_log2FC_calc$tissue_id %in% keep,], formula = target_id + tissue_id ~ sex, value.var = "meanz")

# calculate log2FC
sex_split$b <- log2(sex_split$Female / sex_split$Male)

# remove excess columns
df_MA_sex_short <- dplyr::select(sex_split, -Female, -Male)

# Read in the Tukiainen log2FC data, transform it to long format.
tuki_df <- fread("97_indexes_annotations/Suppl.Table.2.csv")
tuki_df[tuki_df$`Gene name` == "6-Sep"]$`Gene name` <- "SEPT6"
tuki_log2fc_df <- tuki_df %>% dplyr:: select('Gene name', ends_with("logFC"))
melt_tuki_log2fc_df <- melt(tuki_log2fc_df)
melt_tuki_log2fc_df$variable <- gsub(melt_tuki_log2fc_df$variable, pattern = "_logFC", replacement = "")
names(melt_tuki_log2fc_df)[names(melt_tuki_log2fc_df) == 'Gene name'] <- 'gene'

#read in annotations (PAR and the Tukiainen annotation)
PAR <- rbind(fread("97_indexes_annotations/group-715_PAR1.csv"), fread("97_indexes_annotations/group-716_PAR2.csv"))
esc <- fread("97_indexes_annotations/landscape.Suppl.Table.13.csv")
esc[esc$Gene_name == "6-Sep",]$Gene_name <- "SEPT6"
esc <- esc[,-1]
esc <- esc[!(esc$Gene_name == "IDS" & esc$Reported_XCI_status == "Unknown"),]

# Read in GCF gene annotation
annotation_gene_table <- subset(fread("97_indexes_annotations/GCF_000001405.38_GRCh38.p12_genomic_with_correct_contig_names.tsv"), type == "gene" & is.na(pseudo))

# merge annotations
chrx_annotation <- merge(esc, PAR, by.x = "Gene_name" , by.y = "Approved symbol" , all = T)[-1,]
chrx_annotation$category <- chrx_annotation$Reported_XCI_status
chrx_annotation[chrx_annotation$PAR == "PAR1" | chrx_annotation$PAR == "PAR2",]$category <- "PAR"
chrx_annotation[chrx_annotation$Gene_name == "XG",]$category <- "PAR"

chrx_annotation[chrx_annotation$category == "Unknown",]$category <- "Potential"
chrx_annotation[chrx_annotation$category == "Variable",]$category <- "Potential"

################################################################
###################### Add annotation to df ####################
################################################################

# Merge
df_smash_anno_chrx <- merge(df_MA_sex_short, chrx_annotation, by.x = "target_id", by.y = "Gene_name", all.x = F)

# small fixes
df_smash_anno_chrx$b <- as.numeric(df_smash_anno_chrx$b)

# rename tissue names
df_smash_anno_chrx[df_smash_anno_chrx$tissue_id == "Nerve_Tibial",]$tissue_id <- "nerve"
df_smash_anno_chrx[df_smash_anno_chrx$tissue_id == "Artery_Tibial",]$tissue_id <- "artery-3"
df_smash_anno_chrx[df_smash_anno_chrx$tissue_id == "Heart_Left_Ventricle",]$tissue_id <- "heart-2"
df_smash_anno_chrx[df_smash_anno_chrx$tissue_id == "Whole_Blood",]$tissue_id <- "whole blood"
df_smash_anno_chrx[df_smash_anno_chrx$tissue_id == "Skin_Sun_Exposed_Lower_leg",]$tissue_id <- "skin-2"
df_smash_anno_chrx[df_smash_anno_chrx$tissue_id == "Adipose_Subcutaneous",]$tissue_id <- "adipose-1"
df_smash_anno_chrx[df_smash_anno_chrx$tissue_id == "Brain_Cortex",]$tissue_id <- "brain-6"
df_smash_anno_chrx[df_smash_anno_chrx$tissue_id == "Muscle_Skeletal",]$tissue_id <- "muscle"
df_smash_anno_chrx[df_smash_anno_chrx$tissue_id == "Lung",]$tissue_id <- "lung"

# set order for plotting
orderz <- df_smash_anno_chrx %>% dplyr::group_by(tissue_id, category) %>% dplyr::summarise(medianz = median(b))
PAR_orderz <- orderz[orderz$category == "PAR",]
PAR_orderz <- PAR_orderz[order(PAR_orderz$medianz, decreasing = T),]
Escape_orderz <- orderz[orderz$category == "Escape",]
Escape_orderz <- Escape_orderz[order(Escape_orderz$medianz, decreasing = T),]
Inactive_orderz <- orderz[orderz$category == "Inactive",]
Inactive_orderz <- Inactive_orderz[order(Inactive_orderz$medianz, decreasing = T),]

# plot
PAR <- 
  ggplot(df_smash_anno_chrx[df_smash_anno_chrx$category %in% c("PAR"),], aes(x=factor(tissue_id, levels = PAR_orderz$tissue_id), y=b, col = tissue_id)) + 
  geom_quasirandom(show.legend = F, alpha = .5)+
  geom_boxplot(show.legend = F, outlier.shape = NA) + 
  geom_hline(yintercept = 0, lty=2)+
  theme_AL_box_rotX() +
  scale_y_continuous(breaks=c(-0.5,0,0.5), limits = c(-0.51,0.51))+
  scale_color_manual(values = color_valzz)+
  labs(x="", title = "PAR")

escapez <- 
  ggplot(df_smash_anno_chrx[df_smash_anno_chrx$category %in% c("Escape"),], aes(x=factor(tissue_id, levels = Escape_orderz$tissue_id), y=b, col = tissue_id)) + 
  geom_quasirandom(show.legend = F, alpha = .5)+
  geom_boxplot(show.legend = F, outlier.shape = NA) + 
  geom_hline(yintercept = 0, lty=2)+
  theme_AL_box_rotX() +
  scale_y_continuous(breaks=c(-0.5,0,0.5), limits = c(-0.51,0.51))+
  scale_color_manual(values = color_valzz)+
  labs(x="", title = "Escape")

inacz <- 
  ggplot(df_smash_anno_chrx[df_smash_anno_chrx$category %in% c("Inactive"),], aes(x=factor(tissue_id, levels = Inactive_orderz$tissue_id), y=b, col = tissue_id)) + 
  geom_quasirandom(show.legend = F, alpha = .5)+
  geom_boxplot(show.legend = F, outlier.shape = NA) + 
  geom_hline(yintercept = 0, lty=2)+
  theme_AL_box_rotX() +
  scale_y_continuous(breaks=c(-0.5,0,0.5), limits = c(-0.51,0.51))+
  scale_color_manual(values = color_valzz)+
  labs(x="", title = "Inactive")

ggsave2("04_plots/Suppl_Figure_3e.pdf", width = 16, height = 6,
        plot_grid(PAR, escapez, inacz, ncol = 3)
)


#Export source data
exportz <- df_smash_anno_chrx[df_smash_anno_chrx$category %in% c("PAR", "Escape", "Inactive"),c("target_id",   "tissue_id"      ,      "b", "Chr", "category")]
colnames(exportz) <- c("target_id",   "tissue_id"      ,      "b", "Chr", "category_Tuki_but_with_PAR_separate")
write.table(exportz, file = "06_source_data/Suppl_Fig_3e.tsv", quote = F, row.names = F, sep = "\t")

