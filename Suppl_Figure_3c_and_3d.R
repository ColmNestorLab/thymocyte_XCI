library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggsci)
library(ggbeeswarm)


source("plot_parameters.R")


#read in annotations (PAR and the Tukiainen annotation)
PAR <- rbind(fread("97_indexes_annotations/group-715_PAR1.csv"), fread("97_indexes_annotations/group-716_PAR2.csv"))
esc <- fread("97_indexes_annotations/landscape.Suppl.Table.13.csv")
esc[esc$Gene_name == "6-Sep",]$Gene_name <- "SEPT6"
esc <- esc[,-1]

# merge annotations
chrx_annotation <- merge(esc, PAR, by.x = "Gene_name" , by.y = "Approved symbol" , all = T)[-1,]
chrx_annotation$category <- chrx_annotation$Reported_XCI_status
chrx_annotation[chrx_annotation$PAR == "PAR1" | chrx_annotation$PAR == "PAR2",]$category <- "PAR"

# read in RNA-seq meta
meta_full <- fread("96_metadata/metadata_RNAseq_full_for_tissues.tsv")

# read in the 10vs10 TPM matrix
df_TPM_tissues <- dplyr::select(melt(fread("02_tidy_data/GTEX/tissue_comparison/10vs10_TPM_matrix.tsv")), -V1)

# keep only meta for samples we have included in the 10vs10.
meta_tissues <- subset(meta_full, meta_full$`entity:sample_id` %in% gsub(unique(df_TPM_tissues$variable), pattern = "\\.", replacement = "-"))

# rename tissue
meta_tissues[meta_tissues$tissue_id %in% c("Cells_EBV-transformed_lymphocytes"),]$tissue_id <- "Cells_EBV_transformed_lymphocytes"
meta_tissues[meta_tissues$tissue_id == "Brain_Spinal_cord_cervical_c-1",]$tissue_id <- "Brain_Spinal_cord_cervical_c1"

# change "-" to "."
meta_tissues$sample <- gsub(meta_tissues$`entity:sample_id`, pattern = "-", replacement = "\\.")

# read in the sleuth output after leveraging sex across all the tissues.
df_tissues <- fread("02_tidy_data/Suppl_tables/Supplementary_Data_Sleuth_table_tissue_comparison.tsv")

# make sure b (log2FC) is numeric.
df_tissues$b <- as.numeric(df_tissues$b)

#merge
df_TPM_tissues_anno <- merge(df_TPM_tissues, meta_tissues, by.x = "variable", by.y = "sample")

# summarise expression (TPM) of genes per tissue
df_TPM_tissues_anno_meanz <- df_TPM_tissues_anno %>% dplyr::group_by(tissue_id, gene) %>% dplyr::summarise(meanz = mean(value))


# read in thymoycte sleuth output
thymocytes_across <- data.table(read_excel("Supplementary Data.xlsx", sheet = 4, skip = 1, col_types = c("text","text","numeric","numeric","numeric")))

# read in thymocyte TPMs
tpmz_thymo <- melt(data.table(read_excel("Supplementary Data.xlsx", sheet = 2, skip = 1, col_types = c("text",rep(paste0("numeric"), 46)))))

# get mean expression and add tissue_id
tpmz_thymo_meanz <- tpmz_thymo %>% dplyr::group_by(gene) %>% dplyr::summarise(meanz = mean(value))
tpmz_thymo_meanz$V3 <- "Thymocytes"

# add tissue_id and remove excess columns
thymocytes_across$tissue_id <- "Thymocytes"
thymocytes_across_short <- dplyr::select(thymocytes_across, -contig)


# make data frame, change colnames
thymo_tpmz_meanz <- data.frame(tpmz_thymo_meanz)
colnames(thymo_tpmz_meanz) <- c("gene", "meanz", "tissue_id")

# rbind thymocyte and 10vs10 tissue TPM frames
smash_tpm <- rbind(data.frame(df_TPM_tissues_anno_meanz), thymo_tpmz_meanz)

# change colnames
df_tissues <- df_tissues[,c("target_id", "pval", "qval", "b", "tissue_id")]
colnames(df_tissues) <- c("gene", "pval", "qval", "log2FC", "tissue_id")

# rbind thymocyte and 10vs10 tissue slueth output data frames
df_tissue <- rbind(df_tissues, thymocytes_across_short)

# merge TPM and slueth output data frames
df_tissue_tpmz <- merge(df_tissue, smash_tpm, by = c("gene", "tissue_id"), all.x = T)

# add XCI status annotation
df_tissue_tpmz_anno <- merge(df_tissue_tpmz, chrx_annotation, by.x = "gene", by.y = "Gene_name")

# Rename tissues
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Adipose_Subcutaneous",]$tissue_id <- "adipose-1"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Adipose_Visceral_Omentum",]$tissue_id <- "adipose-2"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Artery_Aorta",]$tissue_id <- "artery-1"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Artery_Coronary",]$tissue_id <- "artery-2"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Artery_Tibial",]$tissue_id <- "artery-3"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Brain_Amygdala",]$tissue_id <- "brain-1"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Brain_Anterior_cingulate_cortex_BA24",]$tissue_id <- "brain-2"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Brain_Caudate_basal_ganglia",]$tissue_id <- "brain-3"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Brain_Cerebellar_Hemisphere",]$tissue_id <- "brain-4"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Brain_Cerebellum",]$tissue_id <- "brain-5"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Brain_Cortex",]$tissue_id <- "brain-6"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Brain_Frontal_Cortex_BA9",]$tissue_id <- "brain-7"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Brain_Hippocampus",]$tissue_id <- "brain-8"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Brain_Hypothalamus",]$tissue_id <- "brain-9"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Brain_Nucleus_accumbens_basal_ganglia",]$tissue_id <- "brain-10"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Brain_Putamen_basal_ganglia",]$tissue_id <- "brain-11"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Brain_Spinal_cord_cervical_c1",]$tissue_id <- "brain-12"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Brain_Substantia_nigra",]$tissue_id <- "brain-13"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Breast_Mammary_Tissue",]$tissue_id <- "breast"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Minor_Salivary_Gland",]$tissue_id <- "salivary gland"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Colon_Sigmoid",]$tissue_id <- "colon-1"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Colon_Transverse",]$tissue_id <- "colon-2"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Esophagus_Gastroesophageal_Junction",]$tissue_id <- "esophagus-1"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Esophagus_Mucosa",]$tissue_id <- "esophagus-2"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Esophagus_Muscularis",]$tissue_id <- "esophagus-3"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Cells_Cultured_fibroblasts",]$tissue_id <- "fibroblasts"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Heart_Atrial_Appendage",]$tissue_id <- "heart-1"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Heart_Left_Ventricle",]$tissue_id <- "heart-2"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Kidney_Cortex",]$tissue_id <- "kidney-1"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Kidney_Medulla",]$tissue_id <- "kidney-2"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Cells_EBV_transformed_lymphocytes",]$tissue_id <- "lymphocytes"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Muscle_Skeletal",]$tissue_id <- "muscle"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Nerve_Tibial",]$tissue_id <- "nerve"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Skin_Not_Sun_Exposed_Suprapubic",]$tissue_id <- "skin-1"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Skin_Sun_Exposed_Lower_leg",]$tissue_id <- "skin-2"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Small_Intestine_Terminal_Ileum",]$tissue_id <- "small intestine"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Adrenal_Gland",]$tissue_id <- "adrenal gland"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Bladder",]$tissue_id <- "bladder"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Liver",]$tissue_id <- "liver"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Lung",]$tissue_id <- "lung"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Pancreas",]$tissue_id <- "pancreas"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Pituitary",]$tissue_id <- "pituitary"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Spleen",]$tissue_id <- "spleen"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Stomach",]$tissue_id <- "stomach"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Thymocytes",]$tissue_id <- "thymocytes"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Thyroid",]$tissue_id <- "thyroid"
df_tissue_tpmz_anno[df_tissue_tpmz_anno$tissue_id == "Whole_Blood",]$tissue_id <- "whole blood"


# split the df into PAR and non-PAR escapees.
df_tissue_PAR_order <- df_tissue_tpmz_anno[df_tissue_tpmz_anno$meanz > 1 & df_tissue_tpmz_anno$gene %in% chrx_annotation[chrx_annotation$category == "PAR",]$Gene_name,] %>% dplyr::group_by(tissue_id) %>% dplyr::summarise(medianz=median(log2FC))
df_tissue_PAR_order <- df_tissue_PAR_order[order(df_tissue_PAR_order$medianz, decreasing = T),]

df_tissue_esc_order <- df_tissue_tpmz_anno[df_tissue_tpmz_anno$meanz > 1 & df_tissue_tpmz_anno$gene %in% chrx_annotation[chrx_annotation$category == "Escape",]$Gene_name,] %>% dplyr::group_by(tissue_id) %>% dplyr::summarise(medianz=median(log2FC))
df_tissue_esc_order <- df_tissue_esc_order[order(df_tissue_esc_order$medianz, decreasing = T),]


# Cap log2FC values for nicer plotting
keklol <- df_tissue_tpmz_anno
keklol[keklol$log2FC > 1,]$log2FC <- 1.02
keklol[keklol$log2FC < -1,]$log2FC <- -1.02

# create plot order separate for each plot df.
keklol_PAR <- keklol
keklol_PAR$tissue_id <- factor(keklol_PAR$tissue_id, levels = df_tissue_PAR_order$tissue_id)

keklol_escape <- keklol
keklol_escape$tissue_id <- factor(keklol_escape$tissue_id, levels = rev(df_tissue_esc_order$tissue_id))

ggsave2("04_plots/Suppl_Figure_3c_3d.pdf", height = 8, width = 14,
        plot_grid(ncol=1,
                  ggplot(keklol_PAR[keklol_PAR$category == "PAR" & keklol_PAR$meanz > 1,], 
                         aes(x=tissue_id, y=log2FC, col=tissue_id))+ 
                    geom_quasirandom(alpha = 0.5, width = 0.3) +
                    geom_boxplot(outlier.shape = NA) + 
                    geom_hline(yintercept = c(0,1,-1),lty=2) + 
                    theme_AL_box_rotX() +
                    labs(x="", y="FC (log2)", title = "PAR")+
                    theme(legend.title = element_blank(), legend.position = "none")+
                    coord_cartesian(ylim=c(-1.02,1.02))+
                    scale_y_continuous(breaks = c(-1,0,1))+
                    scale_color_manual(values = color_valzz)
                  ,
                  
                  ggplot(keklol_escape[keklol_escape$category == "Escape" & keklol_escape$meanz > 1,], 
                         aes(x=tissue_id, y=log2FC, col=tissue_id))+ 
                    geom_quasirandom(alpha = 0.5, width = 0.3) +  
                    geom_boxplot(outlier.shape = NA) + 
                    geom_hline(yintercept = c(0,1,-1),lty=2) + 
                    theme_AL_box_rotX() +
                    labs(x="", y="FC (log2)", title = "Escape")+
                    theme(legend.title = element_blank(), legend.position = "none")+
                    coord_cartesian(ylim=c(-1.02,1.02))+
                    scale_y_continuous(breaks = c(-1,0,1))+
                    scale_color_manual(values = color_valzz)
        )
)

# export source data
sfig_3c_source <- dplyr::select(keklol_PAR[keklol_PAR$category == "PAR" & keklol_PAR$meanz > 1,], gene,   tissue_id  ,             pval      ,         qval     , log2FC, category)
sfig_3d_source <- dplyr::select(keklol_escape[keklol_escape$category == "Escape" & keklol_escape$meanz > 1,], gene,   tissue_id  ,             pval      ,         qval     , log2FC, category)

write.table(sfig_3c_source, file = "06_source_data/Suppl_fig_3c.tsv", quote = F, row.names = F, sep = "\t")
write.table(sfig_3d_source, file = "06_source_data/Suppl_fig_3d.tsv", quote = F, row.names = F, sep = "\t")
