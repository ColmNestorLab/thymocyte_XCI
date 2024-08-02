library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggbeeswarm)
source("03_r_scripts/plot_parameters.R")


#read in annotations (PAR and the Tukiainen annotation)
PAR <- rbind(fread("97_indexes_annotations/group-715_PAR1.csv"), fread("97_indexes_annotations/group-716_PAR2.csv"))
esc <- fread("97_indexes_annotations/landscape.Suppl.Table.13.csv")
esc[esc$Gene_name == "6-Sep",]$Gene_name <- "SEPT6"
esc <- esc[,-1]

# merge annotations
chrx_annotation <- merge(esc, PAR, by.x = "Gene_name" , by.y = "Approved symbol" , all = T)[-1,]
chrx_annotation$category <- chrx_annotation$Reported_XCI_status
chrx_annotation[chrx_annotation$PAR == "PAR1" | chrx_annotation$PAR == "PAR2",]$category <- "PAR"



meta_full <- fread("96_metadata/metadata_RNAseq_full_for_tissues.tsv")
df_TPM_tissues <- dplyr::select(melt(fread("02_tidy_data/GTEX/tissue_comparison/10vs10_TPM_matrix.tsv")), -V1)

meta_tissues <- subset(meta_full, meta_full$`entity:sample_id` %in% gsub(unique(df_TPM_tissues$variable), pattern = "\\.", replacement = "-"))

meta_tissues[meta_tissues$tissue_id %in% c("Cells_EBV-transformed_lymphocytes"),]$tissue_id <- "Cells_EBV_transformed_lymphocytes"

meta_tissues$sample <- gsub(meta_tissues$`entity:sample_id`, pattern = "-", replacement = "\\.")


df_tissues <- fread("Supplementary Table 17")

asdfffffffff <- df_tissues[!df_tissues$target_id %in% c("XIST"),] %>% dplyr::group_by(tissue_id) %>% dplyr::summarise(meanzz=mean(b))
asdfffffffff <- asdfffffffff[order(asdfffffffff$meanzz, decreasing = T),]


## counts ##
meta_tissues_count <- meta_tissues %>% dplyr::group_by(tissue_id, sex) %>% dplyr::count()

ggplot(meta_tissues_count, aes(x=tissue_id, y=n, fill=sex)) + geom_bar(stat="identity") + geom_text(aes(label=n)) + theme_AL_box_rotX()


unique(meta_tissues_count$tissue_id)
length(meta_tissues$`entity:sample_id`)


# remove NAs and headers left while merging
df_tissues <- df_tissues[!is.na(df_tissues$b),]
df_tissues <- subset(df_tissues, target_id != "target_id")
df_tissues$b <- as.numeric(df_tissues$b)
df_tissues <- df_tissues[df_tissues$tissue_id != "All_tissues",]


df_TPM_tissues_anno <- merge(df_TPM_tissues, meta_tissues, by.x = "variable", by.y = "sample")

df_TPM_tissues_anno[df_TPM_tissues_anno$tissue_id == "Brain_Spinal_cord_cervical_c-1",]$tissue_id <- "Brain_Spinal_cord_cervical_c1"

test <- df_TPM_tissues_anno
test$category <- "toss"
test[test$gene %in% chrx_annotation[chrx_annotation$category == "PAR",]$Gene_name,]$category <- "PAR"
test[test$gene %in% "XIST",]$category <- "XIST"

category_meanz <- test %>% dplyr::group_by(tissue_id, category, sex) %>% dplyr::summarise(meanz = mean(value))

test_order <- category_meanz[category_meanz$category == "PAR" & category_meanz$sex == "Female",]
test_order <- test_order[order(test_order$meanz, decreasing = T),]

ggplot(category_meanz[category_meanz$category != "toss",], aes(x=tissue_id, y = meanz, col=sex)) + geom_line(aes(group=sex)) + facet_wrap(~category, scales = "free_y")

ggplot(df_TPM_tissues_anno[df_TPM_tissues_anno$gene == "XIST",], aes(x=tissue_id, y=value, col=sex)) + geom_boxplot() + theme_AL_box_rotX()

df_TPM_tissues_anno_meanz <- df_TPM_tissues_anno %>% dplyr::group_by(tissue_id, gene) %>% dplyr::summarise(meanz = mean(value))



thymocytes_across <- fread("Supplementary Table 5")

tpmz_thymo <- melt(fread("Supplementary Table 3"))

tpmz_thymo_meanz <- tpmz_thymo %>% dplyr::group_by(gene) %>% dplyr::summarise(meanz = mean(value))
tpmz_thymo_meanz$V3 <- "Thymocytes"

thymocytes_across$tissue_id <- "Thymocytes"
thymocytes_across_short <- dplyr::select(thymocytes_across, -seqnames)



thymo_tpmz_meanz <- data.frame(tpmz_thymo_meanz)

colnames(thymo_tpmz_meanz) <- c("gene", "meanz", "tissue_id")

smash_tpm <- rbind(data.frame(df_TPM_tissues_anno_meanz), thymo_tpmz_meanz)




df_tissue <- rbind(df_tissues, thymocytes_across_short)

df_tissue_tpmz <- merge(df_tissue, smash_tpm, by.x = c("target_id", "tissue_id"), by.y = c("gene", "tissue_id"), all.x = T)

df_tissue_tpmz_anno <- merge(df_tissue_tpmz, chrx_annotation, by.x = "target_id", by.y = "Gene_name")


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



df_tissue_PAR_order <- df_tissue_tpmz_anno[df_tissue_tpmz_anno$meanz > 1 & df_tissue_tpmz_anno$target_id %in% chrx_annotation[chrx_annotation$category == "PAR",]$Gene_name,] %>% dplyr::group_by(tissue_id) %>% dplyr::summarise(medianz=median(b))

df_tissue_PAR_order <- df_tissue_PAR_order[order(df_tissue_PAR_order$medianz, decreasing = T),]


df_tissue_esc_order <- df_tissue_tpmz_anno[df_tissue_tpmz_anno$meanz > 1 & df_tissue_tpmz_anno$target_id %in% chrx_annotation[chrx_annotation$category == "Escape",]$Gene_name,] %>% dplyr::group_by(tissue_id) %>% dplyr::summarise(medianz=median(b))

df_tissue_esc_order <- df_tissue_esc_order[order(df_tissue_esc_order$medianz, decreasing = T),]




keklol <- df_tissue_tpmz_anno
keklol[keklol$b > 1,]$b <- 1.02
keklol[keklol$b < -1,]$b <- -1.02


colorado <- function(src, boulder) {
  if (!is.factor(src)) src <- factor(src)                   # make sure it's a factor
  src_levels <- levels(src)                                 # retrieve the levels in their order
  brave <- boulder %in% src_levels                          # make sure everything we want to make bold is actually in the factor levels
  if (all(brave)) {                                         # if so
    b_pos <- purrr::map_int(boulder, ~which(.==src_levels)) # then find out where they are
    b_vec <- rep("plain", length(src_levels))               # make'm all plain first
    b_vec[b_pos] <- "bold"                                  # make our targets bold
    b_vec                                                   # return the new vector
  } else {
    stop("All elements of 'boulder' must be in src")
  }
}

keklol_PAR <- keklol
keklol_PAR$tissue_id <- factor(keklol_PAR$tissue_id, levels = df_tissue_PAR_order$tissue_id)

keklol_escape <- keklol
keklol_escape$tissue_id <- factor(keklol_escape$tissue_id, levels = rev(df_tissue_esc_order$tissue_id))

length(unique(keklol_escape$tissue_id))

library(ggsci)
kek <- c(pal_npg("nrc")(10), pal_nejm(palette = "default")(8), pal_jama(palette = "default")(7), pal_jco(palette = "default")(7))

ggsave2("04_plots/Suppl_Figure_7B_7C.pdf", height = 8, width = 14,
        plot_grid(ncol=1,
                  ggplot(keklol_PAR[keklol_PAR$category == "PAR" & keklol_PAR$meanz > 1,], 
                         aes(x=tissue_id, y=b, col=tissue_id))+ 
                    geom_quasirandom(alpha = 0.5, width = 0.3) +
                    geom_boxplot(outlier.shape = NA) + 
                    geom_hline(yintercept = c(0,1,-1),lty=2) + 
                    theme_AL_box_rotX() +
                    labs(x="", y="FC (log2)", title = "PAR")+
                    theme(legend.title = element_blank(), legend.position = "none")+
                    coord_cartesian(ylim=c(-1.02,1.02))+
                    scale_y_continuous(breaks = c(-1,0,1))+
                    theme(axis.text.x=element_text(face=colorado(keklol_PAR$tissue_id, c("thymocytes"))))+
                    scale_color_manual(values = color_valzz)
                  ,
                  
                  ggplot(keklol_escape[keklol_escape$category == "Escape" & keklol_escape$meanz > 1,], 
                         aes(x=tissue_id, y=b, col=tissue_id))+ 
                    geom_quasirandom(alpha = 0.5, width = 0.3) +  
                    geom_boxplot(outlier.shape = NA) + 
                    geom_hline(yintercept = c(0,1,-1),lty=2) + 
                    theme_AL_box_rotX() +
                    labs(x="", y="FC (log2)", title = "Escape")+
                    theme(legend.title = element_blank(), legend.position = "none")+
                    coord_cartesian(ylim=c(-1.02,1.02))+
                    scale_y_continuous(breaks = c(-1,0,1))+
                    theme(axis.text.x=element_text(face=colorado(keklol_escape$tissue_id, c("thymocytes"))))+
                    scale_color_manual(values = color_valzz)
        )
)


