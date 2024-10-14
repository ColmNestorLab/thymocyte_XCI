library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)

source("plot_parameters.R")

# Read in 10vs10 sleuth table
df_tissues <- fread("02_tidy_data/Suppl_tables/Supplementary_Data_Sleuth_table_tissue_comparison.tsv")

# get MSL3 and KDM6a values, make sure b (log2FC) is numeric
msl3 <- df_tissues[df_tissues$target_id %in% c("MSL3", "KDM6A"),]
msl3$b <- as.numeric(msl3$b)

# rename tissues
msl3[msl3$tissue_id == "Adipose_Subcutaneous",]$tissue_id <- "adipose-1"
msl3[msl3$tissue_id == "Adipose_Visceral_Omentum",]$tissue_id <- "adipose-2"
msl3[msl3$tissue_id == "Artery_Aorta",]$tissue_id <- "artery-1"
msl3[msl3$tissue_id == "Artery_Coronary",]$tissue_id <- "artery-2"
msl3[msl3$tissue_id == "Artery_Tibial",]$tissue_id <- "artery-3"
msl3[msl3$tissue_id == "Brain_Amygdala",]$tissue_id <- "brain-1"
msl3[msl3$tissue_id == "Brain_Anterior_cingulate_cortex_BA24",]$tissue_id <- "brain-2"
msl3[msl3$tissue_id == "Brain_Caudate_basal_ganglia",]$tissue_id <- "brain-3"
msl3[msl3$tissue_id == "Brain_Cerebellar_Hemisphere",]$tissue_id <- "brain-4"
msl3[msl3$tissue_id == "Brain_Cerebellum",]$tissue_id <- "brain-5"
msl3[msl3$tissue_id == "Brain_Cortex",]$tissue_id <- "brain-6"
msl3[msl3$tissue_id == "Brain_Frontal_Cortex_BA9",]$tissue_id <- "brain-7"
msl3[msl3$tissue_id == "Brain_Hippocampus",]$tissue_id <- "brain-8"
msl3[msl3$tissue_id == "Brain_Hypothalamus",]$tissue_id <- "brain-9"
msl3[msl3$tissue_id == "Brain_Nucleus_accumbens_basal_ganglia",]$tissue_id <- "brain-10"
msl3[msl3$tissue_id == "Brain_Putamen_basal_ganglia",]$tissue_id <- "brain-11"
msl3[msl3$tissue_id == "Brain_Spinal_cord_cervical_c1",]$tissue_id <- "brain-12"
msl3[msl3$tissue_id == "Brain_Substantia_nigra",]$tissue_id <- "brain-13"
msl3[msl3$tissue_id == "Breast_Mammary_Tissue",]$tissue_id <- "breast"
msl3[msl3$tissue_id == "Minor_Salivary_Gland",]$tissue_id <- "salivary gland"
msl3[msl3$tissue_id == "Colon_Sigmoid",]$tissue_id <- "colon-1"
msl3[msl3$tissue_id == "Colon_Transverse",]$tissue_id <- "colon-2"
msl3[msl3$tissue_id == "Esophagus_Gastroesophageal_Junction",]$tissue_id <- "esophagus-1"
msl3[msl3$tissue_id == "Esophagus_Mucosa",]$tissue_id <- "esophagus-2"
msl3[msl3$tissue_id == "Esophagus_Muscularis",]$tissue_id <- "esophagus-3"
msl3[msl3$tissue_id == "Cells_Cultured_fibroblasts",]$tissue_id <- "fibroblasts"
msl3[msl3$tissue_id == "Heart_Atrial_Appendage",]$tissue_id <- "heart-1"
msl3[msl3$tissue_id == "Heart_Left_Ventricle",]$tissue_id <- "heart-2"
msl3[msl3$tissue_id == "Kidney_Cortex",]$tissue_id <- "kidney-1"
msl3[msl3$tissue_id == "Kidney_Medulla",]$tissue_id <- "kidney-2"
msl3[msl3$tissue_id == "Cells_EBV_transformed_lymphocytes",]$tissue_id <- "lymphocytes"
msl3[msl3$tissue_id == "Muscle_Skeletal",]$tissue_id <- "muscle"
msl3[msl3$tissue_id == "Nerve_Tibial",]$tissue_id <- "nerve"
msl3[msl3$tissue_id == "Skin_Not_Sun_Exposed_Suprapubic",]$tissue_id <- "skin-1"
msl3[msl3$tissue_id == "Skin_Sun_Exposed_Lower_leg",]$tissue_id <- "skin-2"
msl3[msl3$tissue_id == "Small_Intestine_Terminal_Ileum",]$tissue_id <- "small intestine"
msl3[msl3$tissue_id == "Adrenal_Gland",]$tissue_id <- "adrenal gland"
msl3[msl3$tissue_id == "Bladder",]$tissue_id <- "bladder"
msl3[msl3$tissue_id == "Liver",]$tissue_id <- "liver"
msl3[msl3$tissue_id == "Lung",]$tissue_id <- "lung"
msl3[msl3$tissue_id == "Pancreas",]$tissue_id <- "pancreas"
msl3[msl3$tissue_id == "Pituitary",]$tissue_id <- "pituitary"
msl3[msl3$tissue_id == "Spleen",]$tissue_id <- "spleen"
msl3[msl3$tissue_id == "Stomach",]$tissue_id <- "stomach"
msl3[msl3$tissue_id == "Thymocytes",]$tissue_id <- "thymocytes"
msl3[msl3$tissue_id == "Thyroid",]$tissue_id <- "thyroid"
msl3[msl3$tissue_id == "Whole_Blood",]$tissue_id <- "whole blood"

# get plot order, # level tissue_id accordingly
orderz <- msl3[msl3$target_id == "MSL3",]
orderz <- orderz[order(orderz$b, decreasing = T)]
msl3$tissue_id <- factor(msl3$tissue_id, levels = orderz$tissue_id)

#plot
ggsave2(filename = "04_plots/Supplementary_fig_7f.pdf", width = 8, height = 4,
        plot_grid(
ggplot(msl3, aes(x=tissue_id,y=b, col = target_id)) + 
  geom_point() + 
  theme_AL_box_rotX() + 
  geom_hline(yintercept = 0, lty=2)+
  labs(x="", y="FC(log2)")+
  theme(legend.title = element_blank())+
  coord_cartesian(ylim=c(-0.26, 0.76))+
  scale_color_npg()
  )
)

# Export source data
write.table(msl3, file = "06_source_data/Suppl_Fig_7f.tsv", quote = F, row.names = F, sep = "\t")

