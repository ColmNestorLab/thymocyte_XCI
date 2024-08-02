options(stringsAsFactors = F)

######## read in data #############
message("loading sleuth object")
source("03_r_scripts/load_scripts/Load_thymocyte_sleuth_expression_data_with_Turner.R")

# Filter out genes not expressed in any subtype
tpm_filter <- txi_thymo$abundance[apply(txi_thymo$abundance,1,function(x) any(tapply(x,meta_thymo_transition$Celltype,function(y) mean(y>=1) )==1) ),]

##PLOT PCA
library(ggfortify)
pca_thymo <- prcomp(t(txi_thymo$counts))

sample_order <- colnames(txi_thymo$counts)
meta_thymo_transition$Sample <- factor(meta_thymo_transition$Sample, levels = sample_order)

meta_thymo_transition <- meta_thymo_transition[order(meta_thymo_transition$Sample),]
meta_thymo_transition$Sex2 <- ifelse(meta_thymo_transition$Sex == "T", yes = "Turner", no = "Other")

meta_thymo_transition$Celltype <- factor(meta_thymo_transition$Celltype, levels = c("ETP", "DN", "DPearly", "DPlate", "CD4SP", "CD8SP"))

##PLOT PCA
library(ggfortify)

PCA_plot <- autoplot(pca_thymo, data=meta_thymo_transition, col = "Celltype", shape = "Sex", size=3, stroke=NA, x = 1, y=2) + 
  theme_AL_box() + 
  annotate("text", x=-0.1, y=0.2, label= "ETP/DN")  +  
  annotate("text", x=0.2, y=0.066, label= "DPearly/DPlate")  +  
  annotate("text", x=-0.1, y=-0.1, label= "CD4SP/CD8SP")


library(cowplot)
ggsave2("04_plots/Figure_4B.pdf",width = 10,height = 8,
       plot_grid(PCA_plot))
