options(stringsAsFactors = F)

######## read in data #############
message("loading sleuth object")
source("Load_thymocyte_sleuth_expression_data_with_Turner.R")
library(ggfortify)
library(cowplot)

# Filter out genes not expressed in any subtype
tpm_filter <- txi_thymo$abundance[apply(txi_thymo$abundance,1,function(x) any(tapply(x,meta_thymo_transition$Celltype,function(y) mean(y>=1) )==1) ),]

# run prcomp
pca_thymo <- prcomp(t(txi_thymo$counts))

# get sample order
sample_order <- colnames(txi_thymo$counts)
# make faactor
meta_thymo_transition$Sample <- factor(meta_thymo_transition$Sample, levels = sample_order)

# order meta, make sure celltype is ordered factor.
meta_thymo_transition <- meta_thymo_transition[order(meta_thymo_transition$Sample),]
meta_thymo_transition$Celltype <- factor(meta_thymo_transition$Celltype, levels = c("ETP", "DN", "DPearly", "DPlate", "CD4SP", "CD8SP"))

## plot
PCA_plot <- autoplot(pca_thymo, data=meta_thymo_transition, col = "Celltype", shape = "Sex", size=3, stroke=NA, x = 1, y=2) + 
  theme_AL_box() + 
  annotate("text", x=-0.1, y=0.2, label= "ETP/DN")  +  
  annotate("text", x=0.2, y=0.066, label= "DPearly/DPlate")  +  
  annotate("text", x=-0.1, y=-0.1, label= "CD4SP/CD8SP")

ggsave2("04_plots/Figure_8b.pdf",width = 10,height = 8,
       plot_grid(PCA_plot))


# Export source data
write.table(PCA_plot$data[,c("PC1", "PC2", "Patient_ID", "Sex", "Celltype")], file = "06_source_data/fig_8b.tsv", quote = F, row.names = F, sep = "\t")
