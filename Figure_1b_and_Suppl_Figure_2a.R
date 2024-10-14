#Load sleuth object
source("Load_thymocyte_sleuth_expression_data.R")

library(cowplot)

## PLOT PCA
library(ggfortify)

# function to plot which genes contribute to the principal components.
plot.prcomp.loadings <- function(x,pc,n){
  require(ggplot2)
  lds <- sort(x$rotation[,pc],decreasing = T)
  lds_filt <- c(head(lds,n),tail(lds,n))
  qplot(x=reorder(names(lds_filt),-lds_filt,mean),y=lds_filt,geom="col") + labs(y="Contribution", title=paste("Principle component",pc)) + theme(axis.title.x = element_blank())
}

# run prcomp
pca_thymo <- prcomp(t(txi_thymo_transition$counts[,1:40]))

# plot and store
PCA_plot <- autoplot(pca_thymo,data=meta_thymo_transition[1:40,],colour="Celltype",shape = "Sex" ,size=3,stroke=NA, x = 1, y=2) + theme_AL_box()

# Export PCA 1 and PCA 2
lel <- data.frame(pca_thymo$rotation[,1:2])
lel$gene <- rownames(lel)

ggsave2(filename = "04_plots/figure_1B.pdf", width = 150, 
        height = 100, units = "mm", 
        plot_grid(PCA_plot)
)

# Export
ggsave("04_plots/Suppl_Figure_2a.pdf",width = 15,height = 8,
       gridExtra::grid.arrange(
         plot.prcomp.loadings(pca_thymo,1,10) + theme_AL_box_rotX(axis.title.x = element_blank()) + geom_hline(yintercept = 0, lty=2),
         plot.prcomp.loadings(pca_thymo,2,10) + theme_AL_box_rotX(axis.title.x = element_blank()) + geom_hline(yintercept = 0, lty=2)
       )
)

# Export source data
PC1_contr <- sort(pca_thymo$rotation[,1], decreasing = T)
PC1_contr_filt <- data.frame(c(head(PC1_contr,10),tail(PC1_contr,10)))
PC2_contr <- sort(pca_thymo$rotation[,2], decreasing = T)
PC2_contr_filt <- data.frame(c(head(PC2_contr,10),tail(PC2_contr,10)))

colnames(PC1_contr_filt) <- "Contribution"
PC1_contr_filt$gene <- rownames(PC1_contr_filt)
PC1_contr_filt$PC <- 1

colnames(PC2_contr_filt) <- "Contribution"
PC2_contr_filt$gene <- rownames(PC2_contr_filt)
PC2_contr_filt$PC <- 2

write.table(data.frame(rbind(PC1_contr_filt, PC2_contr_filt)), file = "06_source_data/Suppl_fig_2a.tsv", quote = F, row.names = F, sep = "\t")

write.table(PCA_plot$data[,c("PC1", "PC2", "Sample", "Patient_ID", "Sex", "Celltype")], file = "06_source_data/fig_1b.tsv", quote = F, row.names = F, sep = "\t")
