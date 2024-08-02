#Load sleuth object
source("03_r_scripts/load_scripts/Load_thymocyte_sleuth_expression_data.R")

library(cowplot)

##PLOT PCA
library(ggfortify)
plot.prcomp.loadings <- function(x,pc,n){
  require(ggplot2)
  lds <- sort(x$rotation[,pc],decreasing = T)
  lds_filt <- c(head(lds,n),tail(lds,n))
  qplot(x=reorder(names(lds_filt),-lds_filt,mean),y=lds_filt,geom="col") + labs(y="Contribution", title=paste("Principle component",pc)) + theme(axis.title.x = element_blank())
}

pca_thymo <- prcomp(t(txi_thymo_transition$counts[,1:40]))

meta_thymo_transition$Celltype <- factor(meta_thymo_transition$Celltype, levels = c("ETP", "DN", "DPearly", "DPlate", "CD4SP", "CD8SP"))
PCA_plot <- autoplot(pca_thymo,data=meta_thymo_transition[1:40,],colour="Celltype",shape = "Sex" ,size=3,stroke=NA, x = 1, y=2) + theme_AL_box()

# Export PCA 1 and PCA 2
lel <- data.frame(pca_thymo$rotation[,1:2])
lel$gene <- rownames(lel)

ggsave2(filename = "04_plots/figure_1B.pdf", width = 150, 
        height = 100, units = "mm", 
        plot_grid(PCA_plot)
)

# Export as A4.
ggsave("04_plots/Suppl_Figure_1B.pdf",width = 15,height = 8,
       gridExtra::grid.arrange(
         plot.prcomp.loadings(pca_thymo,1,10) + theme_AL_box_rotX(axis.title.x = element_blank()) + geom_hline(yintercept = 0, lty=2),
         plot.prcomp.loadings(pca_thymo,2,10) + theme_AL_box_rotX(axis.title.x = element_blank()) + geom_hline(yintercept = 0, lty=2)
       )
)