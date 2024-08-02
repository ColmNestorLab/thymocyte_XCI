library(data.table)
library(ggplot2)
library(cowplot)
source("03_r_scripts/plot_parameters.R")

df_tpmz <- fread("Supplementary Table 3")

goiz <- c("MYCN", "LMO2", "LYL1", "RAG2", "RORC", "TCF12", "ZBTB7B", "GATA3", "RUNX3")

df_tpmz_filt <- melt(df_tpmz[df_tpmz$gene %in% goiz,])

df_tpmz_filt <- df_tpmz_filt %>% tidyr::separate(variable, into = c("sample", "celltype", NA, NA), sep = "_")

ggsave2("04_plots/figure_1C.pdf", width = 12, height = 6,
ggplot(df_tpmz_filt, aes(x=factor(celltype,levels=cell_order), y=value)) + geom_boxplot(outlier.shape = NA) + geom_point(aes(col = celltype)) + facet_wrap(~factor(gene, levels = goiz), scales = "free_y", nrow=1) + theme_AL_box_rotX() + coord_cartesian(ylim=c(0, NA))
)
