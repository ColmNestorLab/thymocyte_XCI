library(data.table)
library(ggplot2)
library(cowplot)
library(readxl)
library(tidyverse)
source("plot_parameters.R")

# read in thymocyte TPM matrix
df_tpmz <- data.table(read_excel("02_tidy_data/Suppl_tables/Supplementary Tables.xlsx", sheet = 3, skip = 1))

# select genes to plot
goiz <- c("MYCN", "LMO2", "LYL1", "RAG2", "RORC", "TCF12", "ZBTB7B", "GATA3", "RUNX3")
df_tpmz_filt <- melt(df_tpmz[df_tpmz$gene %in% goiz,], id.vars = "gene")

# separate variable into sample and celltype
df_tpmz_filt <- df_tpmz_filt %>% tidyr::separate(variable, into = c("sample", "celltype", NA, NA), sep = "_")

# plot
ggsave2("04_plots/figure_1c.pdf", width = 12, height = 6,
ggplot(df_tpmz_filt, aes(x=factor(celltype,levels=cell_order), y=as.numeric(value))) + geom_boxplot(outlier.shape = NA) + geom_point(aes(col = celltype)) + facet_wrap(~factor(gene, levels = goiz), scales = "free_y", nrow=1) + theme_AL_box_rotX() + coord_cartesian(ylim=c(0, NA))
)

# export source data
colnames(df_tpmz_filt) <- c("gene", "sample", "celltype", "TPM")
write.table(df_tpmz_filt, file = "06_source_data/Fig_1c.tsv", quote = F, row.names = F, sep = "\t")