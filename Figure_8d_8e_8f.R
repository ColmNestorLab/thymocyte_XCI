options(stringsAsFactors = F)
library(data.table)
library(immunarch)
library(cowplot)
source("plot_parameters.R")

# Read in TRUST4 data
immdata_sans_merged <- repLoad(.path = "01_processed_data/Bulk_RNAseq_TCR_analysis/for_immunarch/")

# for the downsampling, keep only cells that have actually gone through TCR rearrangment (i.e. only cells after DN).
cells_of_interest <- names(immdata_sans_merged$data)[grepl(names(immdata_sans_merged$data), pattern = paste0(c("DPearly","DPlate","CD4SP","CD8SP"), collapse = "|"))]
down_samp <- repSample(
  immdata_sans_merged$data[c(cells_of_interest)],
  .method = c("downsample"),
  .n = NA
)

ggsave2(filename = "04_plots/Figure_8d.pdf", width = 5.5,
        plot_grid(repExplore(immdata_sans_merged$data, "volume") %>% vis(.test = F, .by = c("Cell", "Sample_id"), .meta = immdata_sans_merged$meta) + theme(text = element_text(size=8), axis.text.x = element_text(size=6))))

ggsave2(filename = "04_plots/Figure_8e.pdf", width = 4,
        plot_grid(repDiversity(down_samp, "d50") %>% vis(.test = F, .by = c("Sample_id"), .meta = immdata_sans_merged$meta) + theme(text = element_text(size=8), axis.text.x = element_text(size=6)) + labs(x="")))

ggsave2(filename = "04_plots/Figure_8f.pdf", width = 8,
        plot_grid(plot_grid(repClonality(down_samp, .method = "rare") %>% vis(.test=T, .by = c("Sample_id"), .meta = immdata_sans_merged$meta) + theme(text = element_text(size=8), axis.text.x = element_text(size=6)))))


#export source data. Excluding actual TCR sequences as they are sensitive data.
fig_8d_source <- map_df(names(immdata_sans_merged$data), function(name) {
  df <- immdata_sans_merged$data[[name]]
  df$Sample <- name
  df$Index <- 1:nrow(df)
  return(df)
})

fig_8d_source <- fig_8d_source[,c("Clones", "Proportion","Index", "Sample"),]
fig_8d_source <- fig_8d_source %>% tidyr::separate(Sample, into = c("sample", "celltype"))
write.table(fig_8d_source, file = "06_source_data/Fig_8d.tsv", quote = F, row.names = F, sep = "\t")


fig_8fe_source <- map_df(names(down_samp), function(name) {
  df <- down_samp[[name]]
  df$Sample <- name
  df$Index <- 1:nrow(df)
  return(df)
})

fig_8fe_source <- fig_8fe_source[,c("Clones", "Proportion","Index", "Sample"),]
fig_8fe_source <- fig_8fe_source %>% tidyr::separate(Sample, into = c("sample", "celltype"))
write.table(fig_8fe_source, file = "06_source_data/Fig_8e.tsv", quote = F, row.names = F, sep = "\t")
write.table(fig_8fe_source, file = "06_source_data/Fig_8f.tsv", quote = F, row.names = F, sep = "\t")
