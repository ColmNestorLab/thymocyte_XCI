options(stringsAsFactors = F)
library(data.table)
library(immunarch)
library(cowplot)
source("03_r_scripts/plot_parameters.R")

# Read in TRUST4 data
immdata_sans_merged <- repLoad(.path = "01_processed_data/Bulk_RNAseq_TCR_analysis/for_immunarch/")

# for the downsampling, keep only cells that have actually gone through TCR rearrangment (i.e. only cells after DN).
cells_of_interest <- names(immdata_sans_merged$data)[grepl(names(immdata_sans_merged$data), pattern = paste0(c("DPearly","DPlate","CD4SP","CD8SP"), collapse = "|"))]

down_samp <- repSample(
  immdata_sans_merged$data[c(cells_of_interest)],
  .method = c("downsample"),
  .n = NA
)

ggsave2(filename = "04_plots/Figure_4D.pdf", width = 5.5,
        plot_grid(repExplore(immdata_sans_merged$data, "volume") %>% vis(.test = F, .by = c("Cell", "Sample_id"), .meta = immdata_sans_merged$meta) + theme(text = element_text(size=8), axis.text.x = element_text(size=6))))

ggsave2(filename = "04_plots/Figure_4E.pdf", width = 4,
        plot_grid(repDiversity(down_samp, "d50") %>% vis(.test = F, .by = c("Sample_id"), .meta = immdata_sans_merged$meta) + theme(text = element_text(size=8), axis.text.x = element_text(size=6)) + labs(x="")))

ggsave2(filename = "04_plots/Figure_4F_sample.pdf", width = 8,
        plot_grid(plot_grid(repClonality(down_samp, .method = "rare") %>% vis(.test=T, .by = c("Sample_id"), .meta = immdata_sans_merged$meta) + theme(text = element_text(size=8), axis.text.x = element_text(size=6)))))
