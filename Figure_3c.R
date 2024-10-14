
library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(ggrastr)

source("plot_parameters.R")

point_size <- 1.5

######## Load data ###########

ase_unfilt <- fread("02_tidy_data/Suppl_tables/Supplementary_data_ASE_table_all_chroms_and_samples.tsv")

# add chrX and autosomal tag
ase_unfilt$chromosome <- "LOL"
ase_unfilt[ase_unfilt$contig == "chrX",]$chromosome <- "chrX"
ase_unfilt[ase_unfilt$contig != "chrX",]$chromosome <- "Autosome"

# calculate means
ase_unfilt_means <- ase_unfilt[keep == T & sample != "M1" & sample != "M2" & sample != "M3"] %>% dplyr::group_by(gene, sample, chromosome) %>% dplyr::summarise(alt_mean = mean(altCount), ref_mean = mean(refCount))

# make sure tag is factor for plotting order.
ase_unfilt_means$chromosome <- factor(ase_unfilt_means$chromosome, levels = c("Autosome","chrX"))

# set genes to label
genes_to_label <- c("XIST", "TSIX", "CD99", "NOTCH1")

# Plot
p.allelic.genomic.all_means <-
  ase_unfilt_means %>% arrange(chromosome) %>%
  ggplot(aes(y= log10(ref_mean + 1), x= log10(alt_mean + 1), col=chromosome)) +
  geom_abline(slope=1, intercept=0) +
  geom_point_rast(data=ase_unfilt_means[!ase_unfilt_means$gene %in% genes_to_label,],size = point_size) +
  geom_point_rast(data=ase_unfilt_means[!ase_unfilt_means$gene %in% genes_to_label,],size = point_size) +
  geom_point_rast(data=ase_unfilt_means[!ase_unfilt_means$gene %in% genes_to_label & ase_unfilt_means$chromosome == "Autosome",],size = point_size, color = "grey")+
  geom_point_rast(data=ase_unfilt_means[!ase_unfilt_means$gene %in% genes_to_label & ase_unfilt_means$chromosome == "chrX",],size = point_size, color = "red")+
  geom_point_rast(data=ase_unfilt_means[ase_unfilt_means$gene %in% genes_to_label,],size = point_size, color = "black") +
  facet_wrap(sample~., scale = "free", ncol = 1) +
  labs(y="log10(ref+1)", x="log10(alt+1)") +
  theme_AL_simple() +
  theme( legend.position = "none", strip.background = element_blank(), axis.title.x = element_text(size=8), axis.title.y = element_text(size=8), strip.text.x = element_text(size = 8), axis.text = element_text(size=8))+
  scale_y_continuous(limits = c(0, 4), expand = c(0.05, 0), breaks = c(0,2,4))+
  scale_x_continuous(limits = c(0, 4), expand = c(0.05, 0), breaks = c(0,2,4))+
  geom_text_repel(data=ase_unfilt_means[ase_unfilt_means$gene %in% genes_to_label,], aes(label=gene), col = "black", nudge_x = 0.75, min.segment.length = 0.00000000001,  size=3)+
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"), strip.text.x = element_text(margin = margin(.01, 0, .01, 0, "cm")))


ggsave2("04_plots/Figure_3c.pdf",
        plot_grid(ncol=3, p.allelic.genomic.all_means, NULL, NULL))


# Export source data
write.table(ase_unfilt_means, file = "06_source_data/Fig_3c.tsv", quote = F, row.names = F, sep = "\t")

