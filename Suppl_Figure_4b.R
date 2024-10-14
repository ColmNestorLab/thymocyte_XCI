library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(ggrastr)
library(readxl)

source("plot_parameters.R")

point_size <- 1.5

####### Read in data ########
ase_unfilt <- fread("02_tidy_data/Suppl_tables/Supplementary_data_ASE_table_all_chroms_and_samples.tsv")

ase.f3_unfilt_with_chrx_info <- ase_unfilt[!is.na(ase_unfilt$gene) & sample == "F3" & contig == "chrX"]

###### add PAR info #####
PAR <- rbind(fread("97_indexes_annotations/group-715_PAR1.csv"), fread("97_indexes_annotations/group-716_PAR2.csv"))
PAR$GENES <- PAR$`Approved symbol`
ase.f3_unfilt_with_chrx_info <- merge(ase.f3_unfilt_with_chrx_info, PAR, by.x = "gene", by.y = "GENES", all.x = T)
ase.f3_unfilt_with_chrx_info[is.na(PAR)]$PAR <- "NON_PAR"
ase.f3_unfilt_with_chrx_info[ase.f3_unfilt_with_chrx_info$PAR == "PAR1" | ase.f3_unfilt_with_chrx_info$PAR == "PAR2"]$PAR <- "PAR"

##### add Esc info ######
x.esc <- dplyr::select(fread("97_indexes_annotations/landscape.Suppl.Table.13.csv"), Gene_name, Reported_XCI_status, Sex_bias_in_GTEx, XCI_across_tissues, XCI_in_single_cells)
x.esc[x.esc$Gene_name == "6-Sep"]$Gene_name <- "SEPT6"
x.esc <- subset(x.esc, !(Reported_XCI_status == "Unknown" & Gene_name == "IDS"))
x.esc <- x.esc[x.esc$Gene_name != ""]

##### add info to df #########
ase.f3_unfilt_with_chrx_info <- merge(ase.f3_unfilt_with_chrx_info, x.esc, by.x = "gene", by.y = "Gene_name", all.x = T)

# remove overlapping hSNPs
ase_unfilt <- ase_unfilt[!grepl(ase_unfilt$gene, pattern = ";")]
ase.f3_unfilt_with_chrx_info <- ase.f3_unfilt_with_chrx_info[!grepl(ase.f3_unfilt_with_chrx_info$gene, pattern = ";")]

# change cell order
ase_unfilt$cell <- factor(ase_unfilt$cell, levels = cell_order)
ase.f3_unfilt_with_chrx_info$cell <- factor(ase.f3_unfilt_with_chrx_info$cell, levels = cell_order)

################# Small fixes ####################
ase.f3_unfilt_with_chrx_info[ase.f3_unfilt_with_chrx_info$Reported_XCI_status == "Unknown"]$Reported_XCI_status <- "Potential"
ase.f3_unfilt_with_chrx_info[is.na(ase.f3_unfilt_with_chrx_info$Reported_XCI_status)]$Reported_XCI_status <- "Potential"
ase.f3_unfilt_with_chrx_info$status <- ifelse(ase.f3_unfilt_with_chrx_info$effectSize > 0.4, yes = "Monoallelic", no = "Biallelic")

ase_unfilt$cell <- factor(ase_unfilt$cell, levels = cell_order)
ase.f3_unfilt_with_chrx_info$cell <- factor(ase.f3_unfilt_with_chrx_info$cell, levels = cell_order)

chrx_gene_order <- fread("97_indexes_annotations/chrx_gene_order.tsv")

## Idenfity 'escape' genes that does not escape in thymocytes
id_genez <- ase.f3_unfilt_with_chrx_info[ase.f3_unfilt_with_chrx_info$keep == T & ase.f3_unfilt_with_chrx_info$Reported_XCI_status == "Escape",]

# add chrX tag
ase_unfilt$chrx <- ifelse(ase_unfilt$contig == "chrX", yes = "chrX", no = "autosome")

######### ASE for F3 ############
p.allelic.genomic.f3 <- 
  ggplot(ase_unfilt[keep == T & sample == "F3"], aes(y= log10(refCount + 1), x= log10(altCount + 1), col=chrx)) +
  geom_abline(slope=1, intercept=0) +
  geom_point_rast(size = point_size) +
  facet_wrap(~cell, ncol = 3, scale = "free") +
  labs(y="log10(ref+1)", x="log10(alt+1)") +
  scale_color_brewer(palette="Set1", name=NULL, labels=c("Autosome", "ChrX")) +
  theme_AL_simple() +
  theme(legend.position = "none", strip.background = element_blank(), axis.title.x = element_text(size=8), axis.title.y = element_text(size=8), strip.text.x = element_text(size = 8), axis.text = element_text(size=8))+
  scale_y_continuous(limits = c(0, 4), expand = c(0.05, 0), breaks = c(0,2,4))+
  scale_x_continuous(limits = c(0, 4), expand = c(0.05, 0), breaks = c(0,2,4))+
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"), strip.text.x = element_text(margin = margin(.01, 0, .01, 0, "cm")))


ggsave2("04_plots/Suppl_Figure_4b.pdf", width = 8,
        plot_grid(p.allelic.genomic.f3))

# Export source data
write.table(ase_unfilt[keep == T & sample == "F3"], file = "06_source_data/Suppl_Fig_4b.tsv", quote = F, row.names = F, sep = "\t")
