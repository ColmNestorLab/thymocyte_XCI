######### Figure 2C #############
library(dplyr)
library(data.table)

source("03_r_scripts/plot_parameters.R")

####### Read in data ########
ase.f3_unfilt <- fread("Supplementary Table 8")

ase.f3_unfilt_with_chrx_info <- ase.f3_unfilt[!is.na(ase.f3_unfilt$gene) & sample == "F3" & contig == "chrX"]

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
ase.f3_unfilt <- ase.f3_unfilt[!grepl(ase.f3_unfilt$gene, pattern = ";")]
ase.f3_unfilt_with_chrx_info <- ase.f3_unfilt_with_chrx_info[!grepl(ase.f3_unfilt_with_chrx_info$gene, pattern = ";")]

# change cell order
ase.f3_unfilt$cell <- factor(ase.f3_unfilt$cell, levels = cell_order)
ase.f3_unfilt_with_chrx_info$cell <- factor(ase.f3_unfilt_with_chrx_info$cell, levels = cell_order)

################# Small fixes ####################
ase.f3_unfilt_with_chrx_info[ase.f3_unfilt_with_chrx_info$Reported_XCI_status == "Unknown"]$Reported_XCI_status <- "Potential"
ase.f3_unfilt_with_chrx_info[is.na(ase.f3_unfilt_with_chrx_info$Reported_XCI_status)]$Reported_XCI_status <- "Potential"
ase.f3_unfilt_with_chrx_info$status <- ifelse(ase.f3_unfilt_with_chrx_info$effectSize > 0.4, yes = "Monoallelic", no = "Biallelic")

ase.f3_unfilt$cell <- factor(ase.f3_unfilt$cell, levels = cell_order)
ase.f3_unfilt_with_chrx_info$cell <- factor(ase.f3_unfilt_with_chrx_info$cell, levels = cell_order)

chrx_gene_order <- fread("97_indexes_annotations/chrx_gene_order.tsv")

## Idenfity 'escape' genes that does not escape in thymocytes
id_genez <- ase.f3_unfilt_with_chrx_info[ase.f3_unfilt_with_chrx_info$keep == T & ase.f3_unfilt_with_chrx_info$Reported_XCI_status == "Escape",]

######### plot #########

library(ggplot2)
library(cowplot)
library(ggrepel)
library(ggrastr)


point_size <- 1.5

######## ASE for F1, F2 and F3 samples ###########

ase.f3_unfilt$chrx <- ifelse(ase.f3_unfilt$contig == "chrX", yes = TRUE, no = FALSE)
ase.f3_unfilt$chromosome <- "LOL"
ase.f3_unfilt[ase.f3_unfilt$contig == "chrX",]$chromosome <- "chrX"
ase.f3_unfilt[ase.f3_unfilt$contig != "chrX",]$chromosome <- "Autosome"

tezzt <- ase.f3_unfilt[keep == T & sample != "M1" & sample != "M2" & sample != "M3"] %>% dplyr::group_by(gene, sample, chrx, chromosome) %>% dplyr::summarise(alt_mean = mean(altCount), ref_mean = mean(refCount))




tezzt$chromosome <- factor(tezzt$chromosome, levels = c("Autosome","chrX"))
ase.f3_unfilt$chromosome <- "LOL"
ase.f3_unfilt[ase.f3_unfilt$contig == "chrX",]$chromosome <- "chrX"
ase.f3_unfilt[ase.f3_unfilt$contig != "chrX",]$chromosome <- "Autosome"

genes_to_label <- c("XIST", "TSIX", "CD99", "NOTCH1")

p.allelic.genomic.all_means <-
  tezzt %>% arrange(chromosome) %>%
  ggplot(aes(y= log10(ref_mean + 1), x= log10(alt_mean + 1), col=chromosome)) +
  geom_abline(slope=1, intercept=0) +
  geom_point_rast(data=tezzt[!tezzt$gene %in% genes_to_label,],size = point_size) +
  geom_point_rast(data=tezzt[!tezzt$gene %in% genes_to_label,],size = point_size) +
  geom_point_rast(data=tezzt[!tezzt$gene %in% genes_to_label & tezzt$chromosome == "Autosome",],size = point_size, color = "grey")+
  geom_point_rast(data=tezzt[!tezzt$gene %in% genes_to_label & tezzt$chromosome == "chrX",],size = point_size, color = "red")+
  geom_point_rast(data=tezzt[tezzt$gene %in% genes_to_label,],size = point_size, color = "black") +
  facet_wrap(sample~., scale = "free", ncol = 1) +
  labs(y="log10(ref+1)", x="log10(alt+1)") +
  theme_AL_simple() +
  theme( legend.position = "none", strip.background = element_blank(), axis.title.x = element_text(size=8), axis.title.y = element_text(size=8), strip.text.x = element_text(size = 8), axis.text = element_text(size=8))+
  scale_y_continuous(limits = c(0, 4), expand = c(0.05, 0), breaks = c(0,2,4))+
  scale_x_continuous(limits = c(0, 4), expand = c(0.05, 0), breaks = c(0,2,4))+
  geom_text_repel(data=tezzt[tezzt$gene %in% genes_to_label,], aes(label=gene), col = "black", nudge_x = 0.75, min.segment.length = 0.00000000001,  size=3)+
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"), strip.text.x = element_text(margin = margin(.01, 0, .01, 0, "cm")))


#standard a4 is width 210, height 297 <- with 20 mm margins -> width 170, height 257.

ggsave2("04_plots/Figure_2C.pdf",
        plot_grid(ncol=3, p.allelic.genomic.all_means, NULL, NULL))
