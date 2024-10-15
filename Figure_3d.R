library(ggplot2)
library(cowplot)
library(ggrepel)
library(ggrastr)
library(dplyr)
library(data.table)
library(readxl)

source("plot_parameters.R")

point_size <- 1.5

####### Read in data ########
ase.f3_unfilt <- data.table(read_excel("Supplementary Data.xlsx", sheet = 8, skip = 1, col_types = c("text", "numeric", "text", "text", "text", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "logical", "logical")))

# remove NA genes
ase.f3_unfilt_with_chrx_info <- ase.f3_unfilt[!is.na(ase.f3_unfilt$gene)]

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
ase.f3_unfilt_with_chrx_info <- ase.f3_unfilt_with_chrx_info[!grepl(ase.f3_unfilt_with_chrx_info$gene, pattern = ";")]

# change cell order
ase.f3_unfilt_with_chrx_info$cell <- factor(ase.f3_unfilt_with_chrx_info$cell, levels = cell_order)

################# Small fixes ####################
ase.f3_unfilt_with_chrx_info[ase.f3_unfilt_with_chrx_info$Reported_XCI_status == "Unknown"]$Reported_XCI_status <- "Potential"
ase.f3_unfilt_with_chrx_info[is.na(ase.f3_unfilt_with_chrx_info$Reported_XCI_status)]$Reported_XCI_status <- "Potential"
ase.f3_unfilt_with_chrx_info$status <- ifelse(ase.f3_unfilt_with_chrx_info$effectSize > 0.4, yes = "Monoallelic", no = "Biallelic")

# make sure cell is ordered factor
ase.f3_unfilt_with_chrx_info$cell <- factor(ase.f3_unfilt_with_chrx_info$cell, levels = cell_order)

# get X-linked gene order
chrx_gene_order <- fread("97_indexes_annotations/chrx_gene_order.tsv")

############## AE across genomic positions on chromosome X ################

# make sure its filtered on highest covered SNP and only X-linked genes
genomic_positions_figure <- ase.f3_unfilt_with_chrx_info[keep == T & contig == "chrX"]

# get Tukiainens XCI status classificatio but separate out PAR genes
genomic_positions_figure$category <- genomic_positions_figure$Reported_XCI_status
genomic_positions_figure[genomic_positions_figure$PAR == "PAR",]$category <- "PAR"

# plot
p.genome.chrx.f3 <-
  ggplot(genomic_positions_figure, aes(x=position/1e6, y=effectSize, col = factor(category, levels = c("PAR", "Escape", "Variable", "Potential", "Inactive")))) +
  geom_hline(yintercept = 0.4, lty=2) +
  geom_point(size = point_size) +
  annotate("rect", xmin=c(10001,155701383)/1e6, xmax=c(2781479, 156030895)/1e6, ymin=c(-0.01,-0.01),ymax=c(0,0)) +
  labs(x="Chromosome X (Mbp)") +
  facet_wrap(~cell, nrow = 2, scales = "free_x") +
  scale_color_manual(values=c("PAR"="#00A087FF", "Escape"="#DC0000FF","Variable"="purple","Potential"="blue","Inactive"="#869a9a"))+
  theme_AL_simple() +
  theme(legend.position = "top",strip.background = element_blank(), legend.title = element_blank())+
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"), strip.text.x = element_text(margin = margin(.01, 0, .01, 0, "cm")))

ggsave2("04_plots/Figure_3d.pdf",
        plot_grid(p.genome.chrx.f3))


# Export source data
write.table(genomic_positions_figure, file = "06_source_data/Fig_3d.tsv", quote = F, row.names = F, sep = "\t")

