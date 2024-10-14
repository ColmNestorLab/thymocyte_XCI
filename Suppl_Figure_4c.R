library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(ggrastr)
library(ggsci)

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
ase.f3_unfilt_with_chrx_info <- ase.f3_unfilt_with_chrx_info[!grepl(ase.f3_unfilt_with_chrx_info$gene, pattern = ";")]

# change cell order
ase.f3_unfilt_with_chrx_info$cell <- factor(ase.f3_unfilt_with_chrx_info$cell, levels = cell_order)

################# Small fixes ####################
ase.f3_unfilt_with_chrx_info[ase.f3_unfilt_with_chrx_info$Reported_XCI_status == "Unknown"]$Reported_XCI_status <- "Potential"
ase.f3_unfilt_with_chrx_info[is.na(ase.f3_unfilt_with_chrx_info$Reported_XCI_status)]$Reported_XCI_status <- "Potential"
ase.f3_unfilt_with_chrx_info$status <- ifelse(ase.f3_unfilt_with_chrx_info$effectSize > 0.4, yes = "Monoallelic", no = "Biallelic")

ase.f3_unfilt_with_chrx_info$cell <- factor(ase.f3_unfilt_with_chrx_info$cell, levels = cell_order)

chrx_gene_order <- fread("97_indexes_annotations/chrx_gene_order.tsv")

## Idenfity 'escape' genes that does not escape in thymocytes
id_genez <- ase.f3_unfilt_with_chrx_info[ase.f3_unfilt_with_chrx_info$keep == T & ase.f3_unfilt_with_chrx_info$Reported_XCI_status == "Escape",]

############## ASE count across development ################

ggsave2(filename = "04_plots/Suppl_Figure_4c.pdf",
        ggplot(ase.f3_unfilt_with_chrx_info[keep == T & contig == "chrX"], aes(x=cell, y=after_stat(count), fill = status)) + geom_bar(stat="count", position = "dodge") + geom_text(stat="count", position = position_dodge(width=.85), aes(label=..count..)) + labs(x="",title = "No. of mono/bial genes") + scale_fill_npg() + theme_AL_box() + theme(legend.title = element_blank())
)

# Export source data
write.table(ase.f3_unfilt_with_chrx_info[keep == T & contig == "chrX"], file = "06_source_data/Suppl_Fig_4c.tsv", quote = F, row.names = F, sep = "\t")

