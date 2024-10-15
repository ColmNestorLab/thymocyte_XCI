library(dplyr)
library(data.table)
library(readxl)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(ggrastr)

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
ase.f3_unfilt <- ase.f3_unfilt[!grepl(ase.f3_unfilt$gene, pattern = ";")]
ase.f3_unfilt_with_chrx_info <- ase.f3_unfilt_with_chrx_info[!grepl(ase.f3_unfilt_with_chrx_info$gene, pattern = ";")]

# change celltype order
ase.f3_unfilt$celltype <- factor(ase.f3_unfilt$celltype, levels = cell_order)
ase.f3_unfilt_with_chrx_info$celltype <- factor(ase.f3_unfilt_with_chrx_info$celltype, levels = cell_order)

################# Small fixes ####################
ase.f3_unfilt_with_chrx_info[ase.f3_unfilt_with_chrx_info$Reported_XCI_status == "Unknown"]$Reported_XCI_status <- "Potential"
ase.f3_unfilt_with_chrx_info[is.na(ase.f3_unfilt_with_chrx_info$Reported_XCI_status)]$Reported_XCI_status <- "Potential"
ase.f3_unfilt_with_chrx_info$status <- ifelse(ase.f3_unfilt_with_chrx_info$effectSize > 0.4, yes = "Monoallelic", no = "Biallelic")

# make sure cell is ordered factor
ase.f3_unfilt$celltype <- factor(ase.f3_unfilt$celltype, levels = cell_order)
ase.f3_unfilt_with_chrx_info$celltype <- factor(ase.f3_unfilt_with_chrx_info$celltype, levels = cell_order)

# get X-linked gene order
chrx_gene_order <- fread("97_indexes_annotations/chrx_gene_order.tsv")


####### Check unfiltered inactive genes throughout development ##########
# get inactive genes
unfilt_inac <- ase.f3_unfilt_with_chrx_info[ase.f3_unfilt_with_chrx_info$Reported_XCI_status == "Inactive" & ase.f3_unfilt_with_chrx_info$gene != "SEPT6" & ase.f3_unfilt_with_chrx_info$variant_totalCount > 20 & ase.f3_unfilt_with_chrx_info$totalCount >= 10,]

# select only columns we want
unfilt_inac_short <- dplyr::select(unfilt_inac, altCount, refCount, gene, celltype, position)

# convert altCount & refCount to minor and major allele (minor allele = fewest reads, major allele = most reads).
unfilt_inac_short_min <- transform(unfilt_inac_short, min = pmin(altCount, refCount))
unfilt_inac_short_min_max <- transform(unfilt_inac_short_min, max = pmax(altCount, refCount))

# melt, remove columns not used
melt_unfilt_inac_short_min_max <- melt(dplyr::select(unfilt_inac_short_min_max, -altCount, -refCount), id.vars = c("gene", "celltype", "position"), measure.vars = c("min", "max"))

# get sums
melt_unfilt_inac_short_min_max_countz <- melt_unfilt_inac_short_min_max %>% dplyr::group_by(celltype, variable) %>% dplyr::summarise(asdf = sum(value))

# dcast
melt_unfilt_inac_short_min_max_countz_dcast <- dcast(melt_unfilt_inac_short_min_max_countz, celltype ~ variable, value.var = "asdf")

# calculate total read count
melt_unfilt_inac_short_min_max_countz_dcast$tot <- melt_unfilt_inac_short_min_max_countz_dcast$min + melt_unfilt_inac_short_min_max_countz_dcast$max


# plot
fig_4c_left <- ggplot(melt_unfilt_inac_short_min_max[melt_unfilt_inac_short_min_max$variable == "min" & melt_unfilt_inac_short_min_max$value != 0,], aes(x=gene,y=value, col = celltype)) + geom_point() + theme_AL_box_rotX() + facet_grid(~celltype, scales = "free_x", space = "free_x") + coord_cartesian(ylim=c(0,NA)) + labs(x="inactive genes from figure 4b with >0 reads from Xi", y="Xi read count")






####### Check immune genes throughout development ##########
# Read in compiled X-linked immune gene lists (see methods for which articles they are from).
gois_1 <- fread("97_indexes_annotations/immune_genes_nri2815.tsv", header = F)
gois_2 <- fread("97_indexes_annotations/immune_genes_s13293-019-0278.tsv", header = T)
gois_3 <- fread("97_indexes_annotations/Immune_chrX_genes_fimmu.2017.01455.tsv", header = T)

# merge the gene names
genez <- unique(sort(c(gois_1$V1, gois_2$Gene_name, gois_3$gene)))

# make df with only the genes above.
unfilt_immune <- ase.f3_unfilt_with_chrx_info[ase.f3_unfilt_with_chrx_info$gene %in% genez & ase.f3_unfilt_with_chrx_info$variant_totalCount > 20 & ase.f3_unfilt_with_chrx_info$totalCount >= 10 & ase.f3_unfilt_with_chrx_info$effectSize > 0.4,]

# select only columns we want
unfilt_immune_short <- dplyr::select(unfilt_immune, altCount, refCount, gene, celltype, position)

# convert altCount & refCount to minor and major allele (minor allele = fewest reads, major allele = most reads).
unfilt_immune_short_min <- transform(unfilt_immune_short, min = pmin(altCount, refCount))
unfilt_immune_short_min_max <- transform(unfilt_immune_short_min, max = pmax(altCount, refCount))

# melt, remove columns not used
melt_unfilt_immune_short_min_max <- melt(dplyr::select(unfilt_immune_short_min_max, -altCount, -refCount), id.vars = c("gene", "celltype", "position"), measure.vars = c("min", "max"))

# get sums
melt_unfilt_immune_short_min_max_countz <- melt_unfilt_immune_short_min_max %>% dplyr::group_by(celltype, variable) %>% dplyr::summarise(asdf = sum(value))

# dcast
melt_unfilt_immune_short_min_max_countz_dcast <- dcast(melt_unfilt_immune_short_min_max_countz, celltype ~ variable, value.var = "asdf")

# calculate total read count
melt_unfilt_immune_short_min_max_countz_dcast$tot <- melt_unfilt_immune_short_min_max_countz_dcast$min + melt_unfilt_immune_short_min_max_countz_dcast$max

#plot
fig_4c_right <- ggplot(melt_unfilt_immune_short_min_max[melt_unfilt_immune_short_min_max$variable == "min" & melt_unfilt_immune_short_min_max$value != 0,], aes(x=gene,y=value, col = celltype)) + geom_point() + theme_AL_box_rotX() + facet_grid(~celltype, scales = "free_x", space = "free_x") + coord_cartesian(ylim=c(0,NA)) + labs(x="inactive genes from figure 4b with >0 reads from Xi", y="Xi read count")



# plot
immune_plot <- ggplot(melt_unfilt_immune_short_min_max_countz_dcast, aes(x=celltype, y=min, label = paste0("total: ", tot))) + geom_bar(stat="identity") + labs(x="",y="minor allele count", title = paste0(length(genez) ," immune genes; ",length(unique(melt_unfilt_immune_short_min_max$position))
," hSNPs in ",length(unique(melt_unfilt_immune_short_min_max$gene))
," genes, wes>20, rnaseq>=10")) + theme_AL_box() + geom_text(aes(y=10.5)) + geom_text(data=melt_unfilt_immune_short_min_max_countz_dcast, aes(label=min))

# plot
inac_plot <- ggplot(melt_unfilt_inac_short_min_max_countz_dcast, aes(x=celltype, y=min, label = paste0("total: ", tot))) + geom_bar(stat="identity") + labs(x="",y="minor allele count", title = paste0("inactive (not SEPT6); ", length(unique(unfilt_inac$position)), " hSNPs in ",length(unique(unfilt_inac$gene)), " genes, wes>20, rnaseq>=10")) + theme_AL_box() + geom_text(aes(y=22.5)) + geom_text(data=melt_unfilt_inac_short_min_max_countz_dcast, aes(label=min))


#plot
ggsave2(filename = "04_plots/Figure_4b.pdf", height = 8, width=14,
        plot_grid(ncol=2, inac_plot, immune_plot))

ggsave2(filename = "04_plots/Figure_4c.pdf", height = 4, width=10,
        plot_grid(ncol=2, fig_4c_left, fig_4c_right))


#Export source data
melt_unfilt_inac_short_min_max_countz_dcast$tag <- "inactive_genes"
melt_unfilt_immune_short_min_max_countz_dcast$tag <- "X_linked_immune_genes"

fig4c_left <- melt_unfilt_inac_short_min_max[melt_unfilt_inac_short_min_max$variable == "min" & melt_unfilt_inac_short_min_max$value != 0,]
fig4c_left$tag <- "inactive_genes"

fig4c_right <- melt_unfilt_immune_short_min_max[melt_unfilt_immune_short_min_max$variable == "min" & melt_unfilt_immune_short_min_max$value != 0,]
fig4c_right$tag <- "X_linked_immune_genes"

write.table(rbind(melt_unfilt_inac_short_min_max_countz_dcast, melt_unfilt_immune_short_min_max_countz_dcast), file = "06_source_data/Fig_4b.tsv", quote = F, row.names = F, sep = "\t")
write.table(rbind(fig4c_left, fig4c_right), file = "06_source_data/Fig_4c.tsv", quote = F, row.names = F, sep = "\t")

