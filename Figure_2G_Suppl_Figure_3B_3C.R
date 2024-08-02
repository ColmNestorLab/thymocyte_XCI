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


####### Check unfiltered inactive genes throughout development ##########

unfilt_inac <- ase.f3_unfilt_with_chrx_info[ase.f3_unfilt_with_chrx_info$Reported_XCI_status == "Inactive" & ase.f3_unfilt_with_chrx_info$gene != "SEPT6" & ase.f3_unfilt_with_chrx_info$variant_totalCount > 20 & ase.f3_unfilt_with_chrx_info$totalCount >= 10,]

length(unique(unfilt_inac$position))
length(unique(unfilt_inac$gene))

range(unfilt_inac$totalCount)
range(unfilt_inac$refCount)
range(unfilt_inac$altCount)

leeeel <- dplyr::select(unfilt_inac, altCount, refCount, gene, cell, position)

leeeel_min <- transform(leeeel, min = pmin(altCount, refCount))
leeeel_min_max <- transform(leeeel_min, max = pmax(altCount, refCount))

melt_leeeel_min_max <- melt(dplyr::select(leeeel_min_max, -altCount, -refCount), id.vars = c("gene", "cell", "position"), measure.vars = c("min", "max"))

melt_leeeel_min_max_countz <- melt_leeeel_min_max %>% dplyr::group_by(cell, variable) %>% dplyr::summarise(asdf = sum(value))

melt_leeeel_min_max_countz_dcast <- dcast(melt_leeeel_min_max_countz, cell ~ variable, value.var = "asdf")

melt_leeeel_min_max_countz_dcast$tot <- melt_leeeel_min_max_countz_dcast$min + melt_leeeel_min_max_countz_dcast$max

#######################

asdfff <- ggplot(melt_leeeel_min_max[melt_leeeel_min_max$variable == "min" & melt_leeeel_min_max$value != 0,], aes(x=gene,y=value, col = cell)) + geom_point() + theme_AL_box_rotX() + facet_grid(~cell, scales = "free_x", space = "free_x") + coord_cartesian(ylim=c(0,NA)) + labs(x="inactive genes from figure 2G with >0 reads from Xi", y="Xi read count")






####### Check immune genes throughout development ##########

gois_1 <- fread("97_indexes_annotations/immune_genes_nri2815.tsv", header = F)
gois_2 <- fread("97_indexes_annotations/immune_genes_s13293-019-0278.tsv", header = T)
gois_3 <- fread("97_indexes_annotations/Immune_chrX_genes_fimmu.2017.01455.tsv", header = T)

genez <- unique(sort(c(gois_1$V1, gois_2$Gene_name, gois_3$gene)))

unfilt_immune <- ase.f3_unfilt_with_chrx_info[ase.f3_unfilt_with_chrx_info$gene %in% genez & ase.f3_unfilt_with_chrx_info$variant_totalCount > 20 & ase.f3_unfilt_with_chrx_info$totalCount >= 10 & ase.f3_unfilt_with_chrx_info$effectSize > 0.4,]

ggplot(unfilt_immune, aes(x=cell, y = effectSize)) + facet_wrap(~gene) + geom_line(aes(group=gene)) + coord_cartesian(ylim=c(0,0.5))

range(unfilt_immune$totalCount)
range(unfilt_immune$refCount)
range(unfilt_immune$altCount)

leeeel2 <- dplyr::select(unfilt_immune, altCount, refCount, gene, cell, position)

leeeel_min2 <- transform(leeeel2, min = pmin(altCount, refCount))
leeeel_min_max2 <- transform(leeeel_min2, max = pmax(altCount, refCount))

melt_leeeel_min_max2 <- melt(dplyr::select(leeeel_min_max2, -altCount, -refCount), id.vars = c("gene", "cell", "position"), measure.vars = c("min", "max"))

melt_leeeel_min_max_countz2 <- melt_leeeel_min_max2 %>% dplyr::group_by(cell, variable) %>% dplyr::summarise(asdf = sum(value))


length(unique(melt_leeeel_min_max2$position))
length(unique(melt_leeeel_min_max2$gene))


#stoopid2 <- as.data.frame(dcast(leeeel_min_max2, gene+position~cell, value.var = c("min", "max")))

#stoopid2[is.na(stoopid2)] <- 0

melt_leeeel_min_max_countz_dcast2 <- dcast(melt_leeeel_min_max_countz2, cell ~ variable, value.var = "asdf")

melt_leeeel_min_max_countz_dcast2$tot <- melt_leeeel_min_max_countz_dcast2$min + melt_leeeel_min_max_countz_dcast2$max

asdfff2 <- ggplot(melt_leeeel_min_max2[melt_leeeel_min_max2$variable == "min" & melt_leeeel_min_max2$value != 0,], aes(x=gene,y=value, col = cell)) + geom_point() + theme_AL_box_rotX() + facet_grid(~cell, scales = "free_x", space = "free_x") + coord_cartesian(ylim=c(0,NA)) + labs(x="inactive genes from figure 2G with >0 reads from Xi", y="Xi read count")


#######################
immune_plot <- ggplot(melt_leeeel_min_max_countz_dcast2, aes(x=cell, y=min, label = paste0("total: ", tot))) + geom_bar(stat="identity") + labs(x="",y="minor allele count", title = paste0(length(genez) ," immune genes; ",length(unique(melt_leeeel_min_max2$position))
," hSNPs in ",length(unique(melt_leeeel_min_max2$gene))
," genes, wes>20, rnaseq>=10")) + theme_AL_box() + geom_text(aes(y=10.5)) + geom_text(data=melt_leeeel_min_max_countz_dcast2, aes(label=min))

inac_plot <- ggplot(melt_leeeel_min_max_countz_dcast, aes(x=cell, y=min, label = paste0("total: ", tot))) + geom_bar(stat="identity") + labs(x="",y="minor allele count", title = paste0("inactive (not SEPT6); ", length(unique(unfilt_inac$position)), " hSNPs in ",length(unique(unfilt_inac$gene)), " genes, wes>20, rnaseq>=10")) + theme_AL_box() + geom_text(aes(y=22.5)) + geom_text(data=melt_leeeel_min_max_countz_dcast, aes(label=min))

library(cowplot)
ggsave2(filename = "04_plots/Figure_2G.pdf", height = 8, width=14,
        plot_grid(ncol=2, inac_plot, immune_plot))



ggsave2(filename = "04_plots/SFig_3b.pdf", height = 4, width=10,
        plot_grid(ncol=1, asdfff))

#ggsave2(filename = "04_plots/Figure_2H_immune_genes_points.pdf", height = 4, width=10,
#        plot_grid(ncol=1, asdfff2))


TPMz <- melt(fread("Supplementary Table 3"))

TPMz <- TPMz[TPMz$gene %in% genez]

TPMz <- TPMz[grepl(TPMz$variable, pattern = "F3")] %>% separate(variable, into = c("sample", "cell", NA, NA))

TPMz$expressed <- ifelse(TPMz$value > 1, yes = "TPM > 1", no = "TPM < 1")

TPMz_countz <- TPMz %>% dplyr::group_by(expressed, cell) %>% dplyr::count()

TPMz_countz2 <- TPMz %>% dplyr::group_by(gene) %>% dplyr::count()


ggsave2(filename = "04_plots/SFig_3c.pdf", height = 6, width=5,
ggplot(TPMz_countz, aes(x=factor(cell, levels = cell_order), y=n, fill = expressed, label = n)) + geom_bar(stat="identity", position = "dodge")+ geom_text(position = position_dodge(width=0.9)) + theme_AL_box_rotX() + labs(x="", y = "immune gene count", title = paste0("73 out of the ", length(genez) ," genes is detected"))
)
