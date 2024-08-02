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

######### plot #########

library(ggplot2)
library(cowplot)
library(ggrepel)
library(ggrastr)


point_size <- 1.5

######## Biallelicly expressed genes across development ##########

goiz <- ase.f3_unfilt_with_chrx_info[ase.f3_unfilt_with_chrx_info$pass_all_filters == T,]$gene

test <- ase.f3_unfilt_with_chrx_info[ase.f3_unfilt_with_chrx_info$gene %in% goiz & ase.f3_unfilt_with_chrx_info$pass_all_filters == T,] %>% dplyr::group_by(gene, cell, Reported_XCI_status) %>% rstatix::get_summary_stats(effectSize, type = "common")

test$category <- test$Reported_XCI_status
test[test$gene %in% ase.f3_unfilt_with_chrx_info[ase.f3_unfilt_with_chrx_info$PAR == "PAR",]$gene,]$category <- "PAR"

library(ggsci)
library(scales)
show_col(pal_npg("nrc")(10))

p.biallelic.chrx.f3_main <-
  ggplot(test, aes(x=cell, y=mean, group=gene, col =  category)) +
  geom_hline(yintercept = 0.4, lty=2) +
  geom_line(aes(group=gene)) +
  geom_point(size = point_size)+
  geom_pointrange(aes(ymin=mean-sd,ymax=mean+sd), fatten = 0.75, show.legend = F)+
  facet_wrap(~factor(gene, levels = c(chrx_gene_order$x, "RPS6KA3")), ncol = 12) +
  theme_AL_box_rotX() +
  theme(strip.background = element_blank(), legend.position = "top", axis.title.x = element_blank(), axis.text.x = element_text(size = 8), axis.text.y = element_blank(), strip.text.x = element_text(size = 8), axis.title.y = element_text(size=8), axis.text = element_text(size=8))+
  labs(y="effectSize")+
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"), strip.text.x = element_text(margin = margin(.01, 0, .01, 0, "cm")))+
  scale_color_manual(values=c("Escape"="#DC0000FF","Inactive"="#869a9a", "Potential"="blue","PAR"="#00A087FF", "Variable"="black"))



ggsave2("04_plots/Suppl_Figure_3A.pdf", width =14, height = 14,
        plot_grid(p.biallelic.chrx.f3_main))




