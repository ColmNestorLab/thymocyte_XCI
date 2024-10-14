
library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(ggrastr)
library(ggsci)
library(scales)

source("plot_parameters.R")

point_size <- 1.5

####### Read in data ########
ase_unfilt <- fread("02_tidy_data/Suppl_tables/Supplementary_data_ASE_table_all_chroms_and_samples.tsv")

# get F3 data
ase.f3_unfilt_with_chrx_info <- ase_unfilt[!is.na(ase_unfilt$gene) & sample == "F3" & contig == "chrX",]

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

# change celltype order
ase_unfilt$celltype <- factor(ase_unfilt$celltype, levels = cell_order)
ase.f3_unfilt_with_chrx_info$celltype <- factor(ase.f3_unfilt_with_chrx_info$celltype, levels = cell_order)

################# Small fixes ####################
ase.f3_unfilt_with_chrx_info[ase.f3_unfilt_with_chrx_info$Reported_XCI_status == "Unknown"]$Reported_XCI_status <- "Potential"
ase.f3_unfilt_with_chrx_info[is.na(ase.f3_unfilt_with_chrx_info$Reported_XCI_status)]$Reported_XCI_status <- "Potential"
ase.f3_unfilt_with_chrx_info$status <- ifelse(ase.f3_unfilt_with_chrx_info$effectSize > 0.4, yes = "Monoallelic", no = "Biallelic")

# make sure cell is ordered factor
ase_unfilt$celltype <- factor(ase_unfilt$celltype, levels = cell_order)
ase.f3_unfilt_with_chrx_info$celltype <- factor(ase.f3_unfilt_with_chrx_info$celltype, levels = cell_order)

# get X-linked gene order
chrx_gene_order <- fread("97_indexes_annotations/chrx_gene_order.tsv")



######## Biallelicly expressed genes across development ##########
# get biallelic genes
goiz <- ase.f3_unfilt_with_chrx_info[ase.f3_unfilt_with_chrx_info$effectSize <= 0.4 & ase.f3_unfilt_with_chrx_info$pass_all_filters == T,]$gene

# get means of escapees
means_ase.f3_escapees <- ase.f3_unfilt_with_chrx_info[ase.f3_unfilt_with_chrx_info$gene %in% goiz & ase.f3_unfilt_with_chrx_info$pass_all_filters == T,] %>% dplyr::group_by(gene, celltype, Reported_XCI_status) %>% rstatix::get_summary_stats(effectSize, type = "common")

# Find non-escaping escapees
#asdfff <- ase.f3_unfilt_with_chrx_info[ase.f3_unfilt_with_chrx_info$Reported_XCI_status == "Escape" & ase.f3_unfilt_with_chrx_info$pass_all_filters == T,]
# ARSD, GEMIN8, MSL3 does not escape in thymocytes
non_esc_escapees <- c("ARSD", "GEMIN8", "MSL3")

# find some inactive genes to add
#asdfff <- ase.f3_unfilt_with_chrx_info[ase.f3_unfilt_with_chrx_info$Reported_XCI_status == "Inactive" & ase.f3_unfilt_with_chrx_info$pass_all_filters == T,]
inactive_genezz <- c("ITM2A", "KIF4A", "PHF6")

# get means of non-escapees and inactives
means_ase.f3_non_escapees_inactive_examples <- ase.f3_unfilt_with_chrx_info[ase.f3_unfilt_with_chrx_info$gene %in% c(non_esc_escapees, inactive_genezz) & ase.f3_unfilt_with_chrx_info$pass_all_filters == T,] %>% dplyr::group_by(gene, celltype, Reported_XCI_status) %>% rstatix::get_summary_stats(effectSize, type = "common")

# make sure celltype are factors
means_ase.f3_non_escapees_inactive_examples$celltype <- factor(means_ase.f3_non_escapees_inactive_examples$celltype, levels = cell_order)

# get Tukiainen escape status but single out non-escape and inactive genes
means_ase.f3_non_escapees_inactive_examples$category <- means_ase.f3_non_escapees_inactive_examples$Reported_XCI_status
means_ase.f3_non_escapees_inactive_examples[means_ase.f3_non_escapees_inactive_examples$gene %in% non_esc_escapees,]$category <- "Non-escape"
means_ase.f3_non_escapees_inactive_examples[means_ase.f3_non_escapees_inactive_examples$gene %in% inactive_genezz,]$category <- "Inactive"


# make sure celltype are factors
means_ase.f3_escapees$celltype <- factor(means_ase.f3_escapees$celltype, levels = cell_order)

# get Tukiainen escape status but single out non-escape and inactive genes
means_ase.f3_escapees$category <- means_ase.f3_escapees$Reported_XCI_status
means_ase.f3_escapees[means_ase.f3_escapees$gene %in% PAR$GENES,]$category <- "PAR"

# get RAG1, GATA3, NOTCH1 (autosomal genes)
find_auto <- ase_unfilt[ase_unfilt$sample == "F3",]
autosomal_df <- find_auto[find_auto$gene %in% c("RAG1", "GATA3", "RORC") & find_auto$pass_all_filters == T,] %>% dplyr::group_by(gene, celltype) %>% rstatix::get_summary_stats(effectSize, type = "common")
autosomal_df$Reported_XCI_status <- "Autosomal"
autosomal_df$category <- "Autosomal"

# bind together
test_smash <- rbind(means_ase.f3_escapees, autosomal_df, means_ase.f3_non_escapees_inactive_examples)

# make gene order for the plot
genez_order <- c(unique(autosomal_df$gene), unique(means_ase.f3_escapees[means_ase.f3_escapees$category == "PAR",]$gene), unique(means_ase.f3_escapees[means_ase.f3_escapees$category == "Escape",]$gene), unique(means_ase.f3_escapees[means_ase.f3_escapees$category == "Variable",]$gene), unique(means_ase.f3_escapees[means_ase.f3_escapees$category == "Potential",]$gene), unique(means_ase.f3_escapees[means_ase.f3_escapees$category == "Inactive",]$gene), unique(means_ase.f3_non_escapees_inactive_examples[means_ase.f3_non_escapees_inactive_examples$category == "Inactive",]$gene), unique(means_ase.f3_non_escapees_inactive_examples[means_ase.f3_non_escapees_inactive_examples$category == "Non-escape",]$gene))



# plot
#show_col(pal_npg("nrc")(10))

p.biallelic.chrx.f3_main <-
  ggplot(test_smash, aes(x=celltype, y=mean, group=gene, col =  category)) +
  geom_hline(yintercept = 0.4, lty=2) +
  geom_line(aes(group=gene)) +
  geom_point(size = point_size)+
  geom_pointrange(aes(ymin=mean-sd,ymax=mean+sd), fatten = 0.75, show.legend = F)+
  facet_wrap(~factor(gene, levels = genez_order), ncol = 10) +
  theme_AL_box_rotX() +
  theme(strip.background = element_blank(), legend.position = "top", axis.title.x = element_blank(), axis.text.x = element_text(size = 8), axis.text.y = element_blank(), strip.text.x = element_text(size = 8), axis.title.y = element_text(size=8), axis.text = element_text(size=8))+
  labs(y="effectSize")+
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"), strip.text.x = element_text(margin = margin(.01, 0, .01, 0, "cm")))+
  scale_color_manual(values=c("Escape"="#DC0000FF","Inactive"="#869a9a", "Potential"="blue","PAR"="#00A087FF", "Autosomal"="black", "Non-escape" = "#7E6148FF"))



ggsave2("04_plots/Figure_3e.pdf", width =14,
        plot_grid(p.biallelic.chrx.f3_main))


# export source data
write.table(test_smash, file = "06_source_data/Fig_3e.tsv", quote = F, row.names = F, sep = "\t")
