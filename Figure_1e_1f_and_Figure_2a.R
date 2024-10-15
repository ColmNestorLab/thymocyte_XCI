library(dplyr)
library(RColorBrewer)
library(rtracklayer)
library(cowplot)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(ggrepel)
library(data.table)
library(ggbeeswarm)
library(rstatix)
library(readxl)

source("plot_parameters.R")

## Read in sleuth output tables
df_sleuth_across <- data.table(read_excel("Supplementary Data.xlsx", sheet = 4, skip = 1))
df_sleuth_split <- data.table(read_excel("Supplementary Data.xlsx", sheet = 5, skip = 1))

# make numeric
df_sleuth_across$log2FC <- as.numeric(df_sleuth_across$log2FC)
df_sleuth_split$log2FC <- as.numeric(df_sleuth_split$log2FC)

# read in TPM data
raw_tpm <- data.table(read_excel("Supplementary Data.xlsx", sheet = 2, skip = 1, col_types = c("text",rep(paste0("numeric"), 46))))

#melt
melt_raw_tpm <- melt(raw_tpm, id.vars = "gene")
melt_raw_tpm$value <- as.numeric(melt_raw_tpm$value)

# change identifiers
melt_raw_tpm$sex <- "KEK"
melt_raw_tpm[grepl(melt_raw_tpm$variable, pattern = "F1"),]$sex <- "Female"
melt_raw_tpm[grepl(melt_raw_tpm$variable, pattern = "F2"),]$sex <- "Female"
melt_raw_tpm[grepl(melt_raw_tpm$variable, pattern = "F3"),]$sex <- "Female"
melt_raw_tpm[grepl(melt_raw_tpm$variable, pattern = "F4"),]$sex <- "Female"
melt_raw_tpm[grepl(melt_raw_tpm$variable, pattern = "M1"),]$sex <- "Male"
melt_raw_tpm[grepl(melt_raw_tpm$variable, pattern = "M3"),]$sex <- "Male"
melt_raw_tpm[grepl(melt_raw_tpm$variable, pattern = "M4"),]$sex <- "Male"

melt_raw_tpm$cell <- "KEK"
melt_raw_tpm[grepl(melt_raw_tpm$variable, pattern = "ETP"),]$cell <- "ETP"
melt_raw_tpm[grepl(melt_raw_tpm$variable, pattern = "DN"),]$cell <- "DN"
melt_raw_tpm[grepl(melt_raw_tpm$variable, pattern = "DPearly"),]$cell <- "DPearly"
melt_raw_tpm[grepl(melt_raw_tpm$variable, pattern = "DPlate"),]$cell <- "DPlate"
melt_raw_tpm[grepl(melt_raw_tpm$variable, pattern = "CD4SP"),]$cell <- "CD4SP"
melt_raw_tpm[grepl(melt_raw_tpm$variable, pattern = "CD8SP"),]$cell <- "CD8SP"

# summarize, get means and keep only genes in our sleuth df
melt_raw_tpm_across <- melt_raw_tpm %>% dplyr::group_by(gene) %>% dplyr::summarise(meanz = mean(value))
melt_raw_tpm_across <- melt_raw_tpm_across[melt_raw_tpm_across$gene %in% df_sleuth_across$gene,]

melt_raw_tpm_means_split <- melt_raw_tpm %>% dplyr::group_by(cell, gene) %>% dplyr::summarise(meanz = mean(value))
melt_raw_tpm_means_split <- melt_raw_tpm_means_split[melt_raw_tpm_means_split$gene %in% df_sleuth_across$gene,]

# merge means with sleuth dfs, keeping only chrX genes
df_sleuth_across_with_TPM <- subset(merge(df_sleuth_across, melt_raw_tpm_across, by = c("gene"), all.x =T), contig == "chrX")
df_sleuth_split_with_TPM <- subset(merge(df_sleuth_split, melt_raw_tpm_means_split, by.x = c("gene", "celltype"), by.y = c("gene", "cell"), all.x =T), contig == "chrX")

# add annotaitons
# tukiainen et al annotation
x.esc <- dplyr::select(fread("97_indexes_annotations/landscape.Suppl.Table.13.csv"), Gene_name, Reported_XCI_status, Sex_bias_in_GTEx, XCI_across_tissues, XCI_in_single_cells)
esc <- dplyr::select(subset(x.esc), Gene_name, Reported_XCI_status)
esc <- esc[esc$Gene_name != ""]
esc[esc$Gene_name %in% "6-Sep"]$Gene_name <- "SEPT6"
esc <- esc[!(esc$Gene_name %in% "IDS" & esc$Reported_XCI_status == "Unknown")]

#PAR
PAR_df <- rbind(fread("97_indexes_annotations/group-715_PAR1.csv"), fread("97_indexes_annotations/group-716_PAR2.csv"))
PAR_df$GENES <- PAR_df$`Approved symbol`

# merge annotations
chrx_annotation <- merge(esc, dplyr::select(PAR_df, 7,8), by.x = "Gene_name", by.y = "GENES", all = T)
chrx_annotation$category <- chrx_annotation$Reported_XCI_status
chrx_annotation[chrx_annotation$PAR == "PAR1" | chrx_annotation$PAR == "PAR2",]$category <- "PAR"

# merge our dfs with the chrx annotation.
df_sleuth_across_with_TPM_anno <- merge(df_sleuth_across_with_TPM, chrx_annotation, by.x = c("gene"), by.y = c("Gene_name"), all.x = T)
df_sleuth_split_with_TPM_anno <- merge(df_sleuth_split_with_TPM, chrx_annotation, by.x = c("gene"), by.y = c("Gene_name"), all.x = T)

# make NAs and Unknown into Potential
df_sleuth_across_with_TPM_anno[is.na(df_sleuth_across_with_TPM_anno$category),]$category <- "Potential"
df_sleuth_across_with_TPM_anno[df_sleuth_across_with_TPM_anno$category == "Unknown",]$category <- "Potential"

df_sleuth_split_with_TPM_anno[is.na(df_sleuth_split_with_TPM_anno$category),]$category <- "Potential"
df_sleuth_split_with_TPM_anno[df_sleuth_split_with_TPM_anno$category == "Unknown",]$category <- "Potential"


# Add levels to category
df_sleuth_across_with_TPM_anno$category <- factor(df_sleuth_across_with_TPM_anno$category, levels = c("PAR", "Escape", "Variable","Inactive", "Potential"))
df_sleuth_split_with_TPM_anno$category <- factor(df_sleuth_split_with_TPM_anno$category, levels = c("PAR", "Escape",  "Variable", "Inactive", "Potential"))

# plot

# set cap, save extreme b values
extreme_vals_across <- df_sleuth_across_with_TPM_anno[df_sleuth_across_with_TPM_anno$log2FC > 1 | df_sleuth_across_with_TPM_anno$log2FC < -1]
extreme_vals_split <- df_sleuth_split_with_TPM_anno[df_sleuth_split_with_TPM_anno$log2FC > 1 | df_sleuth_split_with_TPM_anno$log2FC < -1]
df_sleuth_across_with_TPM_anno[df_sleuth_across_with_TPM_anno$log2FC > 1]$log2FC <- 1
df_sleuth_across_with_TPM_anno[df_sleuth_across_with_TPM_anno$log2FC < -1]$log2FC <- -1
df_sleuth_split_with_TPM_anno[df_sleuth_split_with_TPM_anno$log2FC > 1]$log2FC <- 1
df_sleuth_split_with_TPM_anno[df_sleuth_split_with_TPM_anno$log2FC < -1]$log2FC <- -1

# calc stuff plot adding to the plot for selecting, post-squish.
statz_category_squish_across <- df_sleuth_across_with_TPM_anno %>% group_by(category) %>% get_summary_stats(log2FC, type = "common")
setDT(statz_category_squish_across)

statz_category_squish_split <- df_sleuth_split_with_TPM_anno %>%   group_by(category, celltype) %>%  get_summary_stats(log2FC, type = "common")
setDT(statz_category_squish_split)

# get values for easier plotting further on.
median <- statz_category_squish_across[statz_category_squish_across$category == "Inactive"]$median
upper <- statz_category_squish_across[statz_category_squish_across$category == "Inactive"]$median + statz_category_squish_across[statz_category_squish_across$category == "Inactive"]$iqr
lower <- statz_category_squish_across[statz_category_squish_across$category == "Inactive"]$median - statz_category_squish_across[statz_category_squish_across$category == "Inactive"]$iqr

# same as above, store values for easier plotting.
cells_PAR <- statz_category_squish_split[statz_category_squish_split$category == "PAR"]$celltype
median_PAR <- statz_category_squish_split[statz_category_squish_split$category == "PAR"]$median
upper_PAR <- statz_category_squish_split[statz_category_squish_split$category == "PAR"]$median + statz_category_squish_split[statz_category_squish_split$category == "PAR"]$iqr
lower_PAR <- statz_category_squish_split[statz_category_squish_split$category == "PAR"]$median - statz_category_squish_split[statz_category_squish_split$category == "PAR"]$iqr
PAR_lines <- data.frame(celltype = cells_PAR, median = median_PAR, upper = upper_PAR, lower = lower_PAR)

# get the order of X-linked genes for nice plots.
chrx_gene_order <- unique(subset(fread("97_indexes_annotations/GCF_000001405.38_GRCh38.p12_genomic_with_correct_contig_names.tsv"), seqnames == "chrX")$gene_id)


# plot

fig1e <- 
  ggplot(data = df_sleuth_across_with_TPM_anno[df_sleuth_across_with_TPM_anno$meanz > 1], aes(x = gene, y = log2FC)) + 
  facet_grid(~category, scales = "free_x", space = "free") +  
  geom_pointrange(data = df_sleuth_across_with_TPM_anno[(df_sleuth_across_with_TPM_anno$log2FC >= upper | df_sleuth_across_with_TPM_anno$log2FC <= lower)],aes(ymin = log2FC-se_log2FC, ymax = log2FC+se_log2FC,color = category), na.rm = F, size = 0.35, fatten = 0.75) + 
  geom_pointrange(data = df_sleuth_across_with_TPM_anno[(df_sleuth_across_with_TPM_anno$log2FC <= upper | df_sleuth_across_with_TPM_anno$log2FC >= lower)], aes(ymin = log2FC-se_log2FC, ymax = log2FC+se_log2FC,color = category), na.rm = F, size = 0.35, fatten = 0.75, alpha = .25) +
  theme_AL_box_rotX() + 
  theme(axis.text.x = element_text(size = 4), legend.position = "none") + 
  coord_cartesian(ylim=c(-1,1)) + 
  labs(x="", y = "Mean FC (log2)") + 
  scale_color_manual(values = c("PAR" = "#00A087FF", "Inactive" = "#869a9a", "Variable" = "purple","Escape" = "red", "Potential" = "blue")) + 
  geom_hline(yintercept = median, linetype = "dashed", col = "red") + 
  geom_hline(yintercept = upper, col = "red") + 
  geom_hline(yintercept = lower, col = "red")+
  geom_text_repel(data = subset(df_sleuth_across_with_TPM_anno, gene %in% c("XIST", "KDM6A", "DDX3X", "FOXP3", "CD99", "TLR7", "SEPT6", "PUDP", "XG", "P2RY8")), aes(label = gene), direction = "both", min.segment.length = 0.001, nudge_y = 0.1, nudge_x = -1.5, size = 2) +
  geom_text_repel(data = subset(df_sleuth_across_with_TPM_anno, category %in% c("Inactive", "Potential") & df_sleuth_across_with_TPM_anno$pval < 0.00001 & df_sleuth_across_with_TPM_anno$log2FC < lower), aes(label = gene), direction = "both", min.segment.length = 0.001, nudge_y = 0.1, nudge_x = -1.5, size = 2) +
  geom_text_repel(data = subset(df_sleuth_across_with_TPM_anno, category %in% c("Escape", "Variable") & df_sleuth_across_with_TPM_anno$pval < 0.00001 & df_sleuth_across_with_TPM_anno$log2FC < lower), aes(label = gene), direction = "both", min.segment.length = 0.001, nudge_y = 0.1, nudge_x = -1.5, size = 2) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text.x = element_text(margin = margin(.01, 0, .01, 0, "cm")))+
  scale_y_continuous(breaks=c(1,0,-1))


ggsave2(filename = "04_plots/Figure_1e.pdf", width = 250, 
        height = 100, units = "mm", 
        plot_grid(fig1e)
)

# plot
fig1f <- 
  ggplot(data = df_sleuth_across_with_TPM_anno[df_sleuth_across_with_TPM_anno$meanz > 1], aes(x = category, y = log2FC)) + 
  geom_quasirandom(size = 0.75, aes(col = category), alpha = .5) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
  geom_hline(yintercept = 0, color = "red", alpha = 0.75) + 
  theme_AL_simple() + 
  scale_color_manual(values = c("PAR" = "#00A087FF", "Inactive" = "#869a9a", "Variable" = "purple", "Escape" = "red", "Potential" = "blue")) + 
  theme(legend.title = element_blank(), legend.position = "none", axis.text.x = element_text(size=4)) + 
  labs (x = "", y = "FC (log2)") + 
  scale_y_continuous(breaks=c(-1.0,0,1)) + 
  coord_cartesian(ylim=c(-1.1,1.8))+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

# get comparisons for statistical testing
my_comparisons <- list( c("PAR", "Escape"), c("PAR", "Variable"), c("PAR", "Inactive"), c("PAR", "Potential")  )

# test
stat.test <- compare_means(ref.group = "PAR", log2FC ~ category,  data = df_sleuth_across_with_TPM_anno[df_sleuth_across_with_TPM_anno$meanz > 1], method = "t.test",p.adjust.method = "BH")
stat.test$y.position <- c(1.1,1.25,1.4,1.55)

# add pvals
fig1f_with_pval <- fig1f + stat_pvalue_manual(stat.test, label = "p.adj", size = 3.25)


ggsave2(filename = "04_plots/Figure_1f.pdf", width = 75, 
        height = 100, units = "mm", 
        plot_grid(fig1f_with_pval)
)





# cap log2FC to 0.51 for nice plots.
df_sleuth_split_with_TPM_anno_capped <- df_sleuth_split_with_TPM_anno[df_sleuth_split_with_TPM_anno$category == "PAR" & df_sleuth_split_with_TPM_anno$meanz > 1]

df_sleuth_split_with_TPM_anno_capped[df_sleuth_split_with_TPM_anno_capped$log2FC > 0.5]$log2FC <- 0.51
df_sleuth_split_with_TPM_anno_capped[df_sleuth_split_with_TPM_anno_capped$log2FC < -0.5]$log2FC <- -0.51

# plot

fig2a <- 
  ggplot(data = df_sleuth_split_with_TPM_anno_capped, aes(x = factor(gene, levels = chrx_gene_order), y = log2FC)) + 
  geom_point(data=df_sleuth_split_with_TPM_anno_capped[df_sleuth_split_with_TPM_anno_capped$log2FC != 0.51 & df_sleuth_split_with_TPM_anno_capped$log2FC != -0.51,], size = 0.75) +
  geom_point(data=df_sleuth_split_with_TPM_anno_capped[df_sleuth_split_with_TPM_anno_capped$log2FC == 0.51 | df_sleuth_split_with_TPM_anno_capped$log2FC == -0.51,], size = 0.75, col = "red") +
  theme_AL_box_rotX() + 
  geom_hline(data = PAR_lines, aes(yintercept = median), alpha = 0.5, col = "red", linetype = "dashed") +
  geom_hline(data = PAR_lines, aes(yintercept = upper), alpha = 0.5, col = "red", linetype = "dashed") +
  geom_hline(data = PAR_lines, aes(yintercept = lower), alpha = 0.5, col = "red", linetype = "dashed") +
  facet_wrap(~factor(celltype, levels = c("ETP", "DN", "DPearly" ,"DPlate", "CD4SP", "CD8SP")), ncol = 6) + 
  labs(x = "", y = "FC (log2)") + 
  theme(axis.text.x = element_text(size = 6, face = "italic")) + 
  coord_cartesian(ylim = c(-0.5,0.5)) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + 
  theme(strip.text.x = element_text(margin = margin(.01, 0, .01, 0, "cm")))+
  scale_y_continuous(breaks=c(0.5,0,-0.5)) +
  geom_hline(yintercept = 0, lty = 1)



ggsave2(filename = "04_plots/Figure_2a.pdf", width = 250, 
        height = 65, units = "mm", 
        plot_grid(fig2a)
)


# Export source data
write.table(df_sleuth_across_with_TPM_anno[df_sleuth_across_with_TPM_anno$meanz > 1], file = "06_source_data/Fig_1e.tsv", quote = F, row.names = F, sep = "\t")
write.table(df_sleuth_across_with_TPM_anno[df_sleuth_across_with_TPM_anno$meanz > 1], file = "06_source_data/Fig_1f.tsv", quote = F, row.names = F, sep = "\t")
write.table(df_sleuth_split_with_TPM_anno_capped, file = "06_source_data/Fig_2a.tsv", quote = F, row.names = F, sep = "\t")
