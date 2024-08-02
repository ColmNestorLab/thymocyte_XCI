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

source("03_r_scripts/plot_parameters.R")

## Read in sleuth output tables
df_sleuth_across <- fread("Supplementary Table 5")
df_sleuth_split <- fread("Supplementary Table 6")

# read in TPM data
raw_tpm <- fread("Supplementary Table 3")

#melt
melt_raw_tpm <- melt(raw_tpm)

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

melt_raw_tpm <- melt_raw_tpm[melt_raw_tpm$sex != "KEK",]

# summarize, get means and keep only genes in our sleuth df
melt_raw_tpm_across <- melt_raw_tpm %>% dplyr::group_by(gene) %>% dplyr::summarise(meanz = mean(value))
melt_raw_tpm_across <- melt_raw_tpm_across[melt_raw_tpm_across$gene %in% df_sleuth_across$target_id,]

melt_raw_tpm_means_split <- melt_raw_tpm %>% dplyr::group_by(cell, gene) %>% dplyr::summarise(meanz = mean(value))
melt_raw_tpm_means_split <- melt_raw_tpm_means_split[melt_raw_tpm_means_split$gene %in% df_sleuth_across$target_id,]

# merge means with sleuth dfs, keeping only chrX genes
df_sleuth_across_with_TPM <- subset(merge(df_sleuth_across, melt_raw_tpm_across, by.x = c("target_id"), by.y = c("gene"), all.x =T), seqnames == "chrX")
df_sleuth_split_with_TPM <- subset(merge(df_sleuth_split, melt_raw_tpm_means_split, by.x = c("target_id", "Celltype"), by.y = c("gene", "cell"), all.x =T), seqnames == "chrX")

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
df_sleuth_across_with_TPM_anno <- merge(df_sleuth_across_with_TPM, chrx_annotation, by.x = c("target_id"), by.y = c("Gene_name"), all.x = T)
df_sleuth_split_with_TPM_anno <- merge(df_sleuth_split_with_TPM, chrx_annotation, by.x = c("target_id"), by.y = c("Gene_name"), all.x = T)

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
extreme_vals_across <- df_sleuth_across_with_TPM_anno[df_sleuth_across_with_TPM_anno$b > 1 | df_sleuth_across_with_TPM_anno$b < -1]
extreme_vals_split <- df_sleuth_split_with_TPM_anno[df_sleuth_split_with_TPM_anno$b > 1 | df_sleuth_split_with_TPM_anno$b < -1]
df_sleuth_across_with_TPM_anno[df_sleuth_across_with_TPM_anno$b > 1]$b <- 1
df_sleuth_across_with_TPM_anno[df_sleuth_across_with_TPM_anno$b < -1]$b <- -1
df_sleuth_split_with_TPM_anno[df_sleuth_split_with_TPM_anno$b > 1]$b <- 1
df_sleuth_split_with_TPM_anno[df_sleuth_split_with_TPM_anno$b < -1]$b <- -1

# calc stuff plot adding to the plot for selecting, post-squish.
statz_category_squish_across <- df_sleuth_across_with_TPM_anno %>%   group_by(category) %>%  get_summary_stats(b, type = "common")
setDT(statz_category_squish_across)

statz_category_squish_split <- df_sleuth_split_with_TPM_anno %>%   group_by(category, Celltype) %>%  get_summary_stats(b, type = "common")
setDT(statz_category_squish_split)

median <- statz_category_squish_across[statz_category_squish_across$category == "Inactive"]$median
upper <- statz_category_squish_across[statz_category_squish_across$category == "Inactive"]$median + statz_category_squish_across[statz_category_squish_across$category == "Inactive"]$iqr
lower <- statz_category_squish_across[statz_category_squish_across$category == "Inactive"]$median - statz_category_squish_across[statz_category_squish_across$category == "Inactive"]$iqr

cells_PAR <- statz_category_squish_split[statz_category_squish_split$category == "PAR"]$Celltype
median_PAR <- statz_category_squish_split[statz_category_squish_split$category == "PAR"]$median
upper_PAR <- statz_category_squish_split[statz_category_squish_split$category == "PAR"]$median + statz_category_squish_split[statz_category_squish_split$category == "PAR"]$iqr
lower_PAR <- statz_category_squish_split[statz_category_squish_split$category == "PAR"]$median - statz_category_squish_split[statz_category_squish_split$category == "PAR"]$iqr

PAR_lines <- data.frame(Celltype = cells_PAR, median = median_PAR, upper = upper_PAR, lower = lower_PAR)

chrx_gene_order <- unique(subset(fread("97_indexes_annotations/GCF_000001405.38_GRCh38.p12_genomic_with_correct_contig_names.tsv"), seqnames == "chrX")$gene_id)





################ Figure 1E #################

fig1e <- 
  ggplot(data = df_sleuth_across_with_TPM_anno[df_sleuth_across_with_TPM_anno$meanz > 1], aes(x = target_id, y = b)) + 
  facet_grid(~category, scales = "free_x", space = "free") +  
  geom_pointrange(data = df_sleuth_across_with_TPM_anno[(df_sleuth_across_with_TPM_anno$b >= upper | df_sleuth_across_with_TPM_anno$b <= lower)],aes(ymin = b-se_b, ymax = b+se_b,color = category), na.rm = F, size = 0.35, fatten = 0.75) + 
  geom_pointrange(data = df_sleuth_across_with_TPM_anno[(df_sleuth_across_with_TPM_anno$b <= upper | df_sleuth_across_with_TPM_anno$b >= lower)], aes(ymin = b-se_b, ymax = b+se_b,color = category), na.rm = F, size = 0.35, fatten = 0.75, alpha = .25) +
  theme_AL_box_rotX() + 
  theme(axis.text.x = element_text(size = 4), legend.position = "none") + 
  coord_cartesian(ylim=c(-1,1)) + 
  labs(x="", y = "Mean FC (log2)") + 
  scale_color_manual(values = c("PAR" = "#00A087FF", "Inactive" = "#869a9a", "Variable" = "purple","Escape" = "red", "Potential" = "blue")) + 
  geom_hline(yintercept = median, linetype = "dashed", col = "red") + 
  geom_hline(yintercept = upper, col = "red") + 
  geom_hline(yintercept = lower, col = "red")+
  geom_text_repel(data = subset(df_sleuth_across_with_TPM_anno, target_id %in% c("XIST", "KDM6A", "DDX3X", "FOXP3", "CD99", "TLR7", "SEPT6", "PUDP", "XG", "P2RY8")), aes(label = target_id), direction = "both", min.segment.length = 0.001, nudge_y = 0.1, nudge_x = -1.5, size = 2) +
  geom_text_repel(data = subset(df_sleuth_across_with_TPM_anno, category %in% c("Inactive", "Potential") & df_sleuth_across_with_TPM_anno$pval < 0.00001 & df_sleuth_across_with_TPM_anno$b < lower), aes(label = target_id), direction = "both", min.segment.length = 0.001, nudge_y = 0.1, nudge_x = -1.5, size = 2) +
  geom_text_repel(data = subset(df_sleuth_across_with_TPM_anno, category %in% c("Escape", "Variable") & df_sleuth_across_with_TPM_anno$pval < 0.00001 & df_sleuth_across_with_TPM_anno$b < lower), aes(label = target_id), direction = "both", min.segment.length = 0.001, nudge_y = 0.1, nudge_x = -1.5, size = 2) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text.x = element_text(margin = margin(.01, 0, .01, 0, "cm")))+
  scale_y_continuous(breaks=c(1,0,-1))


ggsave2(filename = "04_plots/Figure_1E.pdf", width = 250, 
        height = 100, units = "mm", 
        plot_grid(fig1e)
)

# get gene names based on the longplot
# first find escapees that does not show sex-boas
find_genez <- df_sleuth_across_with_TPM_anno[df_sleuth_across_with_TPM_anno$category == "Escape" & (b <= 0.20 & b > -0.20) & pval > 0.00001,]

#then, find genes that undergo xci (i.e. inactive genes) that have a male-bias
find_genez2 <- df_sleuth_across_with_TPM_anno[df_sleuth_across_with_TPM_anno$category == "Inactive" & (b < -0.2 & pval < 0.00001),]

################ Figure 1F ## ###############


fig1f <- 
  ggplot(data = df_sleuth_across_with_TPM_anno[df_sleuth_across_with_TPM_anno$meanz > 1], aes(x = category, y = b)) + 
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

my_comparisons <- list( c("PAR", "Escape"), c("PAR", "Variable"), c("PAR", "Inactive"), c("PAR", "Potential")  )

stat.test <- compare_means(ref.group = "PAR", b ~ category,  data = df_sleuth_across_with_TPM_anno[df_sleuth_across_with_TPM_anno$meanz > 1], method = "t.test",p.adjust.method = "BH")
stat.test$y.position <- c(1.1,1.25,1.4,1.55)

fig1f_with_pval <- fig1f + stat_pvalue_manual(stat.test, label = "p.adj", size = 3.25)


ggsave2(filename = "04_plots/Figure_1F.pdf", width = 75, 
        height = 100, units = "mm", 
        plot_grid(fig1f_with_pval)
)



################ Suppl. Figure 1G #################
df_sleuth_split_with_TPM_anno_uncapped <- df_sleuth_split_with_TPM_anno[df_sleuth_split_with_TPM_anno$category == "PAR" & df_sleuth_split_with_TPM_anno$meanz > 1]
df_sleuth_split_with_TPM_anno_capped <- df_sleuth_split_with_TPM_anno[df_sleuth_split_with_TPM_anno$category == "PAR" & df_sleuth_split_with_TPM_anno$meanz > 1]
extreme_values_split <- df_sleuth_split_with_TPM_anno[(df_sleuth_split_with_TPM_anno$b > 0.5 | df_sleuth_split_with_TPM_anno$b < -0.5) & df_sleuth_split_with_TPM_anno$category == "PAR" & df_sleuth_split_with_TPM_anno$meanz > 1]

df_sleuth_split_with_TPM_anno_capped[df_sleuth_split_with_TPM_anno_capped$b > 0.5]$b <- 0.51
df_sleuth_split_with_TPM_anno_capped[df_sleuth_split_with_TPM_anno_capped$b < -0.5]$b <- -0.51

suppl_fig2C <- 
  ggplot(data = df_sleuth_split_with_TPM_anno_capped, aes(x = factor(target_id, levels = chrx_gene_order), y = b)) + 
  geom_point(data=df_sleuth_split_with_TPM_anno_capped[df_sleuth_split_with_TPM_anno_capped$b != 0.51 & df_sleuth_split_with_TPM_anno_capped$b != -0.51,], size = 0.75) +
  geom_point(data=df_sleuth_split_with_TPM_anno_capped[df_sleuth_split_with_TPM_anno_capped$b == 0.51 | df_sleuth_split_with_TPM_anno_capped$b == -0.51,], size = 0.75, col = "red") +
  theme_AL_box_rotX() + 
  geom_hline(data = PAR_lines, aes(yintercept = median), alpha = 0.5, col = "red", linetype = "dashed") +
  geom_hline(data = PAR_lines, aes(yintercept = upper), alpha = 0.5, col = "red", linetype = "dashed") +
  geom_hline(data = PAR_lines, aes(yintercept = lower), alpha = 0.5, col = "red", linetype = "dashed") +
  facet_wrap(~factor(Celltype, levels = c("ETP", "DN", "DPearly" ,"DPlate", "CD4SP", "CD8SP")), ncol = 6) + 
  labs(x = "", y = "FC (log2)") + 
  theme(axis.text.x = element_text(size = 6, face = "italic")) + 
  coord_cartesian(ylim = c(-0.5,0.5)) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + 
  theme(strip.text.x = element_text(margin = margin(.01, 0, .01, 0, "cm")))+
  scale_y_continuous(breaks=c(0.5,0,-0.5)) +
  geom_hline(yintercept = 0, lty = 1)



ggsave2(filename = "04_plots/Figure_1G.pdf", width = 250, 
        height = 65, units = "mm", 
        plot_grid(suppl_fig2C)
)
