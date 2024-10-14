library(data.table)
library(ggplot2)
library(cowplot)
source("plot_parameters.R")

#Load sleuth object
source("03_r_scripts/load_scripts/Load_thymocyte_sleuth_expression_data.R")

## get PCA PCs and export them. I will run them through GSEA preranked.
pca_thymo <- prcomp(t(txi_thymo_transition$counts))

write.table(sort(pca_thymo$rotation[abs(pca_thymo$rotation[,1]) > 0.01,1],decreasing = T),"02_tidy_data/GSEA/software_inputs/PCA_fig1b_PC1_and_PC2_preranked/thymocyte_PCA_component1.rnk",quote = F,sep = "\t",col.names = F)
write.table(sort(pca_thymo$rotation[abs(pca_thymo$rotation[,2]) > 0.01,2],decreasing = T),"02_tidy_data/GSEA/software_inputs/PCA_fig1b_PC1_and_PC2_preranked/thymocyte_PCA_component2.rnk",quote = F,sep = "\t",col.names = F)

# Run the above outputs through GSEA software (Version 4.3.3, build: 16) #
# ftp.broadinstitute.org://pub/gsea/msigdb/human/gene_sets/c5.go.bp.v2023.2.Hs.symbols.gmt
# all default settings, except no_collapse

# plot suppl. figure 2b
# Read in PC1 GSEAPreranked output
PCA1_neg_GSEA <- fread("02_tidy_data/GSEA/software_outputs/PC1.GseaPreranked.1715841767284/gsea_report_for_na_neg_1715841767284.tsv")
PCA1_pos_GSEA <- fread("02_tidy_data/GSEA/software_outputs/PC1.GseaPreranked.1715841767284/gsea_report_for_na_pos_1715841767284.tsv")

# make names easier to read.
PCA1_neg_GSEA$NAME <- gsub(PCA1_neg_GSEA$NAME, pattern = "GOBP_", replacement = "")
PCA1_neg_GSEA$NAME <- gsub(PCA1_neg_GSEA$NAME, pattern = "_", replacement = " ")
PCA1_neg_GSEA$NAME <- gsub("(?<=\\b.)(.*?)\\b", "\\L\\1", PCA1_neg_GSEA$NAME, perl=TRUE)
PCA1_pos_GSEA$NAME <- gsub(PCA1_pos_GSEA$NAME, pattern = "GOBP_", replacement = "")
PCA1_pos_GSEA$NAME <- gsub(PCA1_pos_GSEA$NAME, pattern = "_", replacement = " ")
PCA1_pos_GSEA$NAME <- gsub("(?<=\\b.)(.*?)\\b", "\\L\\1", PCA1_pos_GSEA$NAME, perl=TRUE)

# Read in PC2 GSEAPreranked output
PCA2_neg_GSEA <- fread("02_tidy_data/GSEA/software_outputs/PC2.GseaPreranked.1715841790221/gsea_report_for_na_neg_1715841790221.tsv")
PCA2_pos_GSEA <- fread("02_tidy_data/GSEA/software_outputs/PC2.GseaPreranked.1715841790221/gsea_report_for_na_pos_1715841790221.tsv")

# make names easier to read.
PCA2_neg_GSEA$NAME <- gsub(PCA2_neg_GSEA$NAME, pattern = "GOBP_", replacement = "")
PCA2_neg_GSEA$NAME <- gsub(PCA2_neg_GSEA$NAME, pattern = "_", replacement = " ")
PCA2_neg_GSEA$NAME <- gsub("(?<=\\b.)(.*?)\\b", "\\L\\1", PCA2_neg_GSEA$NAME, perl=TRUE)
PCA2_pos_GSEA$NAME <- gsub(PCA2_pos_GSEA$NAME, pattern = "GOBP_", replacement = "")
PCA2_pos_GSEA$NAME <- gsub(PCA2_pos_GSEA$NAME, pattern = "_", replacement = " ")
PCA2_pos_GSEA$NAME <- gsub("(?<=\\b.)(.*?)\\b", "\\L\\1", PCA2_pos_GSEA$NAME, perl=TRUE)

# bind together
PCA1_components <- rbind(head(PCA1_neg_GSEA, 5), head(PCA1_pos_GSEA, 5))
PCA2_components <- rbind(head(PCA2_neg_GSEA, 5), head(PCA2_pos_GSEA, 5))

# set order
PCA1_order <- c(head(PCA1_neg_GSEA, 5)$NAME, head(PCA1_pos_GSEA, 5)$NAME)
PCA2_order <- c(head(PCA2_neg_GSEA, 5)$NAME, head(PCA2_pos_GSEA, 5)$NAME)

# plot
ggsave2(filename = "04_plots/Suppl_Figure_2b.pdf",
plot_grid(ncol=1,labels = c("PC1", "PC2"),
ggplot(PCA1_components, aes(x=NES,y=factor(NAME, levels = PCA1_order))) + geom_bar(stat="identity", width = 0.4) + coord_cartesian(xlim=c(-4,4)) + scale_x_continuous(position = "top", breaks = c(-4,-3,-2,-1,0,1,2,3,4)) + theme_AL_simple() + theme(axis.text.y = element_text(size=8)) + labs(y="", x="Normalised enrichment score (NES)") + geom_vline(xintercept = 0, lty=2),

plot_grid(ncol=2, NULL, rel_widths = c(0.125,1),
          ggplot(PCA2_components, aes(x=NES,y=factor(NAME, levels = PCA2_order))) + geom_bar(stat="identity", width = 0.4) + coord_cartesian(xlim=c(-3,3)) + scale_x_continuous(position = "top", breaks = c(-3,-2,-1,0,1,2,3)) + theme_AL_simple() + theme(axis.text.y = element_text(size=8)) + labs(y="", x="Normalised enrichment score (NES)") + geom_vline(xintercept = 0, lty=2))
))

# export source data
PCA1_components$tag <- "PCA1_components"
PCA2_components$tag <- "PCA2_components"
write.table(rbind(PCA1_components, PCA2_components), file = "06_source_data/Suppl_Fig_2b.tsv", quote = F, row.names = F, sep = "\t")


# plot suppl. figure 2d
# read in data
cluster1 <- head(dplyr::select(fread("02_tidy_data/GO_output/Clusters/Pantherdb.org_outputs/cluster1.txt"), 1, 6), 10)
cluster2 <- head(dplyr::select(fread("02_tidy_data/GO_output/Clusters/Pantherdb.org_outputs/cluster2.txt"), 1, 6), 10)
cluster3 <- head(dplyr::select(fread("02_tidy_data/GO_output/Clusters/Pantherdb.org_outputs/cluster3.txt"), 1, 6), 10)
cluster4 <- head(dplyr::select(fread("02_tidy_data/GO_output/Clusters/Pantherdb.org_outputs/cluster4.txt"), 1, 6), 10)
cluster5 <- head(dplyr::select(fread("02_tidy_data/GO_output/Clusters/Pantherdb.org_outputs/cluster5.txt"), 1, 6), 10)
cluster6 <- head(dplyr::select(fread("02_tidy_data/GO_output/Clusters/Pantherdb.org_outputs/cluster6.txt"), 1, 6), 10)
cluster7 <- head(dplyr::select(fread("02_tidy_data/GO_output/Clusters/Pantherdb.org_outputs/cluster7.txt"), 1, 6), 10)

# make names easier to read
cluster1$GO_short <- substring(cluster1$`GO biological process complete`, 1, nchar(cluster1$`GO biological process complete`)-13)
cluster2$GO_short <- substring(cluster2$`GO biological process complete`, 1, nchar(cluster2$`GO biological process complete`)-13)
cluster3$GO_short <- substring(cluster3$`GO biological process complete`, 1, nchar(cluster3$`GO biological process complete`)-13)
cluster4$GO_short <- substring(cluster4$`GO biological process complete`, 1, nchar(cluster4$`GO biological process complete`)-13)
cluster5$GO_short <- substring(cluster5$`GO biological process complete`, 1, nchar(cluster5$`GO biological process complete`)-13)
cluster6$GO_short <- substring(cluster6$`GO biological process complete`, 1, nchar(cluster6$`GO biological process complete`)-13)
cluster7$GO_short <- substring(cluster7$`GO biological process complete`, 1, nchar(cluster7$`GO biological process complete`)-13)

# change colnames
colnames(cluster1) <- c("GO_BP", "Fold_enrichment", "GO_BP_short")
colnames(cluster2) <- c("GO_BP", "Fold_enrichment", "GO_BP_short")
colnames(cluster3) <- c("GO_BP", "Fold_enrichment", "GO_BP_short")
colnames(cluster4) <- c("GO_BP", "Fold_enrichment", "GO_BP_short")
colnames(cluster5) <- c("GO_BP", "Fold_enrichment", "GO_BP_short")
colnames(cluster6) <- c("GO_BP", "Fold_enrichment", "GO_BP_short")
colnames(cluster7) <- c("GO_BP", "Fold_enrichment", "GO_BP_short")

# set fold enrichment as numeric, make sure factor is leveled, add cluster tag.
cluster1$Fold_enrichment <- as.numeric(cluster1$Fold_enrichment)
cluster1$GO_BP_short <- factor(cluster1$GO_BP_short, levels = rev(cluster1$GO_BP_short))
cluster1$cluster <- "Cluster 1"

cluster2$Fold_enrichment <- as.numeric(cluster2$Fold_enrichment)
cluster2$GO_BP_short <- factor(cluster2$GO_BP_short, levels = rev(cluster2$GO_BP_short))
cluster2$cluster <- "Cluster 2"

cluster3$Fold_enrichment <- as.numeric(cluster3$Fold_enrichment)
cluster3$GO_BP_short <- factor(cluster3$GO_BP_short, levels = rev(cluster3$GO_BP_short))
cluster3$cluster <- "Cluster 3"

cluster4$Fold_enrichment <- as.numeric(cluster4$Fold_enrichment)
cluster4$GO_BP_short <- factor(cluster4$GO_BP_short, levels = rev(cluster4$GO_BP_short))
cluster4$cluster <- "Cluster 4"

cluster5$Fold_enrichment <- as.numeric(cluster5$Fold_enrichment)
cluster5$GO_BP_short <- factor(cluster5$GO_BP_short, levels = rev(cluster5$GO_BP_short))
cluster5$cluster <- "Cluster 5"

cluster6$Fold_enrichment <- as.numeric(cluster6$Fold_enrichment)
cluster6$GO_BP_short <- factor(cluster6$GO_BP_short, levels = rev(cluster6$GO_BP_short))
cluster6$cluster <- "Cluster 6"

cluster7$Fold_enrichment <- as.numeric(cluster7$Fold_enrichment)
cluster7$GO_BP_short <- factor(cluster7$GO_BP_short, levels = rev(cluster7$GO_BP_short))
cluster7$cluster <- "Cluster 7"

# Rbind the cluster data
smash_cluster <- rbind(cluster1, cluster2, cluster3, cluster4, cluster5, cluster6, cluster7)

#plot
ggsave2(filename = "04_plots/Suppl_Figure_2d.pdf", height = 20,
ggplot(smash_cluster, aes(x=Fold_enrichment,y=GO_BP_short)) + geom_bar(stat="identity", width = 0.4) +scale_x_continuous(position = "top", breaks = c(0,2,4,6,8,10,12,14,16)) + theme_AL_simple() + theme(axis.text.y = element_text(size=8)) + labs(y="", x="fold enrichment") + geom_vline(xintercept = 1, lty=2, col = "red") + coord_cartesian(xlim=c(0,16)) + facet_wrap(~cluster, ncol = 1, scale = "free")
)


write.table(smash_cluster, file = "06_source_data/Suppl_Fig_2d.tsv", quote = F, row.names = F, sep = "\t")
