library(data.table)
library(cowplot)
library(ggplot2)
library(ggfortify)
library(cowplot)
library(ggrepel)
library(tidyverse)
library(ggbeeswarm)
library(Seurat)
source("plot_parameters.R")


### PLOT PCA ###
# read in data
ss2_df <- data.table(read_excel("02_tidy_data/Suppl_tables/Supplementary Tables.xlsx", sheet = 11, skip = 1, col_types = c("text", "numeric", "text", "numeric", "numeric", "numeric", "text", "text", "numeric", "numeric", "numeric", "numeric", "numeric", "logical")))[contig == "chrX" & pass_all_filters == T & !is.na(gene)]

# make gene separate on position, as to keep them separate.
ss2_df$gene_pos <- paste0(ss2_df$gene, "_", ss2_df$position)

# grep XIST hetSNPs.
XIST_test <- ss2_df[grepl(ss2_df$gene, pattern = "XIST"),]

# make plot order
XIST_test$celltype <- factor(XIST_test$celltype, levels = cell_order)
XIST_test_order <- unique(XIST_test[order(XIST_test$Reference_ratio, -XIST_test$celltype),]$cell_id)

# plot
Suppl_Fig_6a <- 
  ggplot(XIST_test, aes(x=factor(cell_id, levels = rev(XIST_test_order)), y= Reference_ratio)) + 
  geom_point(data=XIST_test[XIST_test$Reference_ratio > 0.5,], shape = 24, fill = "red", col="red") + 
  geom_point(data=XIST_test[XIST_test$Reference_ratio < 0.5,], shape = 21, fill = "blue", col="blue") + 
  geom_hline(yintercept = c(0.01,0.99)) + 
  theme_AL_box_rotX() + 
  theme(legend.position = "top", axis.text.x = element_text(size=4)) + 
  labs(x="")

# Export source data
write.table(XIST_test, file = "06_source_data/Suppl_Fig_6a.tsv", quote = F, row.names = F, sep = "\t")


# Assign cluster, using XIST ASE.
XIST_test$XIST_allele <- "unclear"
XIST_test[XIST_test$Reference_ratio < 0.01,]$XIST_allele <- "cluster1"
XIST_test[XIST_test$Reference_ratio > 0.99,]$XIST_allele <- "cluster2"

# Exclude cells with an unclear XIST ASE (doublets?)
excluded_cells_bc_of_XIST <- XIST_test[XIST_test$XIST_allele == "unclear",]$cell_id

# Also exlude cells which, in the PCA, cluster between cluster1 and cluster2.
excluded_cells_bc_of_clustering <- c("ETP_38", "ETP_71", "ETP_89", "CD8SP_96")

# Also excluded one cell (DPearly_92) in which we only have coverage for ZFX, which makes classification impossible due to it being an escape gene.
# ETP_31, CD8SP_78, CD8SP_38, CD8SP_81, CD8SP_61 and DPlate_49 was excluded bc of disconcordant XIST and inactive-X ASE, most likely due to doublets.
excluded_cells_bc_of_escape_gene_only_or_disconcordant_AE <- c("DPearly_92", c("ETP_31", "CD8SP_78", "CD8SP_38", "CD8SP_81", "CD8SP_61", "DPlate_49"))

# all excluded cells
excluded_cells <- c(excluded_cells_bc_of_clustering, excluded_cells_bc_of_XIST, excluded_cells_bc_of_escape_gene_only_or_disconcordant_AE)

# merge
ss2_df_smash <- merge(ss2_df, XIST_test[,c("XIST_allele", "cell_id")], by = c("cell_id"), all.x =F)

# Exclude outliers
ss2_df_smash_filt <- ss2_df_smash[pass_all_filters == T & !ss2_df_smash$cell_id %in% excluded_cells,]

# make matrix
mat_ss2_df <- dcast(unique(ss2_df_smash_filt[ss2_df_smash_filt$pass_all_filters == T & !ss2_df_smash_filt$cell_id %in% excluded_cells,]), gene_pos~cell_id, value.var = "Reference_ratio")
mat_ss2_df[is.na(mat_ss2_df)] <- 0

# add rownames, remove gene_pos column
rownames(mat_ss2_df) <- mat_ss2_df$gene_pos
mat_ss2_df$gene_pos <- NULL

# transform
t_mat_ss2_df <- base::t(mat_ss2_df)
colnames(t_mat_ss2_df) <- rownames(mat_ss2_df)

# make metaz df for PCA plotting
metaz <- data.frame(unique(ss2_df_smash_filt[ss2_df_smash_filt$cell_id %in% colnames(mat_ss2_df) & pass_all_filters == T,c("cell_id", "celltype", "XIST_allele")]))

# add rownames
rownames(metaz) <- metaz$cell_id

# run prcomp
pca_thymo <- prcomp(t_mat_ss2_df)

# plot
Suppl_Fig_6b <- autoplot(pca_thymo,size=1.5,stroke=NA, x = 1, y=2, data=metaz, colour="celltype", label = F, shape = "XIST_allele") + theme_AL_box() + theme(legend.position = "top") + scale_color_manual(values = colors)

# Export source data
write.table(Suppl_Fig_6b$data[,1:5], file = "06_source_data/Suppl_Fig_6b.tsv", quote = F, row.names = F, sep = "\t")


# cells clustered by XIST ASE.
cells_classified_with_XIST <- Suppl_Fig_6b$data$cell_id

# Read in classification/clustering (file is created in 'cluster_cells_generate_Supplementary_Table_10.R')
final_classification <- data.table(read_excel("02_tidy_data/Suppl_tables/Supplementary Tables.xlsx", sheet = 12, skip = 1, col_types = c("text", "text", "text")))

# Merge the classification with the ss2 dataframe
ss2_df_fin <- merge(ss2_df[ss2_df$pass_all_filters == T,], final_classification, by = "cell_id")

# keep genes we want to use for clustering
genes_of_interest <- c("ITM2A_79360867", "ITM2A_79360874", "MORF4L2_103676032", "BEX4_103216603", "PIN4_72181757", "PIN4_72181764", "PIN4_72196871", "TMSB4X_12976131", "TMSB4X_12975167", "LAMP2_120437686", "LAMP2_120456678", "ATRX_77682471")


# Use the cells we have already clustered (using XIST allelic expression) and tease out other usable genes.
# get full data frame, adding the XIST stuff
ss2_df_with_XIST_allele_unfilt <- merge(ss2_df[!ss2_df$cell_id %in% excluded_cells,], XIST_test[,c("XIST_allele", "cell_id")], by = c("cell_id"), all.x =T)

# add tag where we state whether XIST exclusively was used for clustering or if we used other genes.
ss2_df_with_XIST_allele_unfilt$manual_class <- ifelse(ss2_df_with_XIST_allele_unfilt$cell_id %in% cells_classified_with_XIST, yes = "XIST_clustering", no = "manual_clustering")

# add tag to see whether or not we have already assiged the cell to a cluster.
ss2_df_with_XIST_allele_unfilt$DONE <- ifelse(ss2_df_with_XIST_allele_unfilt$cell_id %in% cells_classified_with_XIST, yes = "done", no = "not_done")


#plot
Suppl_Fig_6c <- 
  plot_grid(ggplot(ss2_df_with_XIST_allele_unfilt[ss2_df_with_XIST_allele_unfilt$gene_pos %in% genes_of_interest,], 
                   aes(x=XIST_allele, y=Reference_ratio, col = factor(position))) + 
              geom_quasirandom(dodge.width = 0.5) + 
              geom_hline(yintercept = c(0.05,0.95), lty=2) + 
              theme_AL_box_rotX() +
              theme(legend.position = "top")+
              facet_wrap(~gene, nrow = 1))

# export source data
write.table(ss2_df_with_XIST_allele_unfilt[ss2_df_with_XIST_allele_unfilt$gene_pos %in% genes_of_interest,], file = "06_source_data/Suppl_Fig_6c.tsv", quote = F, row.names = F, sep = "\t")



# Plot all genes per cluster. Mean function.
Mean <-  function(x) {
  setNames(rep(mean(x), 5), c("ymin", "lower", "middle", "upper", "ymax"))
}

# make plot df
final_df <- ss2_df_fin

# count how many cells we detect each gene in
final_df_count <- final_df %>% dplyr::group_by(gene_pos, final_clustering) %>% dplyr::count()

# add that count to the data frame
final_df_add_n <- merge(final_df, final_df_count, by= c("gene_pos", "final_clustering"))

# make a new variable, highlighting how many cells we cover each gene in.
final_df_add_n$gene_pos_n <- paste0(final_df_add_n$gene_pos, " (n = ", final_df_add_n$n, ")")

# get summary statistics for ordering of the X-axis.
final_df_order_cluster1 <- final_df_add_n[final_df_add_n$final_clustering == "cluster1",] %>% dplyr::group_by(gene_pos_n) %>% dplyr::summarise(medianz=mean(Reference_ratio))
final_df_order_cluster1 <- final_df_order_cluster1[order(final_df_order_cluster1$medianz, decreasing = T),]
final_df_order_cluster1$medianz <- round(final_df_order_cluster1$medianz, 4)
final_df_order_cluster2 <- final_df_add_n[final_df_add_n$final_clustering == "cluster2",] %>% dplyr::group_by(gene_pos_n) %>% dplyr::summarise(medianz=mean(Reference_ratio))
final_df_order_cluster2 <- final_df_order_cluster2[order(final_df_order_cluster2$medianz, decreasing = T),]
final_df_order_cluster2$medianz <- round(final_df_order_cluster2$medianz, 4)


# Plot all X-linked genes, removing genes which we only detect in one cell.
sfig6d <-  plot_grid(ncol=1, 
                            ggplot(final_df_add_n[final_clustering == "cluster1" & n > 1], 
                                   aes(x=factor(gene_pos_n, levels = final_df_order_cluster1$gene_pos_n), y=Reference_ratio)) + 
                              geom_quasirandom(aes(col=final_clustering), show.legend = F, size = 1.5) + 
                              theme_AL_box_rotX() + 
                              theme(axis.text.x = element_text(size=6))+
                              labs(title = "cluster1", x="", y="Reference ratio (refCount / totalCount)") + 
                              geom_hline(yintercept = 0.5, lty=2)  + 
                              stat_summary(fun.data = Mean, geom = "boxplot") + 
                              #geom_text(data=final_df_count[final_df_count$class == "cluster_1",], aes(x=gene_pos_n, y=1, label = n), size = 2)+
                              geom_hline(yintercept = c(0.1,0.9))
                            
                            ,
                            
                            ggplot(final_df_add_n[final_clustering == "cluster2" & n > 1], 
                                   aes(x=factor(gene_pos_n, levels = final_df_order_cluster2$gene_pos_n), y=Reference_ratio)) + 
                              geom_quasirandom(aes(col=final_clustering), show.legend = F, size = 1.5) + 
                              theme_AL_box_rotX() + 
                              theme(axis.text.x = element_text(size=6))+
                              labs(title = "cluster_2", x="", y="Reference ratio (refCount / totalCount)") + 
                              geom_hline(yintercept = 0.5, lty=2)  + 
                              stat_summary(fun.data = Mean, geom = "boxplot") + 
                              #geom_text(data=final_df_count[final_df_count$class == "cluster_2",], aes(x=gene_pos_n, y=1, label = n), size = 2)+
                              scale_color_manual(values = c("cluster2" = "#3A9BDC"))+
                              geom_hline(yintercept = c(0.1,0.9)))
      


# Export source data
write.table(final_df_add_n[n > 1], file = "06_source_data/Suppl_Fig_6d.tsv", quote = F, row.names = F, sep = "\t")

# Export source data for fig. 5c.
cluster1_non_zero_or_one <- final_df_order_cluster1[final_df_order_cluster1$medianz < 0.9918 & final_df_order_cluster1$medianz > 0.0039,]
cluster2_non_zero_or_one <- final_df_order_cluster2[final_df_order_cluster2$medianz < 0.9999 & final_df_order_cluster2$medianz > 0.0002,]
write.table(final_df_add_n[final_df_add_n$n > 1 & final_df_add_n$gene_pos_n %in% c(cluster1_non_zero_or_one$gene_pos_n, cluster2_non_zero_or_one$gene_pos_n),], file = "06_source_data/Fig_5c.tsv", quote = F, row.names = F, sep = "\t")


# plot the minor allele read counts for inactive genes (essentially looking for any signs of 'leaky' Xi expression or reactivation of Xi).
final_df_2 <- final_df

# Get which allele is the minor one
final_df_2$inactive_readcount <- pmin(final_df_2$altCount, final_df_2$refCount)
final_df_2$active_readcount <- pmax(final_df_2$altCount, final_df_2$refCount)

# list of escapees to not include.
escapees <- c("RPS4X", "CD99", "TMSB4X", "SLC25A6", "ASMTL", "ZFX", "ZRSR2", "PRKX", "GTPBP6", "XG")

# summarize the number of reads coming from the Xi per cell and cluster.
final_df_filt_stats <- final_df_2[!final_df_2$gene %in% escapees,]  %>% dplyr::group_by(cell_id, celltype, final_clustering) %>% dplyr::summarize(sum_inactive_readcount = sum(inactive_readcount),sum_active_readcount = sum(active_readcount))

# calculate how big a fraction the minor allele reads make up the total read count.
final_df_filt_stats$fraction_minor_allele_read_count <- (final_df_filt_stats$sum_inactive_readcount / (final_df_filt_stats$sum_inactive_readcount + final_df_filt_stats$sum_active_readcount)) *100

# get the number of reads from the Xi in inactive genes per cell 
stats_df <- final_df_filt_stats %>% dplyr::group_by(sum_inactive_readcount, celltype, final_clustering) %>% dplyr::count()

# Group the number of inactive gene Xi reads per cell into 0, 1 or more than 2 reads.
stats_df$sum_inactive_readcount <- as.numeric(stats_df$sum_inactive_readcount)
stats_df$categoriezzz <- "LOL"
stats_df[stats_df$sum_inactive_readcount == 0,]$categoriezzz <- "0"
stats_df[stats_df$sum_inactive_readcount == 1,]$categoriezzz <- "1"
stats_df[stats_df$sum_inactive_readcount >= 2,]$categoriezzz <- ">=2"

# order.
catz_order <- c("0", "1", ">=2")

# get the distribution per celltype
stats_df_test <- stats_df %>% dplyr::group_by(celltype, categoriezzz) %>% dplyr::summarise(sum_n_cats = sum(n))

# plot em
Fig_5b <- 
  ggplot(stats_df_test, aes(x=factor(celltype, cell_order), y=sum_n_cats, label = sum_n_cats, fill = factor(categoriezzz, catz_order))) + 
  geom_bar(stat="identity") + 
  geom_text(position = "stack") + 
  theme_AL_box_rotX(legend.title=element_blank(), legend.position = "top")



# Export source data
write.table(stats_df_test, file = "06_source_data/Fig_5b.tsv", quote = F, row.names = F, sep = "\t")



# read in the seurat object
seurat_object <- LoadSeuratRds("02_tidy_data/thymocyte/SS2/thymocytes_SS2.Rds")

# plot
fig_5a <- DimPlot(seurat_object, label = T) + scale_color_manual(values = colors)

### Export plots!
ggsave2("04_plots/Fig_5a_5b_5c_Suppl_fig_6a_6b_6c_6d.pdf", height=16, width = 14,
        plot_grid(
          plot_grid(ncol=1, rel_heights = c(0.5,0.4,1),
                    plot_grid(ncol=4, rel_widths = c(1,0.25,0.5,0.6),
                              fig_5a,
                              Fig_5b,
                              Suppl_Fig_6a, 
                              Suppl_Fig_6b
                    ),
                    Suppl_Fig_6c,
                    sfig6d
          )
        )
)

# export source data
umap_coords <- Embeddings(seurat_object, "umap")
fig_5a <- as.data.frame(umap_coords)
fig_5a$cellid <- rownames(fig_5a)
write.table(fig_5a, file = "06_source_data/Fig_5a.tsv", quote = F, row.names = F, sep = "\t")



