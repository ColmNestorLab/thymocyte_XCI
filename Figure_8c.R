options(stringsAsFactors = F)

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

######## read in data #############
message("loading sleuth object")
source("Load_thymocyte_sleuth_expression_data_with_Turner.R")

# read in sleuth results from normal thymocyte development analysis.
res_lrt <- data.table(read_excel("Supplementary Data.xlsx", sheet = 3, skip = 1, col_types = c("text", "numeric", "numeric", "numeric", "numeric")))

# Filter out genes not expressed in any subtype
tpm_filter <- txi_thymo$abundance[apply(txi_thymo$abundance,1,function(x) any(tapply(x,meta_thymo_transition$Celltype,function(y) mean(y>=1) )==1) ),]


# average per sex (i.e. ETP F1, F2, F3 into F_ETP and so on)
Femme <- c("F1", "F2", "F3", "F4")
Homme <- c("M1", "M3", "M4")

# make data frame
df_tpm_filter <- data.frame(tpm_filter)

# find the Turner data
turner_columns <- dplyr::select(df_tpm_filter, T1_ETP_RNA_S1, T1_DN_RNA_S2, T1_DPearly_RNA_S3, T1_DPlate_RNA_S4, T1_CD4SP_RNA_S5, T1_CD8SP_RNA_S6)
colnames(turner_columns) <- c("ETP_Turner", "DN_Turner", "DPearly_Turner", "DPlate_Turner", "CD4SP_Turner", "CD8SP_Turner")

# get female sample, calculate rowmeans per celltype
ETP_F <- dplyr::select(df_tpm_filter, F1_ETP_RNA_S5, F2_ETP_RNA_S1, F3_ETP_RNA_S7, F4_ETP_RNA_S19)
DN_F <- dplyr::select(df_tpm_filter, F1_DN_RNA_S6, F2_DN_RNA_S2, F3_DN_RNA_S8, F4_DN_RNA_S20)
DPearly_F <- dplyr::select(df_tpm_filter, F1_DPearly_RNA_S7, F2_DPearly_RNA_S3, F3_DPearly_RNA_S9, F4_DPearly_RNA_S21)
DPlate_F <- dplyr::select(df_tpm_filter, F1_DPlate_RNA_S8, F2_DPlate_RNA_S4, F3_DPlate_RNA_S10, F4_DPlate_RNA_S22)
CD4SP_F <- dplyr::select(df_tpm_filter, F1_CD4SP_RNA_S9, F2_CD4SP_RNA_S5, F3_CD4SP_RNA_S11, F4_CD4SP_RNA_S23)
CD8SP_F <- dplyr::select(df_tpm_filter, F1_CD8SP_RNA_S10, F2_CD8SP_RNA_S6, F3_CD8SP_RNA_S12, F4_CD8SP_RNA_S24)
ETP_F$meanz <- rowMeans(ETP_F)
DN_F$meanz <- rowMeans(DN_F)
DPearly_F$meanz <- rowMeans(DPearly_F)
DPlate_F$meanz <- rowMeans(DPlate_F)
CD4SP_F$meanz <- rowMeans(CD4SP_F)
CD8SP_F$meanz <- rowMeans(CD8SP_F)

# merge and change colnames
femme_means <- dplyr::select(cbind(ETP_F, DN_F$meanz, DPearly_F$meanz, DPlate_F$meanz, CD4SP_F$meanz, CD8SP_F$meanz), -1, -2, -3,-4)
colnames(femme_means) <- c("ETP_F", "DN_F", "DPearly_F", "DPlate_F", "CD4SP_F", "CD8SP_F")

# get male sample, calculate rowmeans per celltype
ETP_M <- dplyr::select(df_tpm_filter, M3_ETP_RNA_S7, M4_ETP_RNA_S13)
DN_M <- dplyr::select(df_tpm_filter, M3_DN_RNA_S8, M4_DN_RNA_S14)
DPearly_M <- dplyr::select(df_tpm_filter, M1_DPearly_RNA_S1, M3_DPearly_RNA_S9, M4_DPearly_RNA_S15)
DPlate_M <- dplyr::select(df_tpm_filter, M1_DPlate_RNA_S2, M3_DPlate_RNA_S10, M4_DPlate_RNA_S16)
CD4SP_M <- dplyr::select(df_tpm_filter, M1_CD4SP_RNA_S3, M3_CD4SP_RNA_S11, M4_CD4SP_RNA_S17)
CD8SP_M <- dplyr::select(df_tpm_filter, M1_CD8SP_RNA_S4, M3_CD8SP_RNA_S12, M4_CD8SP_RNA_S18)

ETP_M$meanz <- rowMeans(ETP_M)
DN_M$meanz <- rowMeans(DN_M)
DPearly_M$meanz <- rowMeans(DPearly_M)
DPlate_M$meanz <- rowMeans(DPlate_M)
CD4SP_M$meanz <- rowMeans(CD4SP_M)
CD8SP_M$meanz <- rowMeans(CD8SP_M)

# merge and change colnames
homme_means <- dplyr::select(cbind(ETP_M, DN_M$meanz, DPearly_M$meanz, DPlate_M$meanz, CD4SP_M$meanz, CD8SP_M$meanz), -1, -2)
colnames(homme_means) <- c("ETP_M", "DN_M", "DPearly_M", "DPlate_M", "CD4SP_M", "CD8SP_M")

# cbind it all together
TPM_sex_means <- cbind(femme_means, homme_means, turner_columns)

# Calculate Z-scores
df_z <- t(apply(TPM_sex_means[row.names(TPM_sex_means) %in% res_lrt$target_id,],1,function(x) (x-mean(x))/sd(x) ))

# se order
order_colz <- c("ETP_F","ETP_M","ETP_Turner",
                "DN_F","DN_M","DN_Turner",
                "DPearly_F","DPearly_M","DPearly_Turner",
                "DPlate_F","DPlate_M","DPlate_Turner",
                "CD4SP_F","CD4SP_M","CD4SP_Turner",
                "CD8SP_F","CD8SP_M","CD8SP_Turner")

# read in clusters
df_1 <- fread("02_tidy_data/GO_output/Clusters/Thymocyte_timeseries_cluster1.txt", header = F)
df_1$x <- 1
df_2 <- fread("02_tidy_data/GO_output/Clusters/Thymocyte_timeseries_cluster2.txt", header = F)
df_2$x <- 2
df_3 <- fread("02_tidy_data/GO_output/Clusters/Thymocyte_timeseries_cluster3.txt", header = F)
df_3$x <- 3
df_4 <- fread("02_tidy_data/GO_output/Clusters/Thymocyte_timeseries_cluster4.txt", header = F)
df_4$x <- 4
df_5 <- fread("02_tidy_data/GO_output/Clusters/Thymocyte_timeseries_cluster5.txt", header = F)
df_5$x <- 5
df_6 <- fread("02_tidy_data/GO_output/Clusters/Thymocyte_timeseries_cluster6.txt", header = F)
df_6$x <- 6
df_7 <- fread("02_tidy_data/GO_output/Clusters/Thymocyte_timeseries_cluster7.txt", header = F)
df_7$x <- 7

# bind togethter
cl <- rbind(df_1, df_2, df_3, df_4, df_5, df_6, df_7)

# keep only cluster names that are in the Z-values df and vice versa..
heatmap_cl <- cl[cl$V1 %in% rownames(df_z),]
df_z_cut <- df_z[rownames(df_z) %in% cl$V1,]

# order it.
df_z_cut_ordered <- df_z_cut[order(match(rownames(df_z_cut), heatmap_cl$V1)), , drop = FALSE]

# level the factor
heatmap_cl$x <- factor(heatmap_cl$x, levels=c("1","2","3","4","5","6","7"))

# make dataframe for annotation
heatmap_anno <- data.frame(namez = rownames(df_z_cut_ordered), index = 1:nrow(df_z_cut_ordered))

# select names
heatmap_anno_fin <- heatmap_anno[heatmap_anno$namez %in% c( "DNMT3B", "HES1", "HOXA9","LYL1", "MYB",  "MYCN",  "NOTCH1",  "SPI1",  "TLX2",  "BCL11A", "DNMT1",  "LMO2",  "RUNX1",  "STAT4",  "BACH2",  "CCR4",  "CCR5",  "CCR7",  "CD28",  "DNMT3A",  "ETS1",  "FOXP3",  "IKZF1",  "IL2RA",  "IRF4","MAF",  "NFATC1",  "NFATC2",  "RUNX3",  "STAT1",  "STAT2",  "STAT3",  "STAT5A",  "STAT5B",  "STAT6",  "TBX21",  "TNF",  "BCL2L1",  "NFATC3", "RAG1",  "RORC",  "TCF12",  "BCL11B",  "CD3E",  "CD3G",  "E4F1",  "GATA3",  "KDM5A",  "LEF1",  "TET1"),]

# make annotation
ha = rowAnnotation(foo = anno_mark(at = heatmap_anno_fin$index, labels = heatmap_anno_fin$namez))

pdf("04_plots/Figure_8c.pdf",width = 12,height = 16)
set.seed(13)
Heatmap(df_z_cut_ordered, split=heatmap_cl$x, show_row_names = F,column_order = order_colz, col= colorRamp2(-2:2,rev(brewer.pal(5,"RdBu"))), show_row_dend = F, right_annotation = ha, use_raster = T, cluster_row_slices = FALSE, border = T)
dev.off()


# export surce data
exportz <- data.frame(df_z_cut_ordered)
exportz$gene <- rownames(exportz)
colnames(exportz)
exportz <- exportz[,c("gene","ETP_F","DN_F","DPearly_F","DPlate_F","CD4SP_F","CD8SP_F","ETP_M","DN_M","DPearly_M","DPlate_M",
                      "CD4SP_M","CD8SP_M","ETP_Turner","DN_Turner","DPearly_Turner","DPlate_Turner","CD4SP_Turner","CD8SP_Turner")]
write.table(exportz, file = "06_source_data/Fig_8c.tsv", quote = F, row.names = F, sep = "\t")
