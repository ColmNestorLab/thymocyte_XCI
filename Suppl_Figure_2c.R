# Load data
source("Load_thymocyte_sleuth_expression_data.R")
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(factoextra)
library(cluster)

# Filter out genes not expressed in any subtype
tpm_filter <- txi_thymo_transition$abundance[apply(txi_thymo_transition$abundance,1,function(x) any(tapply(x,meta_thymo_transition$Celltype,function(y) mean(y>=1) )==1) ),]

# get p-values
res_lrt <- data.table(read_excel("Supplementary Data.xlsx", sheet = 3, skip = 1, col_types = c("text", "numeric", "numeric", "numeric", "numeric")))

# Calculate Z-scores
df_z <- t(apply(tpm_filter[row.names(tpm_filter) %in% res_lrt$target_id,],1,function(x) (x-mean(x))/sd(x) ))

## Decide cluster n ##
pdf("04_plots/Thymocyte_timeseries_clustering_nclust_all.pdf")
gridExtra::grid.arrange(
  fviz_nbclust(df_z, kmeans, method = "wss",k.max = 25, nstart=15, iter.max=100),
  fviz_nbclust(df_z, kmeans, method = "silhouette",k.max = 15, nstart=25,iter.max=100),
  fviz_nbclust(df_z, kmeans, method = "gap_stat",nboot=100, k.max = 15, iter.max=100, nstart=25)
)
dev.off()

### NOT RUN: Test cluster n
#Heatmap(df_z,km=4, show_row_names=F) 
#Heatmap(df_z,km=5, show_row_names=F)
#Heatmap(df_z,km=6, show_row_names=F)
#Heatmap(df_z,km=7, show_row_names=F)
#Heatmap(df_z,km=8, show_row_names=F)
#Heatmap(df_z,km=9, show_row_names=F)
#Heatmap(df_z,km=10, show_row_names=F)
#Heatmap(df_z,km=11, show_row_names=F)

# get cluster genes
sapply(6:12,function(x){
  cl <- kmeans(df_z,x,nstart = 25,iter.max = 50)$cluster
  capture.output(split(names(cl),cl), file = paste0("02_tidy_data/GO_output/Clusters",x,"/Thymocyte_timeseries_clusters.txt"))
  sapply(seq_along(split(names(cl),cl)),function(i) write.table(split(names(cl),cl)[[i]],paste0("02_tidy_data/GO_output/Clusters",x,"/Thymocyte_timeseries_cluster",i,".txt"),quote = F,sep = "\t",row.names = F,col.names = F) )
  print(plot.cl(cl))
})

### 7 clusters seems to yield the most information

# Cluster data into 7 clusters
set.seed(13) # for reproducible cluster order
cl <- kmeans(df_z,7,nstart = 25,iter.max = 100)$cluster

# split(names(cl),cl)

# Export list
capture.output(split(names(cl),cl), file = "02_tidy_data/GO_output/Thymocyte_timeseries_clusters.txt")
sapply(seq_along(split(names(cl),cl)),function(i) write.table(split(names(cl),cl)[[i]],paste0("02_tidy_data/GO_output/Thymocyte_timeseries_cluster",i,".txt"),quote = F,sep = "\t",row.names = F,col.names = F) )

# make new df
df_z_for_HM <- df_z

# order
order_colz <- c( "F1_ETP_RNA_S5", "F2_ETP_RNA_S1", "F3_ETP_RNA_S7", "F4_ETP_RNA_S19","M3_ETP_RNA_S7", "M4_ETP_RNA_S13" ,
                 "F1_DN_RNA_S6", "F2_DN_RNA_S2" , "F3_DN_RNA_S8" , "F4_DN_RNA_S20"  , "M3_DN_RNA_S8", "M4_DN_RNA_S14", 
                 "F1_DPearly_RNA_S7", "F2_DPearly_RNA_S3", "F3_DPearly_RNA_S9", "F4_DPearly_RNA_S21", "M1_DPearly_RNA_S1", "M3_DPearly_RNA_S9", "M4_DPearly_RNA_S15",
                 "F1_DPlate_RNA_S8", "F2_DPlate_RNA_S4", "F3_DPlate_RNA_S10", "F4_DPlate_RNA_S22", "M1_DPlate_RNA_S2", "M3_DPlate_RNA_S10", "M4_DPlate_RNA_S16", 
                 "F1_CD4SP_RNA_S9" , "F2_CD4SP_RNA_S5", "F3_CD4SP_RNA_S11", "F4_CD4SP_RNA_S23", "M1_CD4SP_RNA_S3", "M3_CD4SP_RNA_S11", "M4_CD4SP_RNA_S17",
                 "F1_CD8SP_RNA_S10", "F2_CD8SP_RNA_S6", "F3_CD8SP_RNA_S12", "F4_CD8SP_RNA_S24", "M1_CD8SP_RNA_S4", "M3_CD8SP_RNA_S12", "M4_CD8SP_RNA_S18")
                 
# genes to annotate
anno <-c(
"DNMT3B",
"HES1",
"HOXA9",
"LYL1",
"MYB",
"MYCN",
"NOTCH1",
"SPI1",
"TLX2",
"BCL11A",
"DNMT1",
"LMO2",
"RUNX1",
"STAT4",
"BACH2",
"CCR4",
"CCR5",
"CCR7",
"CD28",
"DNMT3A",
"ETS1",
"FOXP3",
"IKZF1",
"IL2RA",
"IRF4",
"MAF",
"NFATC1",
"NFATC2",
"RUNX3",
"STAT1",
"STAT2",
"STAT3",
"STAT5A",
"STAT5B",
"STAT6",
"TBX21",
"TNF",
"BCL2L1",
"NFATC3",
"RAG1",
"RORC",
"TCF12",
"BCL11B",
"CD3E",
"CD3G",
"E4F1",
"GATA3",
"KDM5A",
"LEF1",
"TET1")

#As we want to fix the order of the clusters, we have to re-order the gene-to-cluster assignment as a factor:
split <- factor(cl, levels=c("1","2","3","4","5","6","7"))

# make df for annotation
anno_df <- data.frame(namez = rownames(df_z_for_HM), index = 1:nrow(df_z_for_HM))

# keep only genes in the anno list
anno_df_filt <- anno_df[anno_df$namez %in% anno,]
  
# finalize annotation
ha = rowAnnotation(foo = anno_mark(at = anno_df_filt$index, labels = anno_df_filt$namez))

# Plot
pdf("04_plots/Suppl_Figure_2c_right.pdf",width = 12,height = 16)
set.seed(13)
Heatmap(df_z_for_HM, split=split, show_row_names = F,column_order = order_colz, col= colorRamp2(-2:2,rev(brewer.pal(5,"RdBu"))), show_row_dend = F, right_annotation = ha, use_raster = T, cluster_row_slices = FALSE, border = T)
dev.off()

# Export source data
exportz <- data.frame(df_z_for_HM)
exportz$gene <- rownames(exportz)
exportz <- exportz[,c(41, 1:40)]
write.table(exportz, file = "06_source_data/Suppl_fig_2c_right.tsv", quote = F, row.names = F, sep = "\t")

# add pseudo time
meta_thymo_transition$Pseudotime <- 100
meta_thymo_transition[meta_thymo_transition$Celltype == "ETP"]$Pseudotime <- 1.0
meta_thymo_transition[meta_thymo_transition$Celltype == "DN"]$Pseudotime <- 2.0
meta_thymo_transition[meta_thymo_transition$Celltype == "DPearly"]$Pseudotime <- 3.0
meta_thymo_transition[meta_thymo_transition$Celltype == "DPlate"]$Pseudotime <- 4.0
meta_thymo_transition[meta_thymo_transition$Celltype == "CD4SP"]$Pseudotime <- 4.9
meta_thymo_transition[meta_thymo_transition$Celltype == "CD8SP"]$Pseudotime <- 5.1

## function to plot the timeseries and splits SP stages into CD4/CD8
plot.cl <- function(cl){
  require(ggplot2)
  require(reshape2)
  df_cl <- melt(df_z)
  df_cl$cluster <- cl
  df_cl$celltype <- sapply(strsplit(as.character(df_cl$Var2),"_"),"[[",2)
  df_cl$tp <- meta_thymo_transition$Pseudotime[match(
    df_cl$celltype,
    meta_thymo_transition$Celltype)]
  ggplot() +
    stat_summary(data=df_cl[-grep("CD4",df_cl$celltype),],fun.data="mean_sdl",fun.args = list(mult=1),geom="ribbon", alpha=0.1,aes(y=value,x=round(tp),group=cluster,fill="CD8")) +
    stat_summary(data=df_cl[-grep("CD4",df_cl$celltype),],fun.y="mean",geom="line",alpha=0.5,aes(y=value,x=round(tp),group=cluster,col="CD8")) +
    stat_summary(data=df_cl[-grep("CD4",df_cl$celltype),],fun.y="mean",geom="point",alpha=0.5,aes(y=value,x=round(tp),group=cluster,col="CD8")) +
    stat_summary(data=df_cl[-grep("CD8",df_cl$celltype),],fun.data="mean_sdl",fun.args = list(mult=1),geom="ribbon", alpha=0.1,aes(y=value,x=round(tp),group=cluster, fill="CD4")) +
    stat_summary(data=df_cl[-grep("CD8",df_cl$celltype),],fun.y="mean",geom="line",alpha=0.5,aes(y=value,x=round(tp),group=cluster,col="CD4")) +
    stat_summary(data=df_cl[-grep("CD8",df_cl$celltype),],fun.y="mean",geom="point",alpha=0.5,aes(y=value,x=round(tp),group=cluster,col="CD4")) +
    facet_wrap(~cluster, ncol=1) +
    labs(y="Z-score") +
    scale_x_continuous(labels=c("1"="ETP\nCD34+\nCD1a-","2"="DN\nCD34+\nCD1a+","3"="DPearly\nCD3-","4"="DPlate\nCD3+","5"="SP\nCD3+")) +
    theme_AL_box(legend.position = "top", axis.title.x = element_blank())+
    theme(axis.text.x = element_text(size=4))+
    coord_cartesian(ylim=c(-2,2))+
    geom_hline(yintercept = 0, lty=2)
}



ggsave("04_plots/Suppl_Figure_2c_left_CLUSTER.pdf",width = 6,height = 14, 
       plot_grid(rel_widths = c(0.75,1),plot.cl(cl),NULL))


# export source data
exportz2 <- data.frame(df_z)
exportz2$gene <- rownames(exportz2)
exportz2 <- exportz2[,c(41, 1:40)]
write.table(exportz2, file = "06_source_data/Suppl_fig_2c_left.tsv", quote = F, row.names = F, sep = "\t")

