options(stringsAsFactors = F)

library(cowplot)
library(ggsci)
library(ggrepel)
library(rstatix)
library(ggpubr)
library(ggplot2)
library(sleuth)

######## read in data #############
message("loading sleuth object")
source("03_r_scripts/load_scripts/Load_thymocyte_sleuth_expression_data_with_Turner.R")

sleuth_contrast <- function(so,x,ref){
  new_s2c <- so$sample_to_covariates
  new_s2c$Celltype <- as.factor(new_s2c$Celltype)
  new_s2c$Celltype <- relevel(new_s2c$Celltype, ref = ref)
  
  so2 <- so
  so2$sample_to_covariates <- new_s2c
  new_fit_name <- paste0('ref',ref)
  so2 <- sleuth_fit(so2,~Celltype, fit_name = new_fit_name)
  so2 <- sleuth_wt(so2, paste0('Celltype',x), which_model = new_fit_name)
  return(so2)
}

so_thymo_transcripts <- sleuth_contrast(so_thymo_transcripts,"DN","ETP")
so_thymo_transcripts <- sleuth_contrast(so_thymo_transcripts,"DPearly","DN")
so_thymo_transcripts <- sleuth_contrast(so_thymo_transcripts,"DPlate","DPearly")
so_thymo_transcripts <- sleuth_contrast(so_thymo_transcripts,"CD4SP","DPlate")
so_thymo_transcripts <- sleuth_contrast(so_thymo_transcripts,"CD8SP","DPlate")
so_thymo_transcripts <- sleuth_contrast(so_thymo_transcripts,"CD4SP","CD8SP")


res_cntr <- list(
  "DN_vs_ETP" = sleuth_results(so_thymo_transcripts,'CelltypeDN', which_model = "refETP",test_type = "wt",show_all = F, pval_aggregate = F),
  "DPearly_vs_DN" = sleuth_results(so_thymo_transcripts,'CelltypeDPearly', which_model = "refDN",test_type = "wt",show_all = F, pval_aggregate = F),
  "DPlate_vs_DPearly" =sleuth_results(so_thymo_transcripts,'CelltypeDPlate', which_model = "refDPearly",test_type = "wt",show_all = F, pval_aggregate = F),
  "CD4SP_vs_DPlate" = sleuth_results(so_thymo_transcripts,'CelltypeCD4SP', which_model = "refDPlate",test_type = "wt",show_all = F, pval_aggregate = F),
  "CD8SP_vs_DPlate" = sleuth_results(so_thymo_transcripts,'CelltypeCD8SP', which_model = "refDPlate",test_type = "wt",show_all = F, pval_aggregate = F),
  "CD4SP_vs_CD8SP" = sleuth_results(so_thymo_transcripts,'CelltypeCD4SP', which_model = "refCD8SP",test_type = "wt",show_all = F, pval_aggregate = F)
)

names(res_cntr)

##############################
## Differential gene counts ##
##############################

tpm_filter <- txi_thymo$abundance[apply(txi_thymo$abundance,1,function(x) any(tapply(x,meta_thymo_transition$Celltype,function(y) mean(y>=1) )==1) ),]

tpm_avg <- t(apply(tpm_filter,1,function(x) tapply(x,meta_thymo_transition$Celltype, mean) ))

log2r <- function(res_nm,cutoff=1e-3){
  idx <- unlist(strsplit(res_nm,"_"))
  res <- na.omit(res_cntr[[res_nm]])
  nm <- intersect(res$gene_id[res$qval < cutoff],row.names(tpm_avg))
  x <- tpm_avg[nm,idx[1]]
  y <- tpm_avg[nm,idx[3]]
  # is the gene expressed in any condition?
  isexpr <- x>=1|y>=1
  l2r <- log2(x+1) - log2(y+1)
  return(l2r[isexpr])
}


tpm_diff <- sapply(names(res_cntr),log2r)
df_diff <- melt(tpm_diff)

headtail <- function(x,n) c(head(x,n),tail(x,n))

ETP_DN_names <- names(headtail(sort(tpm_diff$DN_vs_ETP),10))
DPearly_vs_DN_names <- names(headtail(sort(tpm_diff$DPearly_vs_DN),10))
DPlate_vs_DPearly_names <- names(headtail(sort(tpm_diff$DPlate_vs_DPearly),10))
CD4SP_vs_DPlate_names <- names(headtail(sort(tpm_diff$CD4SP_vs_DPlate),10))
CD8SP_vs_DPlate_names <- names(headtail(sort(tpm_diff$CD8SP_vs_DPlate),10))
CD4SP_vs_CD8SP_names <- names(headtail(sort(tpm_diff$CD4SP_vs_CD8SP),10))

# plot boxplots
df_ETP_vs_DN <- subset(dplyr::select(data.frame(tpm_avg), 3, 6), rownames(tpm_avg) %in% ETP_DN_names)
df_DPearly_vs_DN <- subset(dplyr::select(data.frame(tpm_avg), 3, 4), rownames(tpm_avg) %in% DPearly_vs_DN_names)
df_DPlate_vs_DPearly <- subset(dplyr::select(data.frame(tpm_avg), 4, 5), rownames(tpm_avg) %in% DPlate_vs_DPearly_names)
df_CD4SP_vs_DPlate <- subset(dplyr::select(data.frame(tpm_avg), 5, 1), rownames(tpm_avg) %in% CD4SP_vs_DPlate_names)
df_CD8SP_vs_DPlate <- subset(dplyr::select(data.frame(tpm_avg), 5, 2), rownames(tpm_avg) %in% CD8SP_vs_DPlate_names)
df_CD4SP_vs_CD8SP <- subset(dplyr::select(data.frame(tpm_avg), 1, 2), rownames(tpm_avg) %in% CD4SP_vs_CD8SP_names)

df_ETP_vs_DN$gene <- rownames(df_ETP_vs_DN)
df_DPearly_vs_DN$gene <- rownames(df_DPearly_vs_DN)
df_DPlate_vs_DPearly$gene <- rownames(df_DPlate_vs_DPearly)
df_CD4SP_vs_DPlate$gene <- rownames(df_CD4SP_vs_DPlate)
df_CD8SP_vs_DPlate$gene <- rownames(df_CD8SP_vs_DPlate)
df_CD4SP_vs_CD8SP$gene <- rownames(df_CD4SP_vs_CD8SP)

df_ETP_vs_DN$direction <- ifelse(df_ETP_vs_DN$gene %in% names(head(sort(tpm_diff$DN_vs_ETP), 10)), yes = "Up", no = "Down")
df_DPearly_vs_DN$direction <- ifelse(df_DPearly_vs_DN$gene %in% names(head(sort(tpm_diff$DPearly_vs_DN), 10)), yes = "Up", no = "Down")
df_DPlate_vs_DPearly$direction <- ifelse(df_DPlate_vs_DPearly$gene %in% names(head(sort(tpm_diff$DPlate_vs_DPearly), 10)), yes = "Up", no = "Down")
df_CD4SP_vs_DPlate$direction <- ifelse(df_CD4SP_vs_DPlate$gene %in%  names(head(sort(tpm_diff$CD4SP_vs_DPlate), 10)), yes = "Up", no = "Down")
df_CD8SP_vs_DPlate$direction <- ifelse(df_CD8SP_vs_DPlate$gene %in%  names(head(sort(tpm_diff$CD8SP_vs_DPlate), 10)), yes = "Up", no = "Down")
df_CD4SP_vs_CD8SP$direction <- ifelse(df_CD4SP_vs_CD8SP$gene %in%  names(tail(sort(tpm_diff$CD4SP_vs_CD8SP), 10)), yes = "Up", no = "Down")

df_ETP_vs_DN$comparison <- "ETP_vs_DN"
df_DPearly_vs_DN$comparison <- "DN_vs_DPearly"
df_DPlate_vs_DPearly$comparison <- "DPearly_vs_DPlate"
df_CD4SP_vs_DPlate$comparison <- "DPlate_vs_CD4SP"
df_CD8SP_vs_DPlate$comparison <- "DPlate_vs_CD8SP"
df_CD4SP_vs_CD8SP$comparison <- "CD4SP_vs_CD8SP"

smash <- rbind(melt(df_ETP_vs_DN), melt(df_DPearly_vs_DN),melt(df_DPlate_vs_DPearly),melt(df_CD4SP_vs_DPlate),melt(df_CD8SP_vs_DPlate),melt(df_CD4SP_vs_CD8SP))



Turner_tpm <- data.frame(txi_thymo$abundance)

df_Turner <- dplyr::select(Turner_tpm, "T1_ETP_RNA_S1","T1_DN_RNA_S2","T1_DPearly_RNA_S3","T1_DPlate_RNA_S4","T1_CD4SP_RNA_S5", "T1_CD8SP_RNA_S6")
df_Turner$gene <- rownames(df_Turner)
df_Turner <- melt(df_Turner)
df_Turner <- df_Turner %>% tidyr::separate(variable, into = c(NA, "cell", NA, NA), sep = "_")
df_Turner$disease <- "Turner"

asdf <- merge(df_Turner, dplyr::select(smash, -value), by.x = c("gene", "cell"), by.y = c("gene", "variable"))

smash$disease <- "Normal"

colnames(smash) <- c("gene", "direction","comparison","cell", "TPM", "disease")
colnames(asdf) <- c("gene", "cell", "TPM", "disease", "direction", "comparison")

smash_asdf <- rbind(smash, asdf)

smash_asdf$direction_corrected <- "KEK"
smash_asdf[smash_asdf$direction == "Down",]$direction_corrected <- "Up"
smash_asdf[smash_asdf$direction == "Up",]$direction_corrected <- "Down"



stat.test <- compare_means(TPM ~ disease, smash_asdf, group.by = c("cell", "direction_corrected", "comparison"), method = "t.test", p.adjust.method = "BH")

stat.test <- stat.test %>%
  mutate(y.position = 400)

comparions_to_include <- c("DN_vs_DPearly", "DPlate_vs_CD4SP", "DPlate_vs_CD8SP", "CD4SP_vs_CD8SP")

stat.test_kek <- subset(dplyr::select(stat.test, -4, -p, -p.format, -p.signif, -method), comparison %in% comparions_to_include)

plot_df <- smash_asdf[smash_asdf$comparison %in% comparions_to_include,]


filtzz <- intersect(unique(c(names(headtail(sort(tpm_diff$DPearly_vs_DN),10)),
                             names(headtail(sort(tpm_diff$CD4SP_vs_DPlate),10)),
                             names(headtail(sort(tpm_diff$CD8SP_vs_DPlate),10)),
                             names(headtail(sort(tpm_diff$CD4SP_vs_CD8SP),10))
)
), rownames(txi_thymo$abundance[apply(txi_thymo$abundance,1,function(x) any(tapply(x,meta_thymo_transition$Celltype,function(y) mean(y>=1) )==1) ),]
))

genez_to_show <- c("TARP", "JCHAIN", "CD8A", "RORC", "CD8B", "RAG1", "RAG2", "CD40LG", "KLF2", "NOTCH3", "CCR7", "KLRK1", "CD4", "ZBTB7B", "HLA-F", "TOP2A", "CDK1")



ggsave2("04_plots/Suppl_Figure_6C_left.pdf",width = 14,height = 8,
       ggplot(plot_df, aes(x=factor(cell, levels = c("ETP", "DN", "DPearly", "DPlate", "CD4SP", "CD8SP")),
                           y=TPM)) + 
         geom_point(aes(col = disease)) + 
         geom_line(aes(group=interaction(gene, disease), col = disease)) + 
         facet_wrap(factor(direction_corrected, levels = c("Up", "Down"))~factor(comparison, levels = c("ETP_vs_DN","DN_vs_DPearly","DPearly_vs_DPlate","DPlate_vs_CD4SP","DPlate_vs_CD8SP","CD4SP_vs_CD8SP")), scales = "free", ncol = 4, nrow =2) + 
         scale_color_aaas() + 
         theme_AL_box() + 
         theme(legend.title = element_blank()) + 
         labs(x="", title = "Top 10 checkpoint genes per contrast")+
         geom_text_repel(
           data    = plot_df[plot_df$gene %in% genez_to_show,],
           mapping = aes(label = gene), min.segment.length = 0.11111111111111, nudge_y = 1, max.overlaps = 100)
)



###################
## Plot heatmaps ##
###################

plot.heat.cntr <- function(df,res_nm,row_nm,cutoff=1e-3, ... ){
  require(ComplexHeatmap)
  require(RColorBrewer)
  require(circlize)
  strs <- strsplit(res_nm,"_")
  nm <- c(sapply(strs,"[[",1),sapply(strs,"[[",3))
  nm_reg <- paste0(nm,collapse="|")
  res <- res_cntr[[res_nm]]
  gn <- res$gene_id[res$qval < cutoff]
  idx <- intersect(gn,row.names(df))
  idx2 <- tpm_avg[idx,nm[1]]>=1|tpm_avg[idx,nm[2]]>=1
  df_red <- df[idx[idx2],grepl(nm_reg,colnames(df))]
  df_z <- t(apply(df_red,1,function(x) (x-mean(x))/sd(x) ))
  if(!missing(row_nm)) row.names(df_z)[!row.names(df_z) %in% row_nm] <- " "
  Heatmap(df_z, ... )
}

order_colz <- c("F1_ETP_RNA_S5", "F2_ETP_RNA_S1", "F3_ETP_RNA_S7", "F4_ETP_RNA_S19","M3_ETP_RNA_S7", "M4_ETP_RNA_S13" ,"T1_ETP_RNA_S1",
                 "F1_DN_RNA_S6", "F2_DN_RNA_S2" , "F3_DN_RNA_S8" , "F4_DN_RNA_S20"  , "M3_DN_RNA_S8", "M4_DN_RNA_S14", "T1_DN_RNA_S2",
                 "F1_DPearly_RNA_S7", "F2_DPearly_RNA_S3", "F3_DPearly_RNA_S9", "F4_DPearly_RNA_S21", "M1_DPearly_RNA_S1",  "M3_DPearly_RNA_S9", "M4_DPearly_RNA_S15", "T1_DPearly_RNA_S3",
                 "F1_DPlate_RNA_S8", "F2_DPlate_RNA_S4", "F3_DPlate_RNA_S10", "F4_DPlate_RNA_S22", "M1_DPlate_RNA_S2", "M3_DPlate_RNA_S10", "M4_DPlate_RNA_S16", "T1_DPlate_RNA_S4",
                 "F1_CD4SP_RNA_S9" , "F2_CD4SP_RNA_S5", "F3_CD4SP_RNA_S11", "F4_CD4SP_RNA_S23", "M1_CD4SP_RNA_S3", "M3_CD4SP_RNA_S11", "M4_CD4SP_RNA_S17", "T1_CD4SP_RNA_S5",
                 "F1_CD8SP_RNA_S10", "F2_CD8SP_RNA_S6", "F3_CD8SP_RNA_S12", "F4_CD8SP_RNA_S24", "M1_CD8SP_RNA_S4", "M3_CD8SP_RNA_S12", "M4_CD8SP_RNA_S18",  "T1_CD8SP_RNA_S6")


order_colz_DN_DPearly <- order_colz[grepl(order_colz, pattern = paste0(c("DN", "DPearly"), collapse = "|"))]

order_colz_DPearly_DPlate <- order_colz[grepl(order_colz, pattern = paste0(c("DPearly", "DPlate"), collapse = "|"))]

order_colz_DPlate_CD4SP <- order_colz[grepl(order_colz, pattern = paste0(c("DPlate", "CD4SP"), collapse = "|"))]

order_colz_DPlate_CD8SP <- order_colz[grepl(order_colz, pattern = paste0(c("DPlate", "CD8SP"), collapse = "|"))]

order_colz_CD4SP_CD8SP <- order_colz[grepl(order_colz, pattern = paste0(c("CD4SP", "CD8SP"), collapse = "|"))]

pdf("04_plots/Suppl_Figure_6C_right.pdf",width = 6,height = 8)

plot.heat.cntr(tpm_filter,"DPearly_vs_DN", genez_to_show, show_row_dend = F,col= colorRamp2(-2:2,rev(brewer.pal(5,"RdBu"))),name="DPearly_vs_DN",use_raster = T,raster_device="CairoPNG", column_order =  order_colz_DN_DPearly)

plot.heat.cntr(tpm_filter,"CD4SP_vs_DPlate", genez_to_show,show_row_dend = F,col= colorRamp2(-2:2,rev(brewer.pal(5,"RdBu"))),name="CD4SP_vs_DPlate",use_raster = T,raster_device="CairoPNG", column_order =  order_colz_DPlate_CD4SP)

plot.heat.cntr(tpm_filter,"CD8SP_vs_DPlate", genez_to_show,show_row_dend = F, col= colorRamp2(-2:2,rev(brewer.pal(5,"RdBu"))),name="CD8SP_vs_DPlate",use_raster = T,raster_device="CairoPNG", column_order =  order_colz_DPlate_CD8SP)

plot.heat.cntr(tpm_filter,"CD4SP_vs_CD8SP", genez_to_show,show_row_dend = F, col= colorRamp2(-2:2,rev(brewer.pal(5,"RdBu"))),name="CD4SP_vs_CD8SP",use_raster = T,raster_device="CairoPNG", column_order =  order_colz_CD4SP_CD8SP)
dev.off()






