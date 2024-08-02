#Load sleuth object for differential expression
library(data.table)
library(sleuth)
library(dplyr)
library(tidyr)
library(cowplot)

source("03_r_scripts/plot_parameters.R")

######## read in data #############
message("loading sleuth object")
source("03_r_scripts/load_scripts/Load_thymocyte_sleuth_expression_data_with_Turner.R")

so_thymo_turner <- so_thymo_genes

df_tpm <- melt(sleuth_to_matrix(so_thymo_turner, 'obs_norm', 'tpm'))

colnames(df_tpm) <- c("target_id", "cell_sample", "TPM")

######## read in annotation #############
# Read in the Tukiainen log2FC data, transform it to long format.
tuki_df <- fread("97_indexes_annotations/Suppl.Table.2.csv")
tuki_df[tuki_df$`Gene name` == "6-Sep"]$`Gene name` <- "SEPT6"
tuki_log2fc_df <- tuki_df %>% dplyr:: select('Gene name', ends_with("logFC"))
melt_tuki_log2fc_df <- melt(tuki_log2fc_df)
melt_tuki_log2fc_df$variable <- gsub(melt_tuki_log2fc_df$variable, pattern = "_logFC", replacement = "")
names(melt_tuki_log2fc_df)[names(melt_tuki_log2fc_df) == 'Gene name'] <- 'gene'

#read in annotations (PAR and the Tukiainen annotation)
PAR <- rbind(fread("97_indexes_annotations/group-715_PAR1.csv"), fread("97_indexes_annotations/group-716_PAR2.csv"))
esc <- fread("97_indexes_annotations/landscape.Suppl.Table.13.csv")
esc[esc$Gene_name == "6-Sep",]$Gene_name <- "SEPT6"
esc <- esc[,-1]
esc <- esc[!(esc$Gene_name == "IDS" & esc$Reported_XCI_status == "Unknown"),]

# merge annotations
chrx_annotation <- merge(esc, PAR, by.x = "Gene_name" , by.y = "Approved symbol" , all = T)[-1,]
chrx_annotation$category <- chrx_annotation$Reported_XCI_status
chrx_annotation[chrx_annotation$PAR == "PAR1" | chrx_annotation$PAR == "PAR2",]$category <- "PAR"
chrx_annotation[chrx_annotation$Gene_name == "XG",]$category <- "PAR"
chrx_annotation[chrx_annotation$category == "Unknown",]$category <- "Potential"
chrx_annotation[chrx_annotation$category == "Variable",]$category <- "Potential"


# keep only chrX genes
df_tpm_filt <- df_tpm[df_tpm$target_id %in% chrx_annotation$Gene_name | df_tpm$target_id == "PUDP",]

# fixes
df_tpm_filt <- df_tpm_filt %>% separate(cell_sample, into = c("sample", "cell", NA, NA), sep = "_")

# add sex
df_tpm_filt$sex <- ifelse(df_tpm_filt$sample %in% c("F1", "F2", "F3", "F4"), yes = "female", no = "male")
df_tpm_filt[df_tpm_filt$sample == "T1",]$sex <- "Turner"

df_tpm_filt_anno <- merge(df_tpm_filt, chrx_annotation, by.x = "target_id", by.y = "Gene_name", all.x = T)

df_tpm_filt_anno[df_tpm_filt_anno$target_id == "PUDP",]$category <- "Escape"
df_tpm_filt_anno[df_tpm_filt_anno$target_id == "SEPT6",]$category <- "Escape"

# find genes not expressed in any subtype
removez <- df_tpm_filt_anno %>% dplyr::group_by(target_id, cell) %>% dplyr::summarise(meanz = mean(TPM))
removez <- removez[removez$meanz > 1,] %>% dplyr::count() %>% subset(n >= 1 )
df_tpm_filt_anno <- df_tpm_filt_anno[df_tpm_filt_anno$target_id %in% removez$target_id,]

# get sex-biased escape genes
lrt_res <- fread("Supplementary Table 5")[seqnames == "chrX" & (target_id %in% df_tpm_filt_anno[df_tpm_filt_anno$category == "Escape",]$target_id | target_id %in% c("PUDP", "SEPT6")) & qval < 0.01 & target_id != "XIST"]

df_tpm_filt_anno <- df_tpm_filt_anno[(df_tpm_filt_anno$category == "Escape" & df_tpm_filt_anno$target_id %in% lrt_res$target_id) | (df_tpm_filt_anno$category == "PAR" | df_tpm_filt_anno$category == "Inactive"),]

df_tpm_filt_anno_stats <- df_tpm_filt_anno %>% dplyr::group_by(category, cell, sex) %>% rstatix::get_summary_stats(TPM, type = "common")

## PLOT ##
top <- ggplot(df_tpm_filt_anno_stats[df_tpm_filt_anno_stats$category != "Potential" ,], aes(x=factor(cell, levels = c("ETP", "DN", "DPearly", "DPlate", "CD4SP", "CD8SP")), y = mean,col = sex)) + geom_pointrange(aes(ymin=mean-se, ymax = mean+se), position = position_dodge(width = .75)) + facet_wrap(~factor(category, levels = c("PAR", "Escape", "Inactive")), ncol = 3) + theme_AL_box_rotX() + labs(x="", y="mean TPM +-SE")# + coord_cartesian(ylim=c(0,300))


ggsave2("04_plots/Figure_4A.pdf",height = 4, width = 8,
        plot_grid(top)
)

