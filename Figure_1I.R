library(data.table)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggsci)

source("03_r_scripts/plot_parameters.R")

# read in datas

## Read in sleuth output tables
thymocytes_across <- fread("Supplementary Table 5")
thymocytes_split <- fread("Supplementary Table 6")

# read in TPM data
TPMz <- fread("Supplementary Table 3")

TPMz$meanz <- rowMeans(TPMz[,-c("gene")])
TPMz_meanz <- dplyr::select(TPMz, gene, meanz)
thymocytes_across_meanz <- merge(thymocytes_across, TPMz_meanz, by.x = c("target_id"), by.y = c("gene"))
thymocytes_across_meanz <- thymocytes_across_meanz[thymocytes_across_meanz$meanz > 1,]


TPMz_melt <- melt(TPMz)
TPMz_melt <- TPMz_melt %>% separate(variable, into = c("sample", "cell", NA, NA), sep = "_")

TPMz_melt_meanz <- TPMz_melt %>% dplyr::group_by(gene, cell) %>% dplyr::summarise(meanz = mean(value))

thymocytes_split_tpmz <- merge(thymocytes_split, TPMz_melt_meanz, by.x = c("target_id", "Celltype"), by.y = c("gene", "cell"))

thymocytes_split_tpmz <- thymocytes_split_tpmz[thymocytes_split_tpmz$meanz > 1,]


#read in annotations (PAR and the Tukiainen annotation)
PAR <- rbind(fread("97_indexes_annotations/group-715_PAR1.csv"), fread("97_indexes_annotations/group-716_PAR2.csv"))
PAR$GENES <- PAR$`Approved symbol`
PAR[GENES %in% c("XG"),]$PAR <- "PAR1"

esc <- fread("97_indexes_annotations/landscape.Suppl.Table.13.csv")
esc[esc$Gene_name == "6-Sep",]$Gene_name <- "SEPT6"
esc <- esc[,-1]

# Read in GCF gene annotation
annotation_gene_table <- subset(fread("97_indexes_annotations/GCF_000001405.38_GRCh38.p12_genomic_with_correct_contig_names.tsv"), type == "gene" & is.na(pseudo))
annotation_gene_table <- annotation_gene_table[grepl(annotation_gene_table$seqnames, pattern = "chr"),]

# merge annotations
chrx_annotation <- merge(esc, PAR, by.x = "Gene_name" , by.y = "Approved symbol" , all = T)[-1,]
chrx_annotation$category <- chrx_annotation$Reported_XCI_status
chrx_annotation[chrx_annotation$PAR == "PAR1" | chrx_annotation$PAR == "PAR2",]$category <- "PAR"

chrx_genes <- annotation_gene_table[annotation_gene_table$seqnames == "chrX",]


tukiainen_data <- fread("97_indexes_annotations/Suppl.Table.2.csv")
tukiainen_data_short <- dplyr::select(tukiainen_data, -1, -category, -region)
melt_tukiainen_data_short <- melt(tukiainen_data_short)
melt_tukiainen_data_short <- melt_tukiainen_data_short %>% separate(variable, into = c("tissue_id", "variable"), sep = "_")
colnames(melt_tukiainen_data_short) <- c("target_id", "tissue_id", "variable", "val")
melt_tukiainen_data_short <- melt_tukiainen_data_short[melt_tukiainen_data_short$target_id != "IDS",]

dcast_tukiainen <- dcast(melt_tukiainen_data_short, target_id + tissue_id ~ variable, value.var = "val")
dcast_tukiainen[dcast_tukiainen$target_id == "6-Sep",]$target_id <- "SEPT6"
dcast_tukiainen <- dcast_tukiainen[!is.na(dcast_tukiainen$logFC),]

short_df_our_unsplit <- dplyr::select(thymocytes_across_meanz, target_id, pval, qval)
short_df_our <- dplyr::select(thymocytes_split_tpmz, target_id, Celltype, pval, qval)
short_df_tuki <- dplyr::select(dcast_tukiainen, target_id, tissue_id, Pval, qvalue)
colnames(short_df_our) <- c("target_id", "tissue_id", "Pval", "qvalue")

short_df_our_unsplit$tissue_id <- "THYMOCYES"
colnames(short_df_our_unsplit) <- c("target_id", "Pval", "qvalue", "tissue_id")

short_df_tuki <- short_df_tuki[!is.na(short_df_tuki$Pval),]

plot_df <- unique(rbind(short_df_our, short_df_tuki, short_df_our_unsplit))
plot_df <- plot_df[plot_df$target_id %in% chrx_genes$gene_id]

plot_df$log10pval <- -log10(plot_df$Pval)

plot_df$sig <- ifelse(plot_df$qvalue < 0.01, yes = "sig", no = "ns")
plot_df$sig2 <- ifelse(plot_df$log10pval > 5, yes = "sig", no = "ns")

plot_df_PAR <- plot_df[plot_df$target_id %in% chrx_annotation[chrx_annotation$category == "PAR",]$Gene_name]
plot_df_escape <- plot_df[plot_df$target_id %in% chrx_annotation[chrx_annotation$category == "Escape",]$Gene_name]

plot_df_order_PAR <- plot_df_PAR %>% dplyr::group_by(tissue_id,sig) %>% dplyr::count(.drop = F) %>% subset(sig=="sig")
plot_df_order_PAR <- plot_df_order_PAR[order(plot_df_order_PAR$n, decreasing = T),]$tissue_id

plot_df_order_PAR_fixed <- c(plot_df_order_PAR, "ETP", "DN", "DPearly", "DPlate", "CD4SP", "CD8SP")

plot_df_order_escape <- plot_df_escape %>% dplyr::group_by(tissue_id,sig) %>% dplyr::count(.drop = F) %>% subset(sig=="sig")
plot_df_order_escape <- plot_df_order_escape[order(plot_df_order_escape$n, decreasing = T),]$tissue_id



PAR_test <- plot_df_PAR %>% dplyr::group_by(tissue_id,sig) %>% dplyr::count(.drop = F)

PAR_test_dcast <- dcast(PAR_test, tissue_id~sig, value.var="n")



par_df <- plot_df_PAR[!plot_df_PAR$tissue_id %in% c("ETP", "DN", "DPearly", "DPlate", "CD4SP", "CD8SP"),]
esc_df <- plot_df_escape[!plot_df_escape$tissue_id %in% c("ETP", "DN", "DPearly", "DPlate", "CD4SP", "CD8SP"),]

par_orderz <- par_df %>% dplyr::group_by(tissue_id, sig) %>% dplyr::count()
par_orderz <- dcast(par_orderz, tissue_id ~ sig, value.var = "n") %>% mutate(ratioz = ns/sig)
par_orderz <- par_orderz[order(par_orderz$ratioz, decreasing = F),]

esc_orderz <- esc_df %>% dplyr::group_by(tissue_id, sig) %>% dplyr::count()
esc_orderz <- dcast(esc_orderz, tissue_id ~ sig, value.var = "n") %>% mutate(ratioz = ns/sig)
esc_orderz <- esc_orderz[order(esc_orderz$ratioz, decreasing = F),]

ggsave2(filename = "04_plots/Figure_1I.pdf", height = 14, width = 14,
        plot_grid(ncol=1,
                  ggplot(par_df, aes(x=factor(tissue_id, levels = par_orderz$tissue_id), y=after_stat(count),fill=sig))+
                    geom_bar(stat="count", position = "fill") + 
                    theme_AL_box_rotX() + 
                    labs(x="proportion", title = "PAR genes with significant sex-bias", subtitle = "signifiance threshold = qvalue < 0.01")+ 
                    scale_fill_manual(values=c("sig"="red","ns"="grey")) + 
                    theme(legend.title = element_blank()) + 
                    geom_text(data=par_df, stat="count", aes(label=..count..), position = "fill"),
                  
                  ggplot(esc_df, aes(x=factor(tissue_id, levels = esc_orderz$tissue_id), y=after_stat(count),fill=sig))+
                    geom_bar(stat="count", position = "fill") + 
                    theme_AL_box_rotX() + 
                    labs(x="proportion", title = "Escape genes with significant sex-bias", subtitle = "signifiance threshold = qvalue < 0.01")+ 
                    scale_fill_manual(values=c("sig"="red","ns"="grey")) + 
                    theme(legend.title = element_blank()) + 
                    geom_text(data=esc_df, stat="count", aes(label=..count..), position = "fill")
        )
)

