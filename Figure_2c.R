library(data.table)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggsci)

source("plot_parameters.R")

## Read in sleuth output tables
thymocytes_across <- data.table(read_excel("02_tidy_data/Suppl_tables/Supplementary Tables.xlsx", sheet = 5, skip = 1))
thymocytes_split <- data.table(read_excel("02_tidy_data/Suppl_tables/Supplementary Tables.xlsx", sheet = 6, skip = 1))

# make sure log2FC is numeric
thymocytes_across$log2FC <- as.numeric(thymocytes_across$log2FC)
thymocytes_split$log2FC <- as.numeric(thymocytes_split$log2FC)

# read in TPM data
TPMz <- data.table(read_excel("02_tidy_data/Suppl_tables/Supplementary Tables.xlsx", sheet = 3, skip = 1, col_types = c("text",rep(paste0("numeric"), 46))))

# calc means and keep only genes with mean TPM > 1
TPMz$meanz <- rowMeans(TPMz[,-c("gene")])
TPMz_meanz <- dplyr::select(TPMz, gene, meanz)
thymocytes_across_meanz <- merge(thymocytes_across, TPMz_meanz, by = c("gene"))
thymocytes_across_meanz <- thymocytes_across_meanz[thymocytes_across_meanz$meanz > 1,]

# melt TPMs
TPMz_melt <- melt(TPMz)

# split variable
TPMz_melt <- TPMz_melt %>% separate(variable, into = c("sample", "cell", NA, NA), sep = "_")

# get means
TPMz_melt_meanz <- TPMz_melt %>% dplyr::group_by(gene, cell) %>% dplyr::summarise(meanz = mean(value))

# merge and filter out genes with TPM < 1.
thymocytes_split_tpmz <- merge(thymocytes_split, TPMz_melt_meanz, by.x = c("gene", "celltype"), by.y = c("gene", "cell"))
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

# get chrX genes
chrx_genes <- annotation_gene_table[annotation_gene_table$seqnames == "chrX",]

# get the Tukiainen sex-bias log2FC data
tukiainen_data <- fread("97_indexes_annotations/Suppl.Table.2.csv")

# remove unused columns
tukiainen_data_short <- dplyr::select(tukiainen_data, -1, -category, -region)

# melt
melt_tukiainen_data_short <- melt(tukiainen_data_short)

# separate variable
melt_tukiainen_data_short <- melt_tukiainen_data_short %>% separate(variable, into = c("tissue_id", "variable"), sep = "_")

# change colnames
colnames(melt_tukiainen_data_short) <- c("gene", "tissue_id", "variable", "val")

# remove double entry gene
melt_tukiainen_data_short <- melt_tukiainen_data_short[melt_tukiainen_data_short$gene != "IDS",]

# dcast
dcast_tukiainen <- dcast(melt_tukiainen_data_short, gene + tissue_id ~ variable, value.var = "val")

# change excel conversion error
dcast_tukiainen[dcast_tukiainen$gene == "6-Sep",]$gene <- "SEPT6"

# remove genes without log2FC values
dcast_tukiainen <- dcast_tukiainen[!is.na(dcast_tukiainen$logFC),]

# remove excess columns
short_df_our_unsplit <- dplyr::select(thymocytes_across_meanz, gene, pval, qval)
short_df_our <- dplyr::select(thymocytes_split_tpmz, gene, celltype, pval, qval)
short_df_tuki <- dplyr::select(dcast_tukiainen, gene, tissue_id, Pval, qvalue)

# change colnames to match
colnames(short_df_our) <- c("gene", "tissue_id", "Pval", "qvalue")

# add tissue_id to our thymocyte df
short_df_our_unsplit$tissue_id <- "THYMOCYES"

# change colnames to match
colnames(short_df_our_unsplit) <- c("gene", "Pval", "qvalue", "tissue_id")

# remove rows without pvalues
short_df_tuki <- short_df_tuki[!is.na(short_df_tuki$Pval),]

# bind the dfs together, keeping only chrX genes
plot_df <- unique(rbind(short_df_our, short_df_tuki, short_df_our_unsplit))
plot_df <- plot_df[plot_df$gene %in% chrx_genes$gene_id]

# add a new -log10pval column
plot_df$log10pval <- -log10(plot_df$Pval)

# add significance tag
plot_df$sig <- ifelse(plot_df$qvalue < 0.01, yes = "sig", no = "ns")

# split df into PAR and non-PAR escapees
plot_df_PAR <- plot_df[plot_df$gene %in% chrx_annotation[chrx_annotation$category == "PAR",]$Gene_name]
plot_df_escape <- plot_df[plot_df$gene %in% chrx_annotation[chrx_annotation$category == "Escape",]$Gene_name]

# Remove split thymocytes data.
par_df <- plot_df_PAR[!plot_df_PAR$tissue_id %in% c("ETP", "DN", "DPearly", "DPlate", "CD4SP", "CD8SP"),]
esc_df <- plot_df_escape[!plot_df_escape$tissue_id %in% c("ETP", "DN", "DPearly", "DPlate", "CD4SP", "CD8SP"),]

# create df for order when plotting
par_orderz <- par_df %>% dplyr::group_by(tissue_id, sig) %>% dplyr::count()
par_orderz <- dcast(par_orderz, tissue_id ~ sig, value.var = "n") %>% mutate(ratioz = ns/sig)
par_orderz <- par_orderz[order(par_orderz$ratioz, decreasing = F),]

# create df for order when plotting
esc_orderz <- esc_df %>% dplyr::group_by(tissue_id, sig) %>% dplyr::count()
esc_orderz <- dcast(esc_orderz, tissue_id ~ sig, value.var = "n") %>% mutate(ratioz = ns/sig)
esc_orderz <- esc_orderz[order(esc_orderz$ratioz, decreasing = F),]

# Plot
ggsave2(filename = "04_plots/Figure_2c.pdf", height = 14, width = 14,
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


# Export source data
par_df$group <- "PAR"
esc_df$group <- "Escape"
write.table(rbind(par_df, esc_df), file = "06_source_data/Fig_2c.tsv", quote = F, row.names = F, sep = "\t")
