source("03_r_scripts/plot_parameters.R")

library(data.table)
library(ggrepel)
library(dplyr)
library(RColorBrewer)
library(karyoploteR)
library(rtracklayer)
library(cowplot)
library(ggplot2)
library(tidyr)
library(plyr)
library(ggsci)
library(scales)

# Filter out genes not expressed
# Load data
TPM <- fread("Supplementary Table 3")
tpm_melt <- melt(TPM)
tpm_melt <- tpm_melt %>% separate(variable, into = c("cell", "sample", NA, NA), sep = "_")

# Filter out genes not expressed in any subtype
tpm_filter <- tpm_melt %>% dplyr::group_by(gene, cell) %>% dplyr::summarise(meanz = mean(value))
tpm_filter <- subset(tpm_filter, meanz >= 1)


## Read in significance data
res_thymo_genes <- fread("Supplementary Table 5")

# counts
nrow(subset(unique(dplyr::select(res_thymo_genes, target_id, pval)), pval < 0.00001))
nrow(subset(unique(dplyr::select(res_thymo_genes, target_id, pval, seqnames)), pval < 0.00001 & seqnames == "chrX"))

# calc -log10(p)
res_thymo_genes$log10pval <- -log10(res_thymo_genes$pval)

# Read in GCF gene annotation
annotation_gene_table <- subset(fread("97_indexes_annotations/GCF_000001405.38_GRCh38.p12_genomic_with_correct_contig_names.tsv"), type == "gene" & is.na(pseudo))
annotation_gene_table <- annotation_gene_table[grepl(annotation_gene_table$seqnames, pattern = "chr"),]
#chrx_genes <- annotation_gene_table[annotation_gene_table$seqnames == "chrX",]

# Read in Tukiainen data frame
x.esc <- dplyr::select(fread("97_indexes_annotations/landscape.Suppl.Table.13.csv"), Gene_name, Reported_XCI_status, Sex_bias_in_GTEx, XCI_across_tissues, XCI_in_single_cells, Start, End)
esc <- dplyr::select(subset(x.esc), Gene_name, Reported_XCI_status)
esc <- esc[esc$Gene_name != ""]
esc[esc$Gene_name %in% "6-Sep"]$Gene_name <- "SEPT6"
esc <- esc[!(esc$Gene_name %in% "IDS" & esc$Reported_XCI_status == "Unknown")]


# merge our results with the Tukiainen df.
kary_res_thymo_genes <- merge(res_thymo_genes, esc, by.x = c("target_id"), by.y = c("Gene_name"), all.x = T)

kary_res_thymo_genes <- merge(kary_res_thymo_genes, annotation_gene_table[,c("start", "end", "strand", "gene_id")], by.x = "target_id", by.y = "gene_id", all.x = T)

# remove rows with NA contig, keep only columns we need.
kary_genes <- subset(dplyr::select(kary_res_thymo_genes, seqnames, start, end, strand, target_id, log10pval, Reported_XCI_status, qval, b), !is.na(seqnames) )

# Make all chrx genes not in Tukiainen paper unknown escape status.
kary_genes[is.na(kary_genes$Reported_XCI_status) & kary_genes$seqnames == "chrX",]$Reported_XCI_status <- "Unknown"

# I need to remove stuff that is not on the normal contigs
kary_genes <- subset(kary_genes, seqnames %in% chromz)

# Exclude NAs
kary_genes <- subset(kary_genes, !is.na(log10pval))

# extreme vals
extreme_vals <- subset(kary_genes, log10pval > 20)

# extreme vals
sign_genes <- subset(kary_genes, log10pval > 5)

# cap to 25
kary_genes[kary_genes$log10pval > 20]$log10pval <- 20

# make GRanges object
kary_genes_gr <- makeGRangesFromDataFrame(df = kary_genes, seqnames.field = "seqnames", start.field = "start", end.field = "end", strand.field = "strand", keep.extra.columns=TRUE, na.rm = TRUE)

# Remove chrY
kary_genes_gr_no_y <- kary_genes_gr[seqnames(kary_genes_gr) != "chrY"]

# check -log10(p) max values and round up to nearest 10 value.
ymax_genes <- round_any(max(kary_genes_gr_no_y$log10pval), 5, f = ceiling)
#ymax_transcripts <- round_any(max(kary_transcripts_gr_no_y$log10pval), 10, f = ceiling)

# add ablines at -log10(5) and -log10(10). Calculate where to put it.
threshold_0.01_genes <- 10/ymax_genes
threshold_0.05_genes <- 5/ymax_genes

# Calculate how many genes are significant per chromosome (out of total genes we have for each chromosome).
kary_genes_count <- kary_genes[kary_genes$target_id %in% tpm_filter$gene,]
kary_genes_count$sig <- ifelse(kary_genes_count$log10pval > 5 & (kary_genes_count$b > 0.2 | kary_genes_count$b < -0.2), yes = "sig", no = "not_sig")
#kary_genes_count$sig <- ifelse(kary_genes_count$qval < 0.001, yes = "sig", no = "not_sig")
kary_genes_count_kek <- kary_genes_count %>% dplyr::group_by(seqnames, sig) %>% dplyr::count()
kary_genes_count_kek_wide <- dcast(kary_genes_count_kek, seqnames ~ sig, value.var="n")
#kary_genes_count_kek_wide[is.na(kary_genes_count_kek_wide$sig),]$sig <- 0
kary_genes_count_kek_wide$perc <- round(100*(kary_genes_count_kek_wide$sig / (kary_genes_count_kek_wide$not_sig + kary_genes_count_kek_wide$sig)), 2)
kary_genes_count_kek_wide$tot <- kary_genes_count_kek_wide$not_sig + kary_genes_count_kek_wide$sig

# read in chromosome lengths.
#chrom_length <- fread("/media/god/data1/R/thymodevel/chromosome_lengths.tsv")
#chrom_length$mid <- chrom_length$V2/2

#kary_genes_count_kek_wide <- merge(kary_genes_count_kek_wide, chrom_length, by.x = c("seqnames"), by.y = c("gene"))

DT.m1 = data.table::melt(kary_genes_count_kek_wide, id.vars = c("seqnames"))


DT.m1$value <- round(DT.m1$value, 1)
DT.m1$seqnames <- gsub(DT.m1$seqnames, pattern = "chr", replacement = "")
DT.m1$seqnames <- factor(DT.m1$seqnames, levels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"))


library(ggplotify)

pp <- getDefaultPlotParams(4)
pp$bottommargin <- 5
pp$rightmargin <- 0.0050
pp$leftmargin <- 0.04950

fig1b_manhattan <- as.ggplot(expression(kp <- plotKaryotype(plot.type=4, ideogram.plotter = NULL, chromosomes = chromz, plot.params = pp, labels.plotter = NULL),
                                        kpAxis(kp, ymin = 0, ymax = ymax_genes, tick.pos = c(0, 5, 10, 15, 20, 25)),
                                        kpPlotManhattan(kp, data = kary_genes_gr_no_y[kary_genes_gr_no_y$log10pval > 5], pval = kary_genes_gr_no_y[kary_genes_gr_no_y$log10pval > 5]$log10pval, logp = FALSE, ymin=0, ymax=ymax_genes, points.col = chrom_colors, suggestiveline = NULL, genomewideline = NULL, points.cex = 0.5),
                                        kpPlotManhattan(kp, data = kary_genes_gr_no_y[kary_genes_gr_no_y$log10pval < 5], pval = kary_genes_gr_no_y[kary_genes_gr_no_y$log10pval < 5]$log10pval, logp = FALSE, ymin=0, ymax=ymax_genes, points.col = chrom_colors, suggestiveline = NULL, genomewideline = NULL, points.cex = 0.5),
                                        kpAbline(kp, h = threshold_0.05_genes, col="black")))


fig1b_barplot <-
  ggplot(subset(DT.m1, variable == "perc" ), aes(x=seqnames,y=value,fill=seqnames, label = value)) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=chrom_colors_gg) + 
  labs(y="", x = "") + 
  theme_AL_box() + 
  theme(legend.position = "none") + 
  geom_text(size = 2.5, vjust=0.15) + 
  coord_cartesian(ylim=c(0,5)) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.text.x = element_text(size=8), axis.text.y = element_text(size=8))



top_row <- plot_grid(fig1b_manhattan, fig1b_barplot, ncol = 1, rel_heights = c(1,0.5))


ggsave2(filename = "04_plots/Suppl_Figure_2A.pdf", width = 150, 
        height = 100, units = "mm", 
        plot_grid(top_row)
)