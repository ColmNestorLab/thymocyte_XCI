options(stringsAsFactors = F)

source("03_r_scripts/plot_parameters.R")

library(karyoploteR)
library(rtracklayer)
library(data.table)
library(tidyr)
library(ggplotify)
library(ggplot2)
library(cowplot)
library(grid)
library(rtracklayer)
library(ggplotify)


readcounts <- fread("01_processed_data/Bulk_RNAseq_hetSNP_calling/T1_CD4SP.hetSNP.readcounts.txt")
readcounts_filt <- subset(readcounts, V5 >= 20 & V6 >= 20)

readcounts_filt$stop <- readcounts_filt$V2+1

readcounts_filt <- readcounts_filt[readcounts_filt$V1 %in% chromz,]

readcounts_filt_chrx <- subset(readcounts_filt, V1 == "chrX")

f3 <- fread("01_processed_data/Bulk_RNAseq_hetSNP_calling/F3_CD4SP.hetSNP.readcounts.txt")

f3_filt <- subset(f3, V5 >= 20 & V6 >= 20)

f3_filt$stop <- f3_filt$V2+1

f3_filt <- f3_filt[f3_filt$V1 %in% chromz,]

colnames(f3_filt) <- c("CHROM", "POSITION", "ref", "alt", "ref_readcount", "alt_readcount", "total_readcount", "stop")



## load gene annotations

gtf_raw <- fread("97_indexes_annotations/GCF_000001405.38_GRCh38.p12_genomic_with_correct_contig_names.tsv")

# remove chrY PAR
gtf_raw <- gtf_raw[!(grepl(x=gtf_raw$gene_id, pattern = "_1") & seqnames == "chrY")]


# keep only NM and NR transcripts
patterns <- c("NM", "NR")
result <- dplyr::filter(gtf_raw, grepl(paste(patterns, collapse="|"), transcript_id))

# Exclude psuedo
pseudo_genes <- unique(subset(gtf_raw, pseudo == "TRUE")$gene_id)
result_sans_pseudo <- subset(result, !(gene_id %in% pseudo_genes))

gtf_raw_filt <- gtf_raw[gtf_raw$seqnames %in% c("chr1" , "chr2" , "chr3" , "chr4",  "chr5"  ,"chr6" , "chr7"  ,"chr8" , "chr9" , "chr10" ,"chr11" ,"chr12","chr13", "chr14" ,"chr15" ,"chr16" ,"chr17" ,"chr18" ,"chr19" ,"chr20" ,"chr21", "chr22", "chrX") & gtf_raw$type == "exon",]

gtf <- makeGRangesFromDataFrame(subset(gtf_raw_filt, gene_id %in% result_sans_pseudo$gene_id), keep.extra.columns = T)

# make Granges objects.

# dbsnp
dbsnps <- fread("97_indexes_annotations/short_list_all_chrs.tsv")

dbsnp_turner_hsnps_green <- merge(readcounts_filt, dbsnps, by = c("V1", "V2"), all.x = T )
dbsnp_f3_hsnps_green <- merge(f3_filt, dbsnps, by.x = c("CHROM", "POSITION"), by.y = c("V1", "V2"), all.x = T )

dbsnp_turner_hsnps_green$known <- "KEK"
dbsnp_turner_hsnps_green[is.na(dbsnp_turner_hsnps_green$V3.y),]$known <- "Unknown"
dbsnp_turner_hsnps_green[!is.na(dbsnp_turner_hsnps_green$V3.y),]$known <- "Known"

dbsnp_f3_hsnps_green$known <- "KEK"
dbsnp_f3_hsnps_green[is.na(dbsnp_f3_hsnps_green$V3),]$known <- "Unknown"
dbsnp_f3_hsnps_green[!is.na(dbsnp_f3_hsnps_green$V3),]$known <- "Known"

dbsnp_turner_hsnps_known <- dplyr::select(dbsnp_turner_hsnps_green[dbsnp_turner_hsnps_green$known == "Known",], 1,2,8)
dbsnp_f3_hsnps_known <- dplyr::select(dbsnp_f3_hsnps_green[dbsnp_f3_hsnps_green$known == "Known",], 1,2,8)

colnames(dbsnp_turner_hsnps_known) <- c("seqnames", "start", "end")
colnames(dbsnp_f3_hsnps_known) <- c("seqnames", "start", "end")

# make Granges objects.
dbsnp_turner_hsnps_known_gr <- makeGRangesFromDataFrame(df = dbsnp_turner_hsnps_known, seqnames.field = "seqnames", start.field = "start", end.field = "end", keep.extra.columns=TRUE)

dbsnp_f3_hsnps_known_gr <- makeGRangesFromDataFrame(df = dbsnp_f3_hsnps_known, seqnames.field = "seqnames", start.field = "start", end.field = "end", keep.extra.columns=TRUE)

dbsnp_turner_hsnps_known_gr <- intersect(gtf,dbsnp_turner_hsnps_known_gr, ignore.strand = TRUE)

dbsnp_f3_hsnps_known_gr <- intersect(gtf,dbsnp_f3_hsnps_known_gr, ignore.strand = TRUE)



chroms_sans_x <- c("chr1" , "chr2" , "chr3" , "chr4",  "chr5"  ,"chr6" , "chr7"  ,"chr8" , "chr9" , "chr10" ,"chr11" ,"chr12","chr13", "chr14" ,"chr15" ,"chr16" ,"chr17" ,"chr18" ,"chr19" ,"chr20" ,"chr21", "chr22", "chrX")

plot_col_turner <- as.ggplot(expression(kp <- plotKaryotype(plot.type=2, chromosomes = chroms_sans_x,  main = "Turner"),
                                        kpPlotRegions(kp, data=dbsnp_turner_hsnps_known_gr, col="#A9A9A9", r0 = 0, r1=0.45, avoid.overlapping = F)))

plot_col_f3 <- as.ggplot(expression(kp <- plotKaryotype(plot.type=2, chromosomes = chroms_sans_x, main = "F3"),
                                    kpPlotRegions(kp, data=dbsnp_f3_hsnps_known_gr, col="#A9A9A9", r0 = 0, r1=0.45, avoid.overlapping = F)))



p1 <- plot_grid(plot_col_f3, plot_col_turner, nrow = 1, ncol = 2, labels=c('A', 'B'))
title1 <- ggdraw() + draw_label("hSNPs in exons (readcount both alleles >= 20)", fontface='bold')

ggsave2(filename = "04_plots/Suppl_Figure_6A.pdf", height = 14, width = 12,
        plot_grid(title1, p1, ncol=1, rel_heights=c(0.1, 1)))