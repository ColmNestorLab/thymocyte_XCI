library(karyoploteR)
library(data.table)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicFeatures)
library(dplyr)
library(cowplot)
library(ggrepel)
library(rtracklayer)
library(GenomicRanges)

source("03_r_scripts/plot_parameters.R")

# read in annotation
anno <- data.frame(IlluminaHumanMethylationEPICanno.ilm10b5.hg38::Other)
anno$probeID <- rownames(anno)

# remove excess columns
anno_short <- dplyr::select(anno[!is.na(anno$CHR_hg38),], CHR_hg38, Start_hg38,  End_hg38, probeID)

# convert to GRAnges
anno_gr <- makeGRangesFromDataFrame(anno_short, ignore.strand=T, seqnames.field=c("CHR_hg38"), start.field=c("Start_hg38"),end.field=c("End_hg38"), keep.extra.columns = T)

# Read in annotation table
annotation_gene_table <- fread("97_indexes_annotations/UCSC_refGene.tsv")

# make short, subset standard contigs
annotation_gene_table_short <- subset(dplyr::select(annotation_gene_table, chrom, txStart, txEnd, name2, name, strand), chrom %in% c(paste0("chr", 1:22), "chrX"))

# correct start or stop based on direction of gene
annotation_gene_table_short$start_fixed <- ifelse(annotation_gene_table_short$strand == "-", yes = annotation_gene_table_short$txEnd, no = annotation_gene_table_short$txStart)

# read in txdb object
txdb <- loadDb("97_indexes_annotations/txdb.hg38.refGene.sqlite")

# get transcript and corresponding gene names
res <- select(txdb, keys(txdb, "TXNAME"),
              columns=c("GENEID","TXNAME"),
              keytype="TXNAME")

# add geneid column
res$GENEID <- rownames(res)

# set names
whatyouwant_transcripts_new <- setNames(res$TXNAME, res$GENEID)

# read in EPIC data
df_beta <- data.frame(fread("Supplementary Table 11", header = T))

# Read in metadata
metadata <- fread("96_metadata/EPIC_thymocyte_meta.csv")

# melt
melt_df_beta <- melt(df_beta)

# rename cell types
melt_df_beta$cell <- "ETP"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "DN"),]$cell <- "DN"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "DPe"),]$cell <- "DPearly"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "DPl"),]$cell <- "DPlate"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "CD4"),]$cell <- "CD4SP"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "CD8"),]$cell <- "CD8SP"

# rename samples
melt_df_beta$sample <- "F1"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "F3"),]$sample <- "F3"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "F4"),]$sample <- "F4"

# get medians
melt_df_beta_median <- melt_df_beta %>% dplyr::group_by(probeID, cell) %>% dplyr::summarise(medianz = median(value))

# merge with annotation
melt_df_beta_median_anno <- merge(melt_df_beta_median, anno_short, by = "probeID" )

# make shorted
cpgs_short <- dplyr::select(melt_df_beta_median_anno, CHR_hg38, Start_hg38,  End_hg38, medianz, cell)
colnames(cpgs_short) <- c("seqnames", "start", "end", "median_meth", "cell")

# convert to GRanges
cpgs_gr <- toGRanges(cpgs_short)

# create markers df
markers <- unique(dplyr::select(melt_df_beta_median_anno, CHR_hg38, Start_hg38, End_hg38, probeID))

# change column names
colnames(markers) <- c("chr", "start", "end", "labels")

# make granges
markers_gr <- toGRanges(markers)


# genes to plot that goes in either direction.
genes_to_plot_sense <- c("MSL3", "RAG1", "CD3D", "CD8A")
genes_to_plot_antisense <- c("ARSD", "GEMIN8", "RORC")

# merge vectors
genes_allez <- c(genes_to_plot_sense, genes_to_plot_antisense)

# loop it
for (i in genes_allez){

kek_test <- ifelse(i %in% genes_to_plot_sense, yes = "left", no = "right")
  
table_short <- annotation_gene_table_short[annotation_gene_table_short$name2 %in% paste0(i),][1,]

table_short$start <- table_short$start_fixed-5000
table_short$end <- table_short$start_fixed+5000

zoom.region <- toGRanges(data.frame(table_short$chrom, table_short$start, table_short$end))
zoom.region.for.kplines <- data.frame(seqnames=table_short$chrom, start = table_short$start, end = table_short$end)

pdf(file = paste("04_plots/", i, ".pdf", sep=""), height = 10, width = 12)
kp <- plotKaryotype(chromosomes = table_short$chrom, zoom = zoom.region, plot.type = 3)
#kpPlotMarkers(kp, data=markers_gr, data.panel = 2, text.orientation = "vertical", max.iter = 1000, cex = 0.5, label.dist = 0.0001, label.margin=0.0000001, adjust.label.position = T, r0 = 0.30)
kpPlotGenes(kp, data=txdb, add.transcript.names = T, plot.transcripts = T, transcript.names = whatyouwant_transcripts_new, r0=0.175, r1 = 0.75, data.panel = 2, add.strand.marks = F, add.gene.names = F, transcript.name.position = kek_test)
kpAddBaseNumbers(kp, tick.dist = 2000, minor.tick.dist = 2000, add.units = TRUE, cex = 0.8, digits = 8)

at <- autotrack(current.track = 6, total.tracks = 6)
kpSegments(kp, data = cpgs_gr[cpgs_gr$cell == "ETP" & cpgs_gr$median_meth <= 0.25], y1=cpgs_gr[cpgs_gr$cell == "ETP"& cpgs_gr$median_meth <= 0.25]$median_meth, y0=0, data.panel = 1, r0 = at$r0, r1 = at$r1, col = "blue", lwd = 5)
kpSegments(kp, data = cpgs_gr[cpgs_gr$cell == "ETP" & cpgs_gr$median_meth > 0.25], y1=cpgs_gr[cpgs_gr$cell == "ETP"& cpgs_gr$median_meth > 0.25]$median_meth, y0=0, data.panel = 1, r0 = at$r0, r1 = at$r1, col = "black", lwd = 5)
kpAxis(kp, ymin=0, ymax=1, r0 = at$r0, r1 = at$r1, data.panel = 1)
kpAxis(kp, ymin=0, ymax=1, r0 = at$r0, r1 = at$r1, data.panel = 1)
kpLines(kp, chr=zoom.region.for.kplines$seqnames, x=zoom.region.for.kplines$start:zoom.region.for.kplines$end, y=0, r0 = at$r0, r1 = at$r1, data.panel = 1)
kpLines(kp, chr=zoom.region.for.kplines$seqnames, x=zoom.region.for.kplines$start:zoom.region.for.kplines$end, y=0.5, r0 = at$r0, r1 = at$r1, data.panel = 1, lty=2)
kpText(kp, chr=table_short$chrom, y=0.5, x=table_short$end-500, labels="ETP", r0 = at$r0, r1 = at$r1, data.panel = 1, pch=".", cex=1)


at <- autotrack(current.track = 5, total.tracks = 6)
kpSegments(kp, data = cpgs_gr[cpgs_gr$cell == "DN" & cpgs_gr$median_meth <= 0.25], y1=cpgs_gr[cpgs_gr$cell == "DN"& cpgs_gr$median_meth <= 0.25]$median_meth, y0=0, data.panel = 1, r0 = at$r0, r1 = at$r1, col = "blue", lwd = 5)
kpSegments(kp, data = cpgs_gr[cpgs_gr$cell == "DN" & cpgs_gr$median_meth > 0.25], y1=cpgs_gr[cpgs_gr$cell == "DN"& cpgs_gr$median_meth > 0.25]$median_meth, y0=0, data.panel = 1, r0 = at$r0, r1 = at$r1, col = "black", lwd = 5)
kpAxis(kp, ymin=0, ymax=1, r0 = at$r0, r1 = at$r1, data.panel = 1)
kpLines(kp, chr=zoom.region.for.kplines$seqnames, x=zoom.region.for.kplines$start:zoom.region.for.kplines$end, y=0, r0 = at$r0, r1 = at$r1, data.panel = 1)
kpLines(kp, chr=zoom.region.for.kplines$seqnames, x=zoom.region.for.kplines$start:zoom.region.for.kplines$end, y=0.5, r0 = at$r0, r1 = at$r1, data.panel = 1, lty=2)
kpText(kp, chr=table_short$chrom, y=0.5, x=table_short$end-500, labels="DN", r0 = at$r0, r1 = at$r1, data.panel = 1, pch=".", cex=1)


at <- autotrack(current.track = 4, total.tracks = 6)
kpSegments(kp, data = cpgs_gr[cpgs_gr$cell == "DPearly" & cpgs_gr$median_meth <= 0.25], y1=cpgs_gr[cpgs_gr$cell == "DPearly"& cpgs_gr$median_meth <= 0.25]$median_meth, y0=0, data.panel = 1, r0 = at$r0, r1 = at$r1, col = "blue", lwd = 5)
kpSegments(kp, data = cpgs_gr[cpgs_gr$cell == "DPearly" & cpgs_gr$median_meth > 0.25], y1=cpgs_gr[cpgs_gr$cell == "DPearly"& cpgs_gr$median_meth > 0.25]$median_meth, y0=0, data.panel = 1, r0 = at$r0, r1 = at$r1, col = "black", lwd = 5)
kpAxis(kp, ymin=0, ymax=1, r0 = at$r0, r1 = at$r1, data.panel = 1)
kpAxis(kp, ymin=0, ymax=1, r0 = at$r0, r1 = at$r1, data.panel = 1)
kpLines(kp, chr=zoom.region.for.kplines$seqnames, x=zoom.region.for.kplines$start:zoom.region.for.kplines$end, y=0, r0 = at$r0, r1 = at$r1, data.panel = 1)
kpLines(kp, chr=zoom.region.for.kplines$seqnames, x=zoom.region.for.kplines$start:zoom.region.for.kplines$end, y=0.5, r0 = at$r0, r1 = at$r1, data.panel = 1, lty=2)
kpText(kp, chr=table_short$chrom, y=0.5, x=table_short$end-500, labels="DPearly", r0 = at$r0, r1 = at$r1, data.panel = 1, pch=".", cex=1)

at <- autotrack(current.track = 3, total.tracks = 6)
kpSegments(kp, data = cpgs_gr[cpgs_gr$cell == "DPlate" & cpgs_gr$median_meth <= 0.25], y1=cpgs_gr[cpgs_gr$cell == "DPlate"& cpgs_gr$median_meth <= 0.25]$median_meth, y0=0, data.panel = 1, r0 = at$r0, r1 = at$r1, col = "blue", lwd = 5)
kpSegments(kp, data = cpgs_gr[cpgs_gr$cell == "DPlate" & cpgs_gr$median_meth > 0.25], y1=cpgs_gr[cpgs_gr$cell == "DPlate"& cpgs_gr$median_meth > 0.25]$median_meth, y0=0, data.panel = 1, r0 = at$r0, r1 = at$r1, col = "black", lwd = 5)
kpAxis(kp, ymin=0, ymax=1, r0 = at$r0, r1 = at$r1, data.panel = 1)
kpAxis(kp, ymin=0, ymax=1, r0 = at$r0, r1 = at$r1, data.panel = 1)
kpLines(kp, chr=zoom.region.for.kplines$seqnames, x=zoom.region.for.kplines$start:zoom.region.for.kplines$end, y=0, r0 = at$r0, r1 = at$r1, data.panel = 1)
kpLines(kp, chr=zoom.region.for.kplines$seqnames, x=zoom.region.for.kplines$start:zoom.region.for.kplines$end, y=0.5, r0 = at$r0, r1 = at$r1, data.panel = 1, lty=2)
kpText(kp, chr=table_short$chrom, y=0.5, x=table_short$end-500, labels="DPlate", r0 = at$r0, r1 = at$r1, data.panel = 1, pch=".", cex=1)

at <- autotrack(current.track = 2, total.tracks = 6)
kpSegments(kp, data = cpgs_gr[cpgs_gr$cell == "CD4SP" & cpgs_gr$median_meth <= 0.25], y1=cpgs_gr[cpgs_gr$cell == "CD4SP"& cpgs_gr$median_meth <= 0.25]$median_meth, y0=0, data.panel = 1, r0 = at$r0, r1 = at$r1, col = "blue", lwd = 5)
kpSegments(kp, data = cpgs_gr[cpgs_gr$cell == "CD4SP" & cpgs_gr$median_meth > 0.25], y1=cpgs_gr[cpgs_gr$cell == "CD4SP"& cpgs_gr$median_meth > 0.25]$median_meth, y0=0, data.panel = 1, r0 = at$r0, r1 = at$r1, col = "black", lwd = 5)
kpAxis(kp, ymin=0, ymax=1, r0 = at$r0, r1 = at$r1, data.panel = 1)
kpAxis(kp, ymin=0, ymax=1, r0 = at$r0, r1 = at$r1, data.panel = 1)
kpLines(kp, chr=zoom.region.for.kplines$seqnames, x=zoom.region.for.kplines$start:zoom.region.for.kplines$end, y=0, r0 = at$r0, r1 = at$r1, data.panel = 1)
kpLines(kp, chr=zoom.region.for.kplines$seqnames, x=zoom.region.for.kplines$start:zoom.region.for.kplines$end, y=0.5, r0 = at$r0, r1 = at$r1, data.panel = 1, lty=2)
kpText(kp, chr=table_short$chrom, y=0.5, x=table_short$end-500, labels="CD4SP", r0 = at$r0, r1 = at$r1, data.panel = 1, pch=".", cex=1)

at <- autotrack(current.track = 1, total.tracks = 6)
kpSegments(kp, data = cpgs_gr[cpgs_gr$cell == "CD8SP" & cpgs_gr$median_meth <= 0.25], y1=cpgs_gr[cpgs_gr$cell == "CD8SP"& cpgs_gr$median_meth <= 0.25]$median_meth, y0=0, data.panel = 1, r0 = at$r0, r1 = at$r1, col = "blue", lwd = 5)
kpSegments(kp, data = cpgs_gr[cpgs_gr$cell == "CD8SP" & cpgs_gr$median_meth > 0.25], y1=cpgs_gr[cpgs_gr$cell == "CD8SP"& cpgs_gr$median_meth > 0.25]$median_meth, y0=0, data.panel = 1, r0 = at$r0, r1 = at$r1, col = "black", lwd = 5)
kpAxis(kp, ymin=0, ymax=1, r0 = at$r0, r1 = at$r1, data.panel = 1)
kpAxis(kp, ymin=0, ymax=1, r0 = at$r0, r1 = at$r1, data.panel = 1)
kpLines(kp, chr=zoom.region.for.kplines$seqnames, x=zoom.region.for.kplines$start:zoom.region.for.kplines$end, y=0, r0 = at$r0, r1 = at$r1, data.panel = 1)
kpLines(kp, chr=zoom.region.for.kplines$seqnames, x=zoom.region.for.kplines$start:zoom.region.for.kplines$end, y=0.5, r0 = at$r0, r1 = at$r1, data.panel = 1, lty=2)
kpText(kp, chr=table_short$chrom, y=0.5, x=table_short$end-500, labels="CD8SP", r0 = at$r0, r1 = at$r1, data.panel = 1, pch=".", cex=1)
dev.off()

}



######################## Box- and lineplots #####################
################# Read in and prepare data ####################
# Load libraries
annotation_gene_table_short <- subset(dplyr::select(annotation_gene_table, chrom, txStart, txEnd, name2, name, strand), chrom %in% c(paste0("chr", 1:22), "chrX"))

annotation_gene_table_short$start_fixed <- ifelse(annotation_gene_table_short$strand == "-", yes = annotation_gene_table_short$txEnd, no  = annotation_gene_table_short$txStart)

annotation_gene_table_short_tss1500 <- annotation_gene_table_short



anno <- data.frame(IlluminaHumanMethylationEPICanno.ilm10b5.hg38::Other)
anno$probeID <- rownames(anno)

anno_short <- dplyr::select(anno[!is.na(anno$CHR_hg38),], CHR_hg38, Start_hg38,  End_hg38, probeID)

anno_gr <- makeGRangesFromDataFrame(anno_short, ignore.strand=T, seqnames.field=c("CHR_hg38"), start.field=c("Start_hg38"),end.field=c("End_hg38"), keep.extra.columns = T)



annotation_gene_table_short_tss1500$start_tss1500 <- annotation_gene_table_short_tss1500$start_fixed - 1500
annotation_gene_table_short_tss1500$end_tss1500 <- annotation_gene_table_short_tss1500$start_fixed + 500

annotation_gene_table_short_tss1500 <- annotation_gene_table_short_tss1500[,-c(2,3,6,7)]

colnames(annotation_gene_table_short_tss1500) <- c("seqnames","gene_id","transcript_id","start_tss1500","end_tss1500")

tss1500_gr <- makeGRangesFromDataFrame(annotation_gene_table_short_tss1500, seqnames.field=c("seqnames"), start.field=c("start_tss1500"), end.field=c("end_tss1500"), keep.extra.columns = T)

ol.f3 <- findOverlaps(anno_gr, tss1500_gr)
nm.f3 <- tapply(tss1500_gr$transcript_id[subjectHits(ol.f3)],queryHits(ol.f3),function(x) paste0(x,collapse=";") )

anno_short_tss1500 <- data.table(anno_short)
anno_short_tss1500[, transcript_id := NA]
anno_short_tss1500$transcript_id[as.numeric(names(nm.f3))] <- nm.f3

anno_short_tss1500 <- anno_short_tss1500[!is.na(anno_short_tss1500$transcript_id),]


testing_short <- dplyr::select(anno_short_tss1500, probeID, transcript_id)


luuul <- testing_short %>% tidyr::separate(transcript_id, sep=";" , into = c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u"))


did_it_really_work <- dplyr::select(subset(melt(luuul, id.vars = "probeID"), !is.na(value)), -variable)


WOW_did_it_really_work <- unique(merge(dplyr::select(annotation_gene_table, name, name2), did_it_really_work, by.x = "name", by.y = "value", allow.cartesian=TRUE))

pois <- c("RAG1", "RORC", "CD3D", "CD3G", "CD8A")


# read in EPIC data
df_beta <- data.frame(fread("Supplementary Table 11", header = T))

# Read in metadata
metadata <- fread("96_metadata/EPIC_thymocyte_meta.csv")

# melt
melt_df_beta <- melt(df_beta)

melt_df_beta$cell <- "ETP"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "DN"),]$cell <- "DN"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "DPe"),]$cell <- "DPearly"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "DPl"),]$cell <- "DPlate"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "CD4"),]$cell <- "CD4SP"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "CD8"),]$cell <- "CD8SP"

melt_df_beta$sample <- "F1"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "F3"),]$sample <- "F3"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "F4"),]$sample <- "F4"



probes_methylation <- WOW_did_it_really_work[WOW_did_it_really_work$name2 %in% pois,]

# exclude RAG1 incorrect transcripts
probes_methylation <- probes_methylation[!probes_methylation$name %in% c("NM_001377277", "NM_001377278", "NM_001377279", "NM_001377280"),]

# exclude RORC incorrect transcript
probes_methylation <- probes_methylation[!probes_methylation$name %in% c("NM_005060"),]

# exclude CD8A incorrect transcripts
probes_methylation <- probes_methylation[!(probes_methylation$name2 == "CD8A" & !probes_methylation$name %in% c("NM_171827", "NM_001768", "NR_027353"))]


melt_df_beta_pois <- merge(x=melt_df_beta, y=unique(probes_methylation[,-1]), by = c("probeID"))


kek <- fread("Supplementary Table 3")
kek_filt <- melt(kek[kek$gene %in% pois,])
kek_filt <- kek_filt %>% tidyr::separate(variable, into = c("sample", "cell", NA, NA), sep = "_")
kek_filt <- kek_filt[kek_filt$sample != "T1",]
kek_filt$sex <- ifelse(kek_filt$sample %in% c("F1", "F2", "F3", "F4"), yes = "female", no = "male")


tpmzzz <- 
  ggplot(kek_filt, aes(x=factor(cell, levels = c("ETP", "DN", "DPearly", "DPlate", "CD4SP" ,"CD8SP")), y=value)) + 
  geom_boxplot() + 
  facet_wrap(~gene, scales = "free_y", ncol = 1) + 
  labs(x="", y="expression (TPM)") + 
  theme_AL_box_rotX() + 
  scale_y_continuous(limits = c(0, NA))


melt_df_beta_pois_median <- melt_df_beta_pois %>% dplyr::group_by(probeID, name2, cell) %>% dplyr::summarise(medianz = median(value))


methhh <- 
  ggplot(melt_df_beta_pois_median, aes(x=factor(cell, levels = c("ETP", "DN", "DPearly", "DPlate", "CD4SP" ,"CD8SP")), y=medianz, col = probeID)) + 
  geom_point() +
  geom_line(aes(group=probeID)) +
  facet_wrap(~name2, ncol = 1) + 
  labs(x="", y="methylaion (median beta)") + 
  scale_y_continuous(limits = c(0, NA))+
  theme_AL_box_rotX(legend.position = "none")

ggsave2(filename = "04_plots/Figure_3B_and_suppl_figure_5A_methylation_and_TPM.pdf", height = 14, width = 6,
        plot_grid(ncol=2,
                  tpmzzz, methhh)
)
