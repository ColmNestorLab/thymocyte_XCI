library(karyoploteR)
library(data.table)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicFeatures)
library(dplyr)
library(cowplot)
library(ggrepel)
library(rtracklayer)
library(GenomicRanges)
library(readxl)

source("plot_parameters.R")

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
res <- AnnotationDbi::select(txdb, keys(txdb, "TXNAME"),
              columns=c("GENEID","TXNAME"),
              keytype="TXNAME")

# add geneid column
res$GENEID <- rownames(res)

# set names
whatyouwant_transcripts_new <- setNames(res$TXNAME, res$GENEID)

# read in EPIC data
df_beta <- data.frame(fread("02_tidy_data/Suppl_tables/Supplementary_Data_methylation_values.tsv", header = T))

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

# genes to plot that goes in either direction.
genes_to_plot_sense <- c("MSL3", "RAG1", "CD3D", "CD8A")
genes_to_plot_antisense <- c("ARSD", "GEMIN8", "RORC")

# merge vectors
genes_allez <- c(genes_to_plot_sense, genes_to_plot_antisense)

# loop it, plotting karyoplots of each gene mentioned in 'genez_allez'.
for (i in genes_allez){
  
  # set orientation
  kek_test <- ifelse(i %in% genes_to_plot_sense, yes = "left", no = "right")
  
  # make short table
  table_short <- annotation_gene_table_short[annotation_gene_table_short$name2 %in% paste0(i),][1,]
  
  # create zoom region
  table_short$start <- table_short$start_fixed-5000
  table_short$end <- table_short$start_fixed+5000
  zoom.region <- toGRanges(data.frame(table_short$chrom, table_short$start, table_short$end))
  zoom.region.for.kplines <- data.frame(seqnames=table_short$chrom, start = table_short$start, end = table_short$end)
  
  #plot
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
# get annotation table, keep standard contigs and keep only columns we want.
annotation_gene_table_short <- subset(dplyr::select(annotation_gene_table, chrom, txStart, txEnd, name2, name, strand), chrom %in% c(paste0("chr", 1:22), "chrX"))

# fix the TSS depending on strand
annotation_gene_table_short$start_fixed <- ifelse(annotation_gene_table_short$strand == "-", yes = annotation_gene_table_short$txEnd, no  = annotation_gene_table_short$txStart)

# make new df
annotation_gene_table_short_tss1500 <- annotation_gene_table_short

# extend range, for better overlap finding later on.
annotation_gene_table_short_tss1500$start_tss1500 <- annotation_gene_table_short_tss1500$start_fixed - 1500
annotation_gene_table_short_tss1500$end_tss1500 <- annotation_gene_table_short_tss1500$start_fixed + 500

# keep only columns we want & rename them.
annotation_gene_table_short_tss1500 <- annotation_gene_table_short_tss1500[,-c(2,3,6,7)]
colnames(annotation_gene_table_short_tss1500) <- c("seqnames","gene_id","transcript_id","start_tss1500","end_tss1500")

# make GRanges object
tss1500_gr <- makeGRangesFromDataFrame(annotation_gene_table_short_tss1500, seqnames.field=c("seqnames"), start.field=c("start_tss1500"), end.field=c("end_tss1500"), keep.extra.columns = T)

# get probe annotation
anno <- data.frame(IlluminaHumanMethylationEPICanno.ilm10b5.hg38::Other)
anno$probeID <- rownames(anno)

# keep columns we want
anno_short <- dplyr::select(anno[!is.na(anno$CHR_hg38),], CHR_hg38, Start_hg38,  End_hg38, probeID)

# make into GRanges object
anno_gr <- makeGRangesFromDataFrame(anno_short, ignore.strand=T, seqnames.field=c("CHR_hg38"), start.field=c("Start_hg38"),end.field=c("End_hg38"), keep.extra.columns = T)

# overlap probe annotation and the TSS1500 object.
ol.f3 <- findOverlaps(anno_gr, tss1500_gr)
nm.f3 <- tapply(tss1500_gr$transcript_id[subjectHits(ol.f3)],queryHits(ol.f3),function(x) paste0(x,collapse=";") )
anno_short_tss1500 <- data.table(anno_short)
anno_short_tss1500[, transcript_id := NA]
anno_short_tss1500$transcript_id[as.numeric(names(nm.f3))] <- nm.f3

# remove entries with NA transcript ID
anno_short_tss1500 <- anno_short_tss1500[!is.na(anno_short_tss1500$transcript_id),]

# keep only columns we want.
anno_short_tss1500_short <- dplyr::select(anno_short_tss1500, probeID, transcript_id)

# separate transcript_id
split_anno_short_tss1500_short <- anno_short_tss1500_short %>% tidyr::separate(transcript_id, sep=";" , into = c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u"))

# see if it worked.
melt_split_anno_short_tss1500_short <- dplyr::select(subset(melt(split_anno_short_tss1500_short, id.vars = "probeID"), !is.na(value)), -variable)

# merge annotation with the overlap df.
melt_split_anno_short_tss1500_short_annotated <- unique(merge(dplyr::select(annotation_gene_table, name, name2), melt_split_anno_short_tss1500_short, by.x = "name", by.y = "value", allow.cartesian=TRUE))

# genes we want to plot
gois <- c("RAG1", "RORC", "CD3D", "CD3G", "CD8A")


# read in EPIC data
df_beta <- data.frame(fread("02_tidy_data/Suppl_tables/Supplementary_Data_methylation_values.tsv", header = T))

# Read in metadata
metadata <- fread("96_metadata/EPIC_thymocyte_meta.csv")

# melt
melt_df_beta <- melt(df_beta)

# add celltype column
melt_df_beta$cell <- "ETP"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "DN"),]$cell <- "DN"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "DPe"),]$cell <- "DPearly"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "DPl"),]$cell <- "DPlate"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "CD4"),]$cell <- "CD4SP"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "CD8"),]$cell <- "CD8SP"

# add sample column
melt_df_beta$sample <- "F1"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "F3"),]$sample <- "F3"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "F4"),]$sample <- "F4"

# Keep only probes that are -500 to +1500 in distance to the genes of interest. 
probes_methylation <- melt_split_anno_short_tss1500_short_annotated[melt_split_anno_short_tss1500_short_annotated$name2 %in% gois,]

# exclude RAG1 incorrect transcripts
probes_methylation <- probes_methylation[!probes_methylation$name %in% c("NM_001377277", "NM_001377278", "NM_001377279", "NM_001377280"),]

# exclude RORC incorrect transcript
probes_methylation <- probes_methylation[!probes_methylation$name %in% c("NM_005060"),]

# exclude CD8A incorrect transcripts
probes_methylation <- probes_methylation[!(probes_methylation$name2 == "CD8A" & !probes_methylation$name %in% c("NM_171827", "NM_001768", "NR_027353"))]

# merge
melt_df_beta_gois <- merge(x=melt_df_beta, y=unique(probes_methylation[,-1]), by = c("probeID"))

# read in thymocyte TPMs, keep only genes of interest (gois).
kek <- data.table(read_excel("Supplementary Data.xlsx", sheet = 2, skip = 1, col_types = c("text",rep(paste0("numeric"), 46))))
kek_filt <- melt(kek[kek$gene %in% gois,])
# Split variable
kek_filt <- kek_filt %>% tidyr::separate(variable, into = c("sample", "cell", NA, NA), sep = "_")
# remove turner sample
kek_filt <- kek_filt[kek_filt$sample != "T1",]
# add sex tag
kek_filt$sex <- ifelse(kek_filt$sample %in% c("F1", "F2", "F3", "F4"), yes = "female", no = "male")

# plot
TPM_plot <- 
  ggplot(kek_filt, aes(x=factor(cell, levels = c("ETP", "DN", "DPearly", "DPlate", "CD4SP" ,"CD8SP")), y=value)) + 
  geom_boxplot() + 
  facet_wrap(~gene, scales = "free_y", ncol = 1) + 
  labs(x="", y="expression (TPM)") + 
  theme_AL_box_rotX() + 
  scale_y_continuous(limits = c(0, NA))

# calculate median methylation
melt_df_beta_gois_median <- melt_df_beta_gois %>% dplyr::group_by(probeID, name2, cell) %>% dplyr::summarise(medianz = median(value))

#plot
methylation_plot <- 
  ggplot(melt_df_beta_gois_median, aes(x=factor(cell, levels = c("ETP", "DN", "DPearly", "DPlate", "CD4SP" ,"CD8SP")), y=medianz, col = probeID)) + 
  geom_point() +
  geom_line(aes(group=probeID)) +
  facet_wrap(~name2, ncol = 1) + 
  labs(x="", y="methylaion (median beta)") + 
  scale_y_continuous(limits = c(0, NA))+
  theme_AL_box_rotX(legend.position = "none")

ggsave2(filename = "04_plots/Figure_6b_and_suppl_figure_7a_methylation_and_TPM.pdf", height = 14, width = 6,
        plot_grid(ncol=2,
                  TPM_plot, methylation_plot)
)

# export source data
Fig_6b_tpms <- kek_filt[kek_filt$gene %in% c("RAG1", "RORC"),]
Fig_6b_meth_lines <- melt_df_beta_gois_median[melt_df_beta_gois_median$name2 %in% c("RAG1", "RORC"),]

Suppl_fig_7a_tpms <- kek_filt[kek_filt$gene %in% c("CD3G", "CD3D", "CD8A"),]
Suppl_fig_7a_meth_lines <- melt_df_beta_gois_median[melt_df_beta_gois_median$name2 %in% c("CD3G", "CD3D", "CD8A"),]

write.table(Fig_6b_tpms, file = "06_source_data/Fig_6b_tpms.tsv", quote = F, row.names = F, sep = "\t")
write.table(Fig_6b_meth_lines, file = "06_source_data/Fig_6b_lines.tsv", quote = F, row.names = F, sep = "\t")
write.table(Suppl_fig_7a_tpms, file = "06_source_data/Suppl_Fig_7a_tpms.tsv", quote = F, row.names = F, sep = "\t")
write.table(Suppl_fig_7a_meth_lines, file = "06_source_data/Suppl_Fig_7a_lines.tsv", quote = F, row.names = F, sep = "\t")

fig_6b_for_karyoplot1 <- annotation_gene_table_short[annotation_gene_table_short$name2 %in% c("RAG1"),][1,]
fig_6b_for_karyoplot1$start <- fig_6b_for_karyoplot1$start_fixed-5000
fig_6b_for_karyoplot1$end <- fig_6b_for_karyoplot1$start_fixed+5000
fig_6b_for_karyoplot1_CPGs <- cpgs_short[cpgs_short$seqnames %in% fig_6b_for_karyoplot1$chrom & (cpgs_short$start > fig_6b_for_karyoplot1$start & cpgs_short$end < fig_6b_for_karyoplot1$end),]
fig_6b_for_karyoplot1_CPGs$gene <- "RAG1"

fig_6b_for_karyoplot2 <- annotation_gene_table_short[annotation_gene_table_short$name2 %in% c("RORC"),][1,]
fig_6b_for_karyoplot2$start <- fig_6b_for_karyoplot2$start_fixed-5000
fig_6b_for_karyoplot2$end <- fig_6b_for_karyoplot2$start_fixed+5000
fig_6b_for_karyoplot2_CPGs <- cpgs_short[cpgs_short$seqnames %in% fig_6b_for_karyoplot2$chrom & (cpgs_short$start > fig_6b_for_karyoplot2$start & cpgs_short$end < fig_6b_for_karyoplot2$end),]
fig_6b_for_karyoplot2_CPGs$gene <- "RORC"

Fig_6b_meth_karyoplot <- rbind(fig_6b_for_karyoplot1_CPGs, fig_6b_for_karyoplot2_CPGs)

write.table(Fig_6b_meth_karyoplot, file = "06_source_data/Fig_6b_karyoplot.tsv", quote = F, row.names = F, sep = "\t")

fig_7c_for_karyoplot1 <- annotation_gene_table_short[annotation_gene_table_short$name2 %in% c("ARSD"),][1,]
fig_7c_for_karyoplot1$start <- fig_7c_for_karyoplot1$start_fixed-5000
fig_7c_for_karyoplot1$end <- fig_7c_for_karyoplot1$start_fixed+5000
fig_7c_for_karyoplot1_CPGs <- cpgs_short[cpgs_short$seqnames %in% fig_7c_for_karyoplot1$chrom & (cpgs_short$start > fig_7c_for_karyoplot1$start & cpgs_short$end < fig_7c_for_karyoplot1$end),]
fig_7c_for_karyoplot1_CPGs$gene <- "ARSD"

fig_7c_for_karyoplot2 <- annotation_gene_table_short[annotation_gene_table_short$name2 %in% c("GEMIN8"),][1,]
fig_7c_for_karyoplot2$start <- fig_7c_for_karyoplot2$start_fixed-5000
fig_7c_for_karyoplot2$end <- fig_7c_for_karyoplot2$start_fixed+5000
fig_7c_for_karyoplot2_CPGs <- cpgs_short[cpgs_short$seqnames %in% fig_7c_for_karyoplot2$chrom & (cpgs_short$start > fig_7c_for_karyoplot2$start & cpgs_short$end < fig_7c_for_karyoplot2$end),]
fig_7c_for_karyoplot2_CPGs$gene <- "GEMIN8"

Fig_7c_meth_karyoplot <- rbind(fig_7c_for_karyoplot1_CPGs, fig_7c_for_karyoplot2_CPGs)

write.table(Fig_7c_meth_karyoplot, file = "06_source_data/Fig_7c_karyoplot.tsv", quote = F, row.names = F, sep = "\t")


suppl_fig_7a_for_karyoplot1 <- annotation_gene_table_short[annotation_gene_table_short$name2 %in% c("CD3D"),][1,]
suppl_fig_7a_for_karyoplot1$start <- suppl_fig_7a_for_karyoplot1$start_fixed-5000
suppl_fig_7a_for_karyoplot1$end <- suppl_fig_7a_for_karyoplot1$start_fixed+5000
suppl_fig_7a_for_karyoplot1_CPGs <- cpgs_short[cpgs_short$seqnames %in% suppl_fig_7a_for_karyoplot1$chrom & (cpgs_short$start > suppl_fig_7a_for_karyoplot1$start & cpgs_short$end < suppl_fig_7a_for_karyoplot1$end),]
suppl_fig_7a_for_karyoplot1_CPGs$gene <- "CD3D/CD3G"

suppl_fig_7a_for_karyoplot2 <- annotation_gene_table_short[annotation_gene_table_short$name2 %in% c("CD8A"),][1,]
suppl_fig_7a_for_karyoplot2$start <- suppl_fig_7a_for_karyoplot2$start_fixed-5000
suppl_fig_7a_for_karyoplot2$end <- suppl_fig_7a_for_karyoplot2$start_fixed+5000
suppl_fig_7a_for_karyoplot2_CPGs <- cpgs_short[cpgs_short$seqnames %in% suppl_fig_7a_for_karyoplot2$chrom & (cpgs_short$start > suppl_fig_7a_for_karyoplot2$start & cpgs_short$end < suppl_fig_7a_for_karyoplot2$end),]
suppl_fig_7a_for_karyoplot2_CPGs$gene <- "CD8A"

suppl_fig_7a_meth_karyoplot <- rbind(suppl_fig_7a_for_karyoplot1_CPGs, suppl_fig_7a_for_karyoplot2_CPGs)

write.table(suppl_fig_7a_meth_karyoplot, file = "06_source_data/Suppl_fig_7a_karyoplot.tsv", quote = F, row.names = F, sep = "\t")


suppl_fig_7d_for_karyoplot1 <- annotation_gene_table_short[annotation_gene_table_short$name2 %in% c("MSL3"),][1,]
suppl_fig_7d_for_karyoplot1$start <- suppl_fig_7d_for_karyoplot1$start_fixed-5000
suppl_fig_7d_for_karyoplot1$end <- suppl_fig_7d_for_karyoplot1$start_fixed+5000
suppl_fig_7d_for_karyoplot1_CPGs <- cpgs_short[cpgs_short$seqnames %in% suppl_fig_7d_for_karyoplot1$chrom & (cpgs_short$start > suppl_fig_7d_for_karyoplot1$start & cpgs_short$end < suppl_fig_7d_for_karyoplot1$end),]
suppl_fig_7d_for_karyoplot1_CPGs$gene <- "MSL3"

suppl_fig_7d_meth_karyoplot <- suppl_fig_7d_for_karyoplot1_CPGs

write.table(suppl_fig_7d_meth_karyoplot, file = "06_source_data/Suppl_fig_7e_karyoplot.tsv", quote = F, row.names = F, sep = "\t")
