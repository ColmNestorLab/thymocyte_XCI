
library(data.table)
library(cowplot)
library(dplyr)
library(ggrepel)
library(rtracklayer)
library(GenomicRanges)
source("plot_parameters.R")

# read in annotation able
annotation_gene_table <- fread("97_indexes_annotations/UCSC_refGene.tsv")

# remove excess columns
annotation_gene_table_short <- subset(dplyr::select(annotation_gene_table, chrom, txStart, txEnd, name2, name, strand), chrom %in% c(paste0("chr", 1:22), "chrX"))

# remove strange transcripts
annotation_gene_table_short <- annotation_gene_table_short[!annotation_gene_table_short$name %in% c("NM_001320753", "NM_000351", "NM_001320754"),]

# correct start or stop based on direction of gene
annotation_gene_table_short$start_fixed <- ifelse(annotation_gene_table_short$strand == "-", yes = annotation_gene_table_short$txEnd, no  = annotation_gene_table_short$txStart)

# make new table
annotation_gene_table_short_tss500 <- annotation_gene_table_short

# extend range to TSS-500.
annotation_gene_table_short_tss500$first_tss500 <- ifelse(annotation_gene_table_short_tss500$strand == "+", yes = annotation_gene_table_short_tss500$start_fixed - 500, no = annotation_gene_table_short_tss500$start_fixed + 500)

# add another column to make sure orientation of gene is correct further down the line,
annotation_gene_table_short_tss500$second_tss500 <- annotation_gene_table_short_tss500$start_fixed + 1
annotation_gene_table_short_tss500 <- transform(annotation_gene_table_short_tss500, start_tss500 = pmin(first_tss500, second_tss500))
annotation_gene_table_short_tss500 <- transform(annotation_gene_table_short_tss500, end_tss500 = pmax(first_tss500, second_tss500))

# keep only columns we want, rename them,
annotation_gene_table_short_tss500 <- dplyr::select(annotation_gene_table_short_tss500, chrom, name2 ,name, start_tss500, end_tss500)
colnames(annotation_gene_table_short_tss500) <- c("seqnames","gene_id","transcript_id","start_tss500","end_tss500")

# make GRanges object.
tss500_gr <- makeGRangesFromDataFrame(annotation_gene_table_short_tss500, seqnames.field=c("seqnames"), start.field=c("start_tss500"), end.field=c("end_tss500"), keep.extra.columns = T)

# read in EPIC probe annotation,
anno <- data.frame(IlluminaHumanMethylationEPICanno.ilm10b5.hg38::Other)
anno$probeID <- rownames(anno)

# remove excess columns.
anno_short <- dplyr::select(anno[!is.na(anno$CHR_hg38),], CHR_hg38, Start_hg38,  End_hg38, probeID)

# make GRanges object.
anno_gr <- makeGRangesFromDataFrame(anno_short, ignore.strand=T, seqnames.field=c("CHR_hg38"), start.field=c("Start_hg38"),end.field=c("End_hg38"), keep.extra.columns = T)


# overlap probe annotation and the TSS-500 object.
ol.f3 <- findOverlaps(anno_gr, tss500_gr)
nm.f3 <- tapply(tss500_gr$transcript_id[subjectHits(ol.f3)],queryHits(ol.f3),function(x) paste0(x,collapse=";") )
anno_short_tss500 <- data.table(anno_short)
anno_short_tss500[, transcript_id := NA]
anno_short_tss500$transcript_id[as.numeric(names(nm.f3))] <- nm.f3

# remove entries with NA transcript ID
anno_short_tss500 <- anno_short_tss500[!is.na(anno_short_tss500$transcript_id),]

# keep only columns we want.
anno_short_tss500_short <- dplyr::select(anno_short_tss500, probeID, transcript_id)

# separate transcript_id
split_anno_short_tss500_short <- anno_short_tss500_short %>% tidyr::separate(transcript_id, sep=";" , into = c("a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u"))

# see if it worked.
melt_split_anno_short_tss500_short <- dplyr::select(subset(melt(split_anno_short_tss500_short, id.vars = "probeID"), !is.na(value)), -variable)

# merge annotation with the overlap df.
melt_split_anno_short_tss500_short_annotated <- unique(merge(dplyr::select(annotation_gene_table, name, name2), melt_split_anno_short_tss500_short, by.x = "name", by.y = "value", allow.cartesian=TRUE))




# get all inactive F3 genes
df_ase <- data.table(read_excel("02_tidy_data/Suppl_tables/Supplementary Tables.xlsx", sheet = 10, skip = 1, col_types = c("text", "numeric", "text", "text", "text", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "logical", "logical")))

df_ase <- subset(df_ase, contig == "chrX" & sample == "F3" & keep == T)

# create classification df.
our_classification <- unique(dplyr::select(df_ase, gene, effectSize))
our_classification <- our_classification %>% dplyr::group_by(gene) %>% dplyr::summarise(mean_ASE = mean(effectSize))
our_classification$category <- ifelse(our_classification$mean_ASE > 0.4, yes = "Inactive", no = "Escape")
our_classification <-our_classification[!grepl(our_classification$gene, pattern = ";"),]
ASE_classification <- subset(our_classification, !(gene == "PRKX" & category == "Inactive") & !(gene == "PUDP" & category == "Inactive") & !(gene == "TXLNG" & category == "Inactive"))

# read in PAR and Tukiainen dfs
PAR <- rbind(fread("97_indexes_annotations/group-715_PAR1.csv"), fread("97_indexes_annotations/group-716_PAR2.csv"))
esc <- fread("97_indexes_annotations/landscape.Suppl.Table.13.csv")
esc[esc$Gene_name == "6-Sep",]$Gene_name <- "SEPT6"
esc <- esc[esc$Gene_ID != "ENSG00000241489.3",]
esc <- esc[,-1]

# merge
chrx_annotation <- merge(esc, PAR, by.x = "Gene_name" , by.y = "Approved symbol" , all = T)[-1,]
chrx_annotation$category <- chrx_annotation$Reported_XCI_status
chrx_annotation[chrx_annotation$PAR == "PAR1" | chrx_annotation$PAR == "PAR2",]$category <- "PAR"

# rename known genes
chrx_annotation[chrx_annotation$Gene_name == "KAL1",]$Gene_name <- "ANOS1"
chrx_annotation[chrx_annotation$Gene_name == "CXorf36",]$Gene_name <- "DIPK2B"
chrx_annotation[chrx_annotation$Gene_name == "HDHD1",]$Gene_name <- "PUDP"
chrx_annotation[chrx_annotation$Gene_name == "RGAG4",]$Gene_name <- "RTL5"


# find ecapees and inactive genes that overlap between our and Tukiainen classification.
pois_esc <- intersect(ASE_classification[ASE_classification$category %in% "Escape",]$gene, chrx_annotation[chrx_annotation$Reported_XCI_status == "Escape",]$Gene_name)
pois_inacs <- intersect(ASE_classification[ASE_classification$category %in% "Inactive",]$gene, chrx_annotation[chrx_annotation$Reported_XCI_status == "Inactive",]$Gene_name)

# Name some other cool genes.
pois_others <- c("KDM6A", "DDX3X", "XIST", "BCOR")
pois_others2 <- c("ITM2A", "CXCR3", "CD40LG", "TLR7")


# Read in methylation data.
df_beta <- data.frame(fread("02_tidy_data/Suppl_tables/Supplementary_Data_methylation_values.tsv", header = T))

# Read in metadata
metadata <- fread("96_metadata/EPIC_thymocyte_meta.csv")

# melt
melt_df_beta <- melt(df_beta)

# add celltype
melt_df_beta$cell <- "ETP"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "DN"),]$cell <- "DN"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "DPe"),]$cell <- "DPearly"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "DPl"),]$cell <- "DPlate"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "CD4"),]$cell <- "CD4SP"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "CD8"),]$cell <- "CD8SP"

# add sample
melt_df_beta$sample <- "F1"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "F3"),]$sample <- "F3"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "F4"),]$sample <- "F4"

# get probes that are within tss-500 for the genes above.
probes_methylation_esc <- melt_split_anno_short_tss500_short_annotated[melt_split_anno_short_tss500_short_annotated$name2 %in% pois_esc,]
probes_methylation_inacs <- melt_split_anno_short_tss500_short_annotated[melt_split_anno_short_tss500_short_annotated$name2 %in% pois_inacs,]
probes_methylation_others <- melt_split_anno_short_tss500_short_annotated[melt_split_anno_short_tss500_short_annotated$name2 %in% pois_others,]
probes_methylation_others2 <- melt_split_anno_short_tss500_short_annotated[melt_split_anno_short_tss500_short_annotated$name2 %in% pois_others2,]

# get methylation of these probes.
melt_df_beta_pois_esc <- merge(x=melt_df_beta, y=unique(probes_methylation_esc[,-1]), by = c("probeID"))
melt_df_beta_pois_inacs <- merge(x=melt_df_beta, y=unique(probes_methylation_inacs[,-1]), by = c("probeID"))
melt_df_beta_others <- merge(x=melt_df_beta, y=unique(probes_methylation_others[,-1]), by = c("probeID"))
melt_df_beta_others2 <- merge(x=melt_df_beta, y=unique(probes_methylation_others2[,-1]), by = c("probeID"))

# calculate median methylation.
melt_df_beta_pois_esc_median <- melt_df_beta_pois_esc %>% dplyr::group_by(probeID, name2, cell) %>% dplyr::summarise(medianz = median(value))
melt_df_beta_pois_inacs_median <- melt_df_beta_pois_inacs %>% dplyr::group_by(probeID, name2, cell) %>% dplyr::summarise(medianz = median(value))
melt_df_beta_others_median <- melt_df_beta_others %>% dplyr::group_by(probeID, name2, cell) %>% dplyr::summarise(medianz = median(value))
melt_df_beta_others2_median <- melt_df_beta_others2 %>% dplyr::group_by(probeID, name2, cell) %>% dplyr::summarise(medianz = median(value))

# plot
meth_esc_summary <-
  ggplot(melt_df_beta_pois_esc, aes(x=factor(cell, levels = c("ETP", "DN", "DPearly", "DPlate", "CD4SP" ,"CD8SP")), y=value)) +
  geom_boxplot() +
  labs(x="", y="methylaion (median beta)") + 
  scale_y_continuous(limits = c(0, 1))+
  theme_AL_box_rotX(legend.position = "none")+
  geom_hline(yintercept = c(0.25,0.75), lty=2)

meth_inac_summary <-
  ggplot(melt_df_beta_pois_inacs, aes(x=factor(cell, levels = c("ETP", "DN", "DPearly", "DPlate", "CD4SP" ,"CD8SP")), y=value)) +
  geom_boxplot() +
  labs(x="", y="methylaion (median beta)") + 
  scale_y_continuous(limits = c(0, 1))+
  theme_AL_box_rotX(legend.position = "none")+
  geom_hline(yintercept = c(0.25,0.75), lty=2)

methhh_inacs <- 
  ggplot(melt_df_beta_pois_inacs_median, aes(x=factor(cell, levels = c("ETP", "DN", "DPearly", "DPlate", "CD4SP" ,"CD8SP")), y=medianz, col = probeID)) + 
  geom_point() +
  geom_line(aes(group=probeID)) +
  facet_wrap(~factor(name2, levels = unique(annotation_gene_table[annotation_gene_table$chrom == "chrX",]$name2))) + 
  labs(x="", y="methylaion (median beta)") + 
  scale_y_continuous(limits = c(0, NA))+
  theme_AL_box_rotX(legend.position = "none")+
  geom_hline(yintercept = c(0.5), lty=2)


methhh_other2 <- 
  ggplot(melt_df_beta_others2_median, aes(x=factor(cell, levels = c("ETP", "DN", "DPearly", "DPlate", "CD4SP" ,"CD8SP")), y=medianz, col = probeID)) + 
  geom_point() +
  geom_line(aes(group=probeID)) +
  facet_wrap(~factor(name2, levels=c("ITM2A", "TLR7", "CD40LG")), ncol = 1) + 
  labs(x="", y="methylation (median beta)") + 
  scale_y_continuous(limits = c(0, 1))+
  theme_AL_box_rotX(legend.position = "none")+
  geom_hline(yintercept = c(0.5), lty=2)


# plot
ggsave2(filename = "04_plots/Figure_6c.pdf", width = 4, height = 6,
plot_grid(ncol=2, meth_esc_summary, meth_inac_summary))

ggsave2(filename = "04_plots/Suppl_Figure_7d.pdf", height = 14, width = 16,
        plot_grid(methhh_inacs))

ggsave2(filename = "04_plots/Figure_6d.pdf",
        plot_grid(ncol=3, methhh_other2))

# Export source data
melt_df_beta_pois_inacs$tag <- "Inactive"
melt_df_beta_pois_esc$tag <- "Escape"
write.table(rbind(melt_df_beta_pois_inacs, melt_df_beta_pois_esc), file = "06_source_data/Fig_6c.tsv", quote = F, row.names = F, sep = "\t")
write.table(melt_df_beta_others2_median, file = "06_source_data/Fig_6d.tsv", quote = F, row.names = F, sep = "\t")
write.table(melt_df_beta_pois_inacs_median, file = "06_source_data/Suppl_Fig_7d.tsv", quote = F, row.names = F, sep = "\t")
