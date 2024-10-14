library(ggalluvial)
library(data.table)
library(dplyr)
library(cowplot)
library(ggrepel)
library(rtracklayer)
library(GenomicRanges)
library(readxl)

source("plot_parameters.R")

# read in the sleuth df
df_sleuth <- data.table(read_excel("02_tidy_data/Suppl_tables/Supplementary Tables.xlsx", sheet = 5, skip = 1))

# keep only columns we want
df_sleuth_short <- dplyr::select(df_sleuth, gene, log2FC, qval)

# load in our ASE data df, this is the starting point of the new classification. In this df I will add genes as escapees that are low in methylation and that are sex-biased.
df_ase <- subset(data.table(read_excel("02_tidy_data/Suppl_tables/Supplementary Tables.xlsx", sheet = 10, skip = 1, col_types = c("text", "numeric", "text", "text", "text", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "logical", "logical"))), contig == "chrX" & sample == "F3" & keep == T)

# create classification df.
our_classification <- unique(dplyr::select(df_ase, gene, effectSize))
our_classification$category <- ifelse(our_classification$effectSize > 0.4, yes = "Inactive", no = "Escape")
our_classification <- unique(dplyr::select(our_classification, -effectSize))
our_classification <-our_classification[!grepl(our_classification$gene, pattern = ";"),]

# minor fixes
ASE_classification <- subset(our_classification, !(gene == "PRKX" & category == "Inactive") & !(gene == "PUDP" & category == "Inactive") & !(gene == "TXLNG" & category == "Inactive"))

# add annotation
smash <- merge(df_sleuth_short, ASE_classification, by.x="gene", by.y = "gene", all.x = T)

# read in Colm's manual curation and add to the df.
pois <- fread("97_indexes_annotations/Xlinked_promoter_methylation.csv", header = F)
pois$methylation_curated <- "Low"
smash1 <- merge(smash, pois, by.x = "gene", by.y = "V1", all.x =T ) 

# store genes of interest
pois <- smash1$gene


# methylation annotation
# read in annotation
annotation_gene_table <- fread("97_indexes_annotations/UCSC_refGene.tsv")

# keep only columns we need
annotation_gene_table_short <- subset(dplyr::select(annotation_gene_table, chrom, txStart, txEnd, name2, name, strand), chrom %in% c(paste0("chr", 1:22), "chrX"))

# remove strange transcripts
annotation_gene_table_short <- annotation_gene_table_short[!annotation_gene_table_short$name %in% c("NM_001320753", "NM_000351", "NM_001320754"),]

# fix start, making it dependent on which strand the gene is located.
annotation_gene_table_short$start_fixed <- ifelse(annotation_gene_table_short$strand == "-", yes = annotation_gene_table_short$txEnd, no  = annotation_gene_table_short$txStart)

# Read in probe methylation
anno <- data.frame(IlluminaHumanMethylationEPICanno.ilm10b5.hg38::Other)
anno$probeID <- rownames(anno)

# keep only columns we need
anno_short <- dplyr::select(anno[!is.na(anno$CHR_hg38),], CHR_hg38, Start_hg38,  End_hg38, probeID)

# make GRanges object
anno_gr <- makeGRangesFromDataFrame(anno_short, ignore.strand=T, seqnames.field=c("CHR_hg38"), start.field=c("Start_hg38"),end.field=c("End_hg38"), keep.extra.columns = T)

# make new df, adding 500 to the TSS.
annotation_gene_table_short_tss500 <- annotation_gene_table_short
annotation_gene_table_short_tss500$first_tss500 <- ifelse(annotation_gene_table_short_tss500$strand == "+", yes = annotation_gene_table_short_tss500$start_fixed - 500, no = annotation_gene_table_short_tss500$start_fixed + 500)

# add this column to make the orientation of the gene correct.
annotation_gene_table_short_tss500$second_tss500 <- annotation_gene_table_short_tss500$start_fixed + 1
annotation_gene_table_short_tss500 <- transform(annotation_gene_table_short_tss500, start_tss500 = pmin(first_tss500, second_tss500))
annotation_gene_table_short_tss500 <- transform(annotation_gene_table_short_tss500, end_tss500 = pmax(first_tss500, second_tss500))

# keep only columns we need, rename them.
annotation_gene_table_short_tss500 <- dplyr::select(annotation_gene_table_short_tss500, chrom, name2 ,name, start_tss500, end_tss500)
colnames(annotation_gene_table_short_tss500) <- c("seqnames","gene_id","transcript_id","start_tss500","end_tss500")

# make GRanges object
tss500_gr <- makeGRangesFromDataFrame(annotation_gene_table_short_tss500, seqnames.field=c("seqnames"), start.field=c("start_tss500"), end.field=c("end_tss500"), keep.extra.columns = T)

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
split_anno_short_tss500_short <- anno_short_tss500_short %>% tidyr::separate(transcript_id, sep=";" , into = c("a","log2FC","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u"))

# see if it worked.
melt_split_anno_short_tss500_short <- dplyr::select(subset(melt(split_anno_short_tss500_short, id.vars = "probeID"), !is.na(value)), -variable)

# merge annotation with the overlap df.
melt_split_anno_short_tss500_short_annotated <- unique(merge(dplyr::select(annotation_gene_table, name, name2), melt_split_anno_short_tss500_short, by.x = "name", by.y = "value", allow.cartesian=TRUE))




# read in EPIC data
df_beta <- data.frame(fread("02_tidy_data/Suppl_tables/Supplementary_Data_methylation_values.tsv", header = T))

# Read in metadata
metadata <- fread("96_metadata/EPIC_thymocyte_meta.csv")

# melt
melt_df_beta <- melt(df_beta)

# add celltype
melt_df_beta$cell <- "ETP"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "DN"),]$cell <- "DN"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "DPl"),]$cell <- "DPearly"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "DPe"),]$cell <- "DPlate"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "CD4"),]$cell <- "CD4SP"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "CD8"),]$cell <- "CD8SP"

# add sample
melt_df_beta$sample <- "F1"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "F3"),]$sample <- "F3"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "F4"),]$sample <- "F4"

# fix name
melt_split_anno_short_tss500_short_annotated[melt_split_anno_short_tss500_short_annotated$name2 == "SEPTIN6",]$name2 <- "SEPT6"

# get probes that are within the tss-500 of genes that we are interested in (pois).
probes_methylation <- melt_split_anno_short_tss500_short_annotated[melt_split_anno_short_tss500_short_annotated$name2 %in% pois,]

# get methylation of these probes. 
melt_df_beta_pois <- merge(x=melt_df_beta, y=unique(probes_methylation[,-1]), by = c("probeID"))

# get means
melt_df_beta_pois_test <- melt_df_beta_pois %>% dplyr::group_by(name2) %>% dplyr::summarise(meanz = mean(value))

# Set the methylation category of these regions per gene
melt_df_beta_pois_test$meth_cat <- "00000"
melt_df_beta_pois_test[melt_df_beta_pois_test$meanz <= 0.25,]$meth_cat <- "0.25"
melt_df_beta_pois_test[melt_df_beta_pois_test$meanz > 0.25 & melt_df_beta_pois_test$meanz < 0.75,]$meth_cat <- "0.26-0.75"
melt_df_beta_pois_test[melt_df_beta_pois_test$meanz >= 0.75,]$meth_cat <- "0.76"

# add tag whether or not its sex-biased.
smash1$sexB <- ifelse(smash1$log2FC >= 0.2 & smash1$qval < 0.05, yes = "sexB", no = "non_sexB")

# merge the data frames
smash2 <- merge(smash1, melt_df_beta_pois_test, by.x = "gene", by.y = "name2", all.x=T)


# get mean ASE.
mean_ASE <- unique(dplyr::select(df_ase, gene, effectSize, celltype)) %>% dplyr::group_by(gene) %>% dplyr::summarise(mean_ase = mean(effectSize))

# merge
smash2_short <- unique(merge(dplyr::select(smash2[!is.na(smash2$meth_cat),], gene, log2FC, meth_cat),mean_ASE, by.x = "gene", by.y = "gene" , all.x =T))

# cap log2FC values to 1 and -1 for easier plotting.
smash2_short[smash2_short$log2FC > 1,]$log2FC <- 1
smash2_short[smash2_short$log2FC < -1,]$log2FC <- -1

# add tag
smash2_short_filt <- smash2_short
smash2_short_filt$stackz <- ifelse(smash2_short_filt$mean_ase > 0.4, yes = "Inactive", no ="Escape")

# plot order
smash2_short_order <- smash2_short_filt[order(-meth_cat,log2FC ), ]

# read in PAR and Tukiainen annotation
PAR <- rbind(fread("97_indexes_annotations/group-715_PAR1.csv"), fread("97_indexes_annotations/group-716_PAR2.csv"))
esc <- fread("97_indexes_annotations/landscape.Suppl.Table.13.csv")
esc[esc$Gene_name == "6-Sep",]$Gene_name <- "SEPT6"
esc <- esc[esc$Gene_ID != "ENSG00000241489.3",]
esc <- esc[,-1]

# merge
chrx_annotation <- merge(esc, PAR, by.x = "Gene_name" , by.y = "Approved symbol" , all = T)[-1,]
chrx_annotation$category <- chrx_annotation$Reported_XCI_status
chrx_annotation[chrx_annotation$PAR == "PAR1" | chrx_annotation$PAR == "PAR2",]$category <- "PAR"

# change known gene names
chrx_annotation[chrx_annotation$Gene_name == "KAL1",]$Gene_name <- "ANOS1"
chrx_annotation[chrx_annotation$Gene_name == "CXorf36",]$Gene_name <- "DIPK2B"
chrx_annotation[chrx_annotation$Gene_name == "HDHD1",]$Gene_name <- "PUDP"
chrx_annotation[chrx_annotation$Gene_name == "RGAG4",]$Gene_name <- "RTL5"

# melt
smash2_short_melt <- melt(smash2_short_filt, id.vars = c("gene", "stackz", "meth_cat"), measure.vars = c("log2FC"))

# remove PAR
smash2_short_melt_sans_PAR <- smash2_short_melt[!(smash2_short_melt$gene %in% PAR$`Approved symbol`),]

# Plot
ggsave2(filename = "04_plots/Figure_7a.pdf", width = 18, height = 6,
        ggplot(smash2_short_melt_sans_PAR, aes(y=factor(meth_cat, levels = rev(c("0.25", "0.26-0.75","0.76"))), x=factor(gene, levels = rev(unique(smash2_short_order$gene)))) ) + 
          geom_tile(aes(fill=value)) + 
          scale_fill_gradientn(colours=c("red","grey","black"), limits =c(-1,1)) + 
          facet_wrap(stackz~., ncol = 1, scales = "free_x") + 
          theme_AL_box_rotX() + 
          theme(axis.text.x =element_text(size=4, face = "italic")) + labs(x="",y="")+
          geom_vline(xintercept = c("F8", "UBA1", "TIMM17B", "OPHN1", "ARR3", "TREX2", "NLGN4X", "NAP1L3", "ITM2A", "CD40LG")))

# make suppl table 13
suppl_table_13 <- smash2_short_melt
colnames(suppl_table_13) <- c("gene", "ASE_F3", "Methylation_category", "log2FC", "log2FC") 

# add classification column, exporting which XCI call we made per gene.
suppl_table_13$category <- "LOL"
suppl_table_13[!is.na(suppl_table_13$ASE_F3),]$category <- suppl_table_13[!is.na(suppl_table_13$ASE_F3),]$ASE_F3
suppl_table_13[suppl_table_13$gene == "XIST"]$category <- "Escape"
suppl_table_13[is.na(suppl_table_13$ASE_F3) & suppl_table_13$Methylation_category == 0.25 & suppl_table_13$log2FC > 0.2,]$category <- "high_confidence_escapee"
suppl_table_13[is.na(suppl_table_13$ASE_F3) & suppl_table_13$Methylation_category == 0.25 & suppl_table_13$log2FC < 0.2 & suppl_table_13$log2FC > 0,]$category <- "low_confidence_escapee"

suppl_table_13[is.na(suppl_table_13$ASE_F3) & suppl_table_13$Methylation_category != 0.25 & suppl_table_13$log2FC < 0.2,]$category <- "Inactive"

suppl_table_13[is.na(suppl_table_13$ASE_F3) & 
                     (
                         (suppl_table_13$Methylation_category == 0.25 & suppl_table_13$log2FC < 0) 
                         |
                         (suppl_table_13$Methylation_category != 0.25 & suppl_table_13$log2FC > 0.2)
                     ) ,]$category <- "Unclassified"

#Genes were assigned as high-confidence escapees if either a) allele-specific expression was lower or equal to 0.4 or log2FC) if the gene’s mean methylation -500 TSS value was below 0.25 and sex-biased expression above 0.2 log2FC. Genes with a mean methylation -500 TSS value below 0.25 and sex-biased expression above 0 but below 0.2 log2FC were assigned low-confidence escapees. Genes were assigned as inactive genes if a) allele-specific expression was higher than 0.4 or log2FC) if the genes mean methylation -500 TSS value was above 0.25 and its sex-biased expression was below 0.2 log2FC. Genes that fell outside of these criteria were assigned “unclassified” status. This category thus includes genes that lack ASE coverage and genes with a methylation value below 0.25 and a sex-bias below 0 log2FC and genes with a methylation value above 0.25 but with a sex-bias above 0 log2FC. 


write.table(dplyr::select(suppl_table_13, -log2FC), file = "05_reports/Suppl_table_13.tsv", quote = F, col.names = T, row.names = F, sep = "\t")


# make alluvial plot
class_for_alluvium <- suppl_table_13

# read in thymocyte TPMs
df_TPM <- data.table(read_excel("02_tidy_data/Suppl_tables/Supplementary Tables.xlsx", sheet = 3, skip = 1, col_types = c("text",rep(paste0("numeric"), 46))))
df_TPM_filt <- melt(df_TPM, id.vars = "gene")

# calculate means
df_TPM_filt <- df_TPM_filt %>% dplyr::group_by(gene) %>% dplyr::summarise(meanz = mean(value))

# select the original Tuki df, this is the original classification of escapees
df_classification <- dplyr::select(chrx_annotation[chrx_annotation$Gene_name %in% class_for_alluvium$gene,], Gene_name, category)
colnames(df_classification) <- c("gene", "category")
df_classification$axes <- "Tukiainen"

# make classification of thymocytes
df_gylemo <- unique(dplyr::select(class_for_alluvium, gene, category))
colnames(df_gylemo) <- c("gene", "category")
df_gylemo$axes <- "Gylemo"

# find genes not covered in Tukiainens classification.
novelz <- data.frame(gene = setdiff(df_gylemo$gene, df_classification$gene),
                             category = "novel")
novelz$axes <- "Unannotated"

# remove genes in our df that is not found in their df, removing PAR genes.
keklol <- df_gylemo[!df_gylemo$gene %in% df_classification$gene,]
df_novelz <- rbind(novelz, keklol)
df_novelz <- df_novelz[!df_novelz$gene %in% PAR$`Approved symbol`,]

# check if format is OK
is_lodes_form(df_novelz, key = "axes", value = "category", id = "gene")

# level the factor
df_novelz$axes <- factor(df_novelz$axes, levels = c("Unannotated", "Gylemo"))


# merge
smash_fin <- rbind(df_classification, df_gylemo)

# level the factor
smash_fin$axes <- factor(smash_fin$axes, levels = c("Tukiainen", "Gylemo"))

# remove PAR genes and genes that are not covered in Tukiainens df
smash_fin <- smash_fin[!smash_fin$gene %in% PAR$`Approved symbol` & !smash_fin$gene %in% df_novelz$gene ,]

# check if format is OK
is_lodes_form(smash_fin, key = "axes", value = "category", id = "gene")


  

# remove genes that are not expressed
smashf_fin_filt <- smash_fin[smash_fin$gene %in% df_TPM_filt[df_TPM_filt$meanz > 1,]$gene,]

# level the factor
smashf_fin_filt$category <- factor(smashf_fin_filt$category, levels = c("Escape","high_confidence_escapee", "low_confidence_escapee", "Inactive", "Variable", "Unknown", "Unclassified"))

# count
smashf_fin_filt %>% dplyr::group_by(category, axes) %>% dplyr::count()



# plot
ggsave2(filename = "04_plots/Figure_7b.pdf",
        plot_grid(ggplot(smashf_fin_filt, aes(alluvium = gene, x = axes, stratum = category)) + 
                    geom_flow(aes(fill=category)) +
                    geom_stratum(aes(fill=category))  + 
                    geom_text(stat="count", aes(label=after_stat(count), group = category), position = "stack", vjust = 2)+
                    theme_AL_simple(legend.title = element_blank())+
                    scale_fill_manual(values = c("Escape" = "red", "high_confidence_escapee" = "#FF474C", "low_confidence_escapee" = "pink", "Inactive" = "#869a9a", "Variable" = "purple","Unknown" = "blue", "Unclassified" = "blue"))
                  )
)


# make suppl table 14
suppl_table_14 <- smashf_fin_filt

suppl_table_14[suppl_table_14$axes == "Tukiainen",]$axes <- "previous assessment"
suppl_table_14[suppl_table_14$axes == "Gylemo",]$axes <- "thymocytes"

suppl_table_14 <- merge(suppl_table_14[suppl_table_14$axes == "thymocytes",c("gene","category")], suppl_table_14[suppl_table_14$axes == "previous assessment",c("gene","category")], by = "gene")

colnames(suppl_table_14) <- c("gene", "thymocytes_classification", "previous_assessment")

suppl_table_14$comparison <- ifelse(suppl_table_14$thymocytes_classification == suppl_table_14$previous_assessment, yes = "same", no = "differs")

suppl_table_14[suppl_table_14$thymocytes_classification == "high_confidence_escapee" & suppl_table_14$previous_assessment == "Escape"]$comparison <- "same"
suppl_table_14[suppl_table_14$thymocytes_classification == "low_confidence_escapee" & suppl_table_14$previous_assessment == "Escape"]$comparison <- "same"

suppl_table_14 %>% dplyr::group_by(comparison) %>% dplyr::count()

write.table(dplyr::select(suppl_table_14, gene, previous_assessment, thymocytes_classification, comparison), file = "05_reports/Suppl_table_13.tsv", quote = F, col.names = T, row.names = F, sep = "\t")


