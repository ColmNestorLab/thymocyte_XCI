library(ggalluvial)
library(data.table)
library(dplyr)
library(cowplot)
library(ggrepel)
library(rtracklayer)
library(GenomicRanges)

source("03_r_scripts/plot_parameters.R")

# starting point
# read in the sleuth df
df_sleuth <- subset(fread("Supplementary Table 5"), seqnames %in% c("chrX"))

df_sleuth_short <- dplyr::select(df_sleuth, target_id, b, qval)

# load in our ASE data df, this is the starting point of the new classification. In this df I will add/change genes as escapees that Colm has eyeballed as low in methylation and that are sex-biased
df_ase <- subset(fread("Supplementary Table 8"), contig == "chrX" & sample == "F3" & keep == T)

# create classification df.
our_classification <- unique(dplyr::select(df_ase, gene, effectSize))
our_classification$category <- ifelse(our_classification$effectSize > 0.4, yes = "Inactive", no = "Escape")
our_classification <- unique(dplyr::select(our_classification, -effectSize))
our_classification <-our_classification[!grepl(our_classification$gene, pattern = ";"),]

ASE_classification <- subset(our_classification, !(gene == "PRKX" & category == "Inactive") & !(gene == "PUDP" & category == "Inactive") & !(gene == "TXLNG" & category == "Inactive"))

# add
smash <- merge(df_sleuth_short, ASE_classification, by.x="target_id", by.y = "gene", all.x = T)

# read in Colm's eyeball results
pois <- fread("97_indexes_annotations/Xlinked_promoter_methylation.csv", header = F)

pois$methylation_curated <- "Low"

smash1 <- merge(smash, pois, by.x = "target_id", by.y = "V1", all.x =T ) 

################# Read in and prepare data ####################
# Load libraries

annotation_gene_table <- fread("97_indexes_annotations/UCSC_refGene.tsv")

annotation_gene_table_short <- subset(dplyr::select(annotation_gene_table, chrom, txStart, txEnd, name2, name, strand), chrom %in% c(paste0("chr", 1:22), "chrX"))

annotation_gene_table_short <- annotation_gene_table_short[!annotation_gene_table_short$name %in% c("NM_001320753", "NM_000351", "NM_001320754"),]

annotation_gene_table_short$start_fixed <- ifelse(annotation_gene_table_short$strand == "-", yes = annotation_gene_table_short$txEnd, no  = annotation_gene_table_short$txStart)

annotation_gene_table_short_tss1500 <- annotation_gene_table_short

anno <- data.frame(IlluminaHumanMethylationEPICanno.ilm10b5.hg38::Other)
anno$probeID <- rownames(anno)

anno_short <- dplyr::select(anno[!is.na(anno$CHR_hg38),], CHR_hg38, Start_hg38,  End_hg38, probeID)

anno_gr <- makeGRangesFromDataFrame(anno_short, ignore.strand=T, seqnames.field=c("CHR_hg38"), start.field=c("Start_hg38"),end.field=c("End_hg38"), keep.extra.columns = T)

annotation_gene_table_short_tss1500$first_tss1500 <- ifelse(annotation_gene_table_short_tss1500$strand == "+", yes = annotation_gene_table_short_tss1500$start_fixed - 500, no = annotation_gene_table_short_tss1500$start_fixed + 500)
annotation_gene_table_short_tss1500$second_tss1500 <- annotation_gene_table_short_tss1500$start_fixed + 1

annotation_gene_table_short_tss1500 <- transform(annotation_gene_table_short_tss1500, start_tss1500 = pmin(first_tss1500, second_tss1500))
annotation_gene_table_short_tss1500 <- transform(annotation_gene_table_short_tss1500, end_tss1500 = pmax(first_tss1500, second_tss1500))

annotation_gene_table_short_tss1500 <- dplyr::select(annotation_gene_table_short_tss1500, chrom, name2 ,name, start_tss1500, end_tss1500)

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

pois <- smash1$target_id


# read in EPIC data
df_beta <- data.frame(fread("Supplementary Table 11", header = T))

# Read in metadata
metadata <- fread("96_metadata/EPIC_thymocyte_meta.csv")

# melt
melt_df_beta <- melt(df_beta)

melt_df_beta$cell <- "ETP"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "DN"),]$cell <- "DN"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "DPl"),]$cell <- "DPearly"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "DPe"),]$cell <- "DPlate"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "CD4"),]$cell <- "CD4SP"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "CD8"),]$cell <- "CD8SP"

melt_df_beta$sample <- "F1"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "F3"),]$sample <- "F3"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "F4"),]$sample <- "F4"

WOW_did_it_really_work[WOW_did_it_really_work$name2 == "SEPTIN6",]$name2 <- "SEPT6"

probes_methylation <- WOW_did_it_really_work[WOW_did_it_really_work$name2 %in% pois,]



melt_df_beta_pois <- merge(x=melt_df_beta, y=unique(probes_methylation[,-1]), by = c("probeID"))

melt_df_beta_pois_test <- melt_df_beta_pois %>% dplyr::group_by(name2) %>% dplyr::summarise(meanz = mean(value))

melt_df_beta_pois_test$meth_cat <- "00000"
melt_df_beta_pois_test[melt_df_beta_pois_test$meanz <= 0.25,]$meth_cat <- "0.25"
melt_df_beta_pois_test[melt_df_beta_pois_test$meanz > 0.25 & melt_df_beta_pois_test$meanz < 0.75,]$meth_cat <- "0.26-0.75"
melt_df_beta_pois_test[melt_df_beta_pois_test$meanz >= 0.75,]$meth_cat <- "0.76"


smash1$sexB <- ifelse(smash1$b >= 0.2 & smash1$qval < 0.05, yes = "sexB", no = "non_sexB")

smash2 <- merge(smash1, melt_df_beta_pois_test, by.x = "target_id", by.y = "name2", all.x=T)


df_ase <- subset(fread("Supplementary Table 8"), contig == "chrX" & sample == "F3" & keep == T)

# create classification df.
our_classification <- unique(dplyr::select(df_ase, gene, effectSize, celltype))
our_classification_meanz <- our_classification %>% dplyr::group_by(gene) %>% dplyr::summarise(mean_ase = mean(effectSize))

smash2_short <- unique(merge(dplyr::select(smash2[!is.na(smash2$meth_cat),], target_id, b, meth_cat),our_classification_meanz, by.x = "target_id", by.y = "gene" , all.x =T))

smash2_short[smash2_short$b > 1,]$b <- 1
smash2_short[smash2_short$b < -1,]$b <- -1

smash2_short_filt <- smash2_short

smash2_short_filt$stackz <- ifelse(smash2_short_filt$mean_ase > 0.4, yes = "Inactive", no ="Escape")

smash2_short_order <- smash2_short_filt[order(-meth_cat,b ), ]

PAR <- rbind(fread("97_indexes_annotations/group-715_PAR1.csv"), fread("97_indexes_annotations/group-716_PAR2.csv"))
esc <- fread("97_indexes_annotations/landscape.Suppl.Table.13.csv")
esc[esc$Gene_name == "6-Sep",]$Gene_name <- "SEPT6"
esc <- esc[esc$Gene_ID != "ENSG00000241489.3",]
esc <- esc[,-1]

chrx_annotation <- merge(esc, PAR, by.x = "Gene_name" , by.y = "Approved symbol" , all = T)[-1,]
chrx_annotation$category <- chrx_annotation$Reported_XCI_status
chrx_annotation[chrx_annotation$PAR == "PAR1" | chrx_annotation$PAR == "PAR2",]$category <- "PAR"

chrx_annotation[chrx_annotation$Gene_name == "KAL1",]$Gene_name <- "ANOS1"
chrx_annotation[chrx_annotation$Gene_name == "CXorf36",]$Gene_name <- "DIPK2B"
chrx_annotation[chrx_annotation$Gene_name == "HDHD1",]$Gene_name <- "PUDP"
chrx_annotation[chrx_annotation$Gene_name == "RGAG4",]$Gene_name <- "RTL5"


smash2_short_melt <- melt(smash2_short_filt, id.vars = c("target_id", "stackz", "meth_cat"), measure.vars = c("b"))

smash2_short_melt_sans_PAR <- smash2_short_melt[!(smash2_short_melt$target_id %in% PAR$`Approved symbol`),]

ggsave2(filename = "04_plots/Figure_3C.pdf", width = 18, height = 6,
        ggplot(smash2_short_melt_sans_PAR, aes(y=factor(meth_cat, levels = rev(c("0.25", "0.26-0.75","0.76"))), x=factor(target_id, levels = rev(unique(smash2_short_order$target_id)))) ) + 
          geom_tile(aes(fill=value)) + 
          scale_fill_gradientn(colours=c("red","grey","black"), limits =c(-1,1)) + 
          facet_wrap(stackz~., ncol = 1, scales = "free_x") + 
          theme_AL_box_rotX() + 
          theme(axis.text.x =element_text(size=4, face = "italic")) + labs(x="",y="")+
          geom_vline(xintercept = c("F8", "UBA1", "TIMM17B", "OPHN1", "ARR3", "TREX2", "NLGN4X", "NAP1L3", "ITM2A", "CD40LG")))

#smash2_short_melt_sans_PAR_NA_Stacks <- smash2_short_melt[!(smash2_short_melt$target_id %in% PAR$`Approved symbol`) & is.na(smash2_short_melt$stackz),]
#smash2_short_melt_sans_PAR_inactive_Stacks <- smash2_short_melt[!(smash2_short_melt$target_id %in% PAR$`Approved symbol`) & smash2_short_melt$stackz == "Inactive",]
#smash2_short_melt_sans_PAR_escape_Stacks <- smash2_short_melt[!(smash2_short_melt$target_id %in% PAR$`Approved symbol`) & smash2_short_melt$stackz == "Escape",]
#selectzzz <- smash2_short_melt[smash2_short_melt$stackz == "Inactive" | (is.na(smash2_short_melt$stackz) & smash2_short_melt$meth_cat != "0.25" & smash2_short_melt$value < 0.2),]
#inactivezzz <- esc[esc$Gene_name %in% selectzzz$target_id,]
#inactivezzz_fin <- inactivezzz[inactivezzz$Reported_XCI_status == "Escape",]


Figure3C_suppl_table <- smash2_short_melt
colnames(Figure3C_suppl_table) <- c("gene", "ASE_F3", "Methylation_category", "b", "log2FC") 

# add classificaiton column, exporting which XCI call we made per gene.
Figure3C_suppl_table$category <- "LOL"
Figure3C_suppl_table[!is.na(Figure3C_suppl_table$ASE_F3),]$category <- Figure3C_suppl_table[!is.na(Figure3C_suppl_table$ASE_F3),]$ASE_F3
Figure3C_suppl_table[Figure3C_suppl_table$gene == "XIST"]$category <- "Escape"
Figure3C_suppl_table[is.na(Figure3C_suppl_table$ASE_F3) & Figure3C_suppl_table$Methylation_category == 0.25 & Figure3C_suppl_table$log2FC > 0.2,]$category <- "high_confidence_escapee"
Figure3C_suppl_table[is.na(Figure3C_suppl_table$ASE_F3) & Figure3C_suppl_table$Methylation_category == 0.25 & Figure3C_suppl_table$log2FC < 0.2 & Figure3C_suppl_table$log2FC > 0,]$category <- "low_confidence_escapee"

Figure3C_suppl_table[is.na(Figure3C_suppl_table$ASE_F3) & Figure3C_suppl_table$Methylation_category != 0.25 & Figure3C_suppl_table$log2FC < 0.2,]$category <- "Inactive"

Figure3C_suppl_table[is.na(Figure3C_suppl_table$ASE_F3) & 
                     (
                         (Figure3C_suppl_table$Methylation_category == 0.25 & Figure3C_suppl_table$log2FC < 0) 
                         |
                         (Figure3C_suppl_table$Methylation_category != 0.25 & Figure3C_suppl_table$log2FC > 0.2)
                     ) ,]$category <- "Unclassified"


#Genes were assigned as high-confidence escapees if either a) allele-specific expression was lower or equal to 0.4 or b) if the gene’s mean methylation -500 TSS value was below 0.25 and sex-biased expression above 0.2 log2FC. Genes with a mean methylation -500 TSS value below 0.25 and sex-biased expression above 0 but below 0.2 log2FC were assigned low-confidence escapees. Genes were assigned as inactive genes if a) allele-specific expression was higher than 0.4 or b) if the genes mean methylation -500 TSS value was above 0.25 and its sex-biased expression was below 0.2 log2FC. Genes that fell outside of these criteria were assigned “unclassified” status. This category thus includes genes that lack ASE coverage and genes with a methylation value below 0.25 and a sex-bias below 0 log2FC and genes with a methylation value above 0.25 but with a sex-bias above 0 log2FC. 


write.table(dplyr::select(Figure3C_suppl_table, -b), file = "05_reports/Suppl_table_12.tsv", quote = F, col.names = T, row.names = F, sep = "\t")



class_for_alluvium <- Figure3C_suppl_table

#class_for_alluvium[class_for_alluvium$meth_cat == "0.25" & class_for_alluvium$value >= 0.2 & is.na(class_for_alluvium$stackz),]$stackz <- "Escape"
#class_for_alluvium[class_for_alluvium$meth_cat != "0.25" & class_for_alluvium$value < 0.2 & is.na(class_for_alluvium$stackz),]$stackz <- "Inactive"
#class_for_alluvium[is.na(class_for_alluvium$stackz),]$stackz <- "Unknown"


df_TPM <- fread("Supplementary Table 3")
df_TPM_filt <- melt(df_TPM)


df_TPM_filt <- df_TPM_filt %>% dplyr::group_by(gene) %>% dplyr::summarise(meanz = mean(value))





# select the original Tuki df, this is the original classification of escapees
df_classification <- dplyr::select(chrx_annotation[chrx_annotation$Gene_name %in% class_for_alluvium$gene,], Gene_name, category)
colnames(df_classification) <- c("gene", "category")
df_classification$axes <- "Tukiainen"


df_gylemo <- unique(dplyr::select(class_for_alluvium, gene, category))
colnames(df_gylemo) <- c("gene", "category")
df_gylemo$axes <- "Gylemo"


novelz <- data.frame(gene = setdiff(df_gylemo$gene, df_classification$gene),
                             category = "novel")
novelz$axes <- "Unannotated"


keklol <- df_gylemo[!df_gylemo$gene %in% df_classification$gene,]
df_novelz <- rbind(novelz, keklol)


is_lodes_form(df_novelz, key = "axes", value = "category", id = "gene")

df_novelz$axes <- factor(df_novelz$axes, levels = c("Unannotated", "Gylemo"))

df_novelz <- df_novelz[!df_novelz$gene %in% PAR$`Approved symbol`,]



smash_fin <- rbind(df_classification, df_gylemo)

smash_fin$axes <- factor(smash_fin$axes, levels = c("Tukiainen", "Gylemo"))

smash_fin <- smash_fin[!smash_fin$gene %in% PAR$`Approved symbol` & !smash_fin$gene %in% df_novelz$gene ,]

is_lodes_form(smash_fin, key = "axes", value = "category", id = "gene")


  


smashf_fin_filt <- smash_fin[smash_fin$gene %in% df_TPM_filt[df_TPM_filt$meanz > 1,]$gene,]

smashf_fin_filt$category <- factor(smashf_fin_filt$category, levels = c("Escape","high_confidence_escapee", "low_confidence_escapee", "Inactive", "Variable", "Unknown", "Unclassified"))

smashf_fin_filt %>% dplyr::group_by(category, axes) %>% dplyr::count()

smashf_fin_filt



ggsave2(filename = "04_plots/Figure_3D.pdf",
        plot_grid(ggplot(smashf_fin_filt, aes(alluvium = gene, x = axes, stratum = category)) + 
                    geom_flow(aes(fill=category)) +
                    geom_stratum(aes(fill=category))  + 
                    geom_text(stat="count", aes(label=after_stat(count), group = category), position = "stack", vjust = 2)+
                    theme_AL_simple(legend.title = element_blank())+
                    scale_fill_manual(values = c("Escape" = "red", "high_confidence_escapee" = "#FF474C", "low_confidence_escapee" = "pink", "Inactive" = "#869a9a", "Variable" = "purple","Unknown" = "blue", "Unclassified" = "blue"))
                  )
)



figure3D_suppl_table <- smashf_fin_filt

figure3D_suppl_table[figure3D_suppl_table$axes == "Tukiainen",]$axes <- "previous assessment"
figure3D_suppl_table[figure3D_suppl_table$axes == "Gylemo",]$axes <- "thymocytes"

figure3D_suppl_table <- merge(figure3D_suppl_table[figure3D_suppl_table$axes == "thymocytes",c("gene","category")], figure3D_suppl_table[figure3D_suppl_table$axes == "previous assessment",c("gene","category")], by = "gene")

colnames(figure3D_suppl_table) <- c("gene", "thymocytes_classification", "previous_assessment")

figure3D_suppl_table$comparison <- ifelse(figure3D_suppl_table$thymocytes_classification == figure3D_suppl_table$previous_assessment, yes = "same", no = "differs")

figure3D_suppl_table[figure3D_suppl_table$thymocytes_classification == "high_confidence_escapee" & figure3D_suppl_table$previous_assessment == "Escape"]$comparison <- "same"
figure3D_suppl_table[figure3D_suppl_table$thymocytes_classification == "low_confidence_escapee" & figure3D_suppl_table$previous_assessment == "Escape"]$comparison <- "same"

figure3D_suppl_table %>% dplyr::group_by(comparison) %>% dplyr::count()

write.table(dplyr::select(figure3D_suppl_table, gene, previous_assessment, thymocytes_classification, comparison), file = "05_reports/Suppl_table_13.tsv", quote = F, col.names = T, row.names = F, sep = "\t")


