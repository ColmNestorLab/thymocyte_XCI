################# Read in and prepare data ####################
# Load libraries
library(data.table)
library(karyoploteR)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(ggplotify)
library(cowplot)
library(corrplot)

source("03_r_scripts/plot_parameters.R")

etp_bigwig.file <- "02_tidy_data/thymocyte/EPIC/bedgraphs_bigwigs/ETP.bw"
dn_bigwig.file <- "02_tidy_data/thymocyte/EPIC/bedgraphs_bigwigs/DN.bw"
dpe_bigwig.file <- "02_tidy_data/thymocyte/EPIC/bedgraphs_bigwigs/DPearly.bw"
dpl_bigwig.file <- "02_tidy_data/thymocyte/EPIC/bedgraphs_bigwigs/DPlate.bw"
cd4_bigwig.file <- "02_tidy_data/thymocyte/EPIC/bedgraphs_bigwigs/CD4SP.bw"
cd8_bigwig.file <- "02_tidy_data/thymocyte/EPIC/bedgraphs_bigwigs/CD8SP.bw"


#plotDefaultPlotParams(plot.type=1)
plot.params <- getDefaultPlotParams(plot.type=1)
plot.params$topmargin <- 2.5
plot.params$bottommargin <- 25
plot.params$data1outmargin <- 5



karyplotz <- as.ggplot(expression(
kp <- plotKaryotype(genome="hg38", chromosomes = "chrX", plot.type = 1, plot.params = plot.params),
kpArrows(kp, chr="chrX", x0=2784257, x1=2784257, y0=-0.4, y1=-0.35, cex=0.5),
kpText(kp, chr="chrX", x=2784257, y=-0.425, labels="2784257", cex=0.75),
kpAddBaseNumbers(kp),
kpPlotBigWig(kp, data=etp_bigwig.file, ymin = 0, ymax=1, r0=0, r1=0.15),
kpPlotBigWig(kp, data=dn_bigwig.file, ymax=1, r0=0.16, r1=0.31),
kpPlotBigWig(kp, data=dpe_bigwig.file, ymax=1, r0=0.32, r1=0.47),
kpPlotBigWig(kp, data=dpl_bigwig.file, ymax=1, r0=0.48, r1=0.63),
kpPlotBigWig(kp, data=cd4_bigwig.file, ymax=1, r0=0.64, r1=0.79),
kpPlotBigWig(kp, data=cd8_bigwig.file, ymax=1, r0=0.8, r1=0.95),

kpAxis(kp, data.panel=1, ymin = 0, ymax=1, r0=0, r1=0.15),
kpAxis(kp, data.panel=1, ymin = 0, ymax=1, r0=0.16, r1=0.31),
kpAxis(kp, data.panel=1, ymin = 0, ymax=1, r0=0.32, r1=0.47),
kpAxis(kp, data.panel=1, ymin = 0, ymax=1, r0=0.48, r1=0.63),
kpAxis(kp, data.panel=1, ymin = 0, ymax=1, r0=0.64, r1=0.79),
kpAxis(kp, data.panel=1, ymin = 0, ymax=1, r0=0.8, r1=0.95)
))

anno <- data.frame(IlluminaHumanMethylationEPICanno.ilm10b5.hg38::Other)
anno$probeID <- rownames(anno)
anno_short <- dplyr::select(anno[!is.na(anno$CHR_hg38),], CHR_hg38, Start_hg38,  End_hg38, probeID)
anno_short_chrx <- anno_short[anno_short$CHR_hg38 == "chrX",]


# Read in methylation
df_beta <- data.frame(fread("Supplementary Table 11", header = T))

# melt
melt_df_beta <- melt(df_beta)

melt_df_beta$celltype <- "ETP"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "DN"),]$celltype <- "DN"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "DPe"),]$celltype <- "DPearly"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "DPl"),]$celltype <- "DPlate"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "CD4"),]$celltype <- "CD4SP"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "CD8"),]$celltype <- "CD8SP"

melt_df_beta$sample_ID <- "F1"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "F3"),]$sample_ID <- "F3"
melt_df_beta[grepl(melt_df_beta$variable, pattern = "F4"),]$sample_ID <- "F4"

# merge keep only chrX
melt_df_beta_chrx <- melt_df_beta[melt_df_beta$probeID %in% anno_short_chrx$probeID,]

# get mean methylation per celltype
melt_df_beta_chrx_meanz <- melt_df_beta_chrx %>% dplyr::group_by(probeID, celltype) %>% dplyr::summarise(meanz = mean(value))

# make matrix
mat <- dcast(melt_df_beta_chrx_meanz, probeID ~ celltype, value.var = "meanz")
rownames(mat) <- mat$probeID
mat$probeID <- NULL

res <- cor(mat, method = "spearman")




ggsave2("04_plots/Suppl_Figure_5B_left.pdf", width = 24,
        plot_grid(karyplotz))

pdf("04_plots/Suppl_Figure_5B_right.pdf")
        corrplot.mixed(res, order = 'AOE')
dev.off()