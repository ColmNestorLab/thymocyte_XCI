################# Read in and prepare data ####################
# Load libraries
library(data.table)
library(pvclust)

source("03_r_scripts/plot_parameters.R")


# Read in data (created and saved above)
df_beta <- data.frame(fread("Supplementary Table 11", header = T))

# Read in metadata
metadata <- fread("96_metadata/EPIC_thymocyte_meta.csv")

# melt
melt_df_beta <- melt(df_beta)

# calculate SD
variable_genes <- melt_df_beta %>% dplyr::group_by(probeID) %>% dplyr::summarise(sdz = sd(value))

# get the 1000 most variable genes
variable_genes_top1000 <- head(variable_genes[order(variable_genes$sdz, decreasing = T),], 1000)

df_pvclust <- df_beta[df_beta$probeID %in% variable_genes_top1000$probeID,]
rownames(df_pvclust) <- df_pvclust$probeID
df_pvclust$probeID <- NULL

# make pvclust object, including on the top 1000 most variable genes
results <- pvclust(df_pvclust, nboot = 10000, method = "euclidean", 
                   method.hclust = "complete", parallel = F)

pdf("04_plots/Figure_3A.pdf", width = 10)
plot(results)
pvrect(results, alpha=0.95)
dev.off()
