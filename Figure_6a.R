# Load libraries
library(data.table)
library(pvclust)

source("plot_parameters.R")


# Read in data (created and saved above)
df_beta <- data.frame(fread("02_tidy_data/Suppl_tables/Supplementary_Data_methylation_values.tsv", header = T))

# Read in metadata
metadata <- fread("96_metadata/EPIC_thymocyte_meta.csv")

# melt
melt_df_beta <- melt(df_beta)

# calculate SD
variable_genes <- melt_df_beta %>% dplyr::group_by(probeID) %>% dplyr::summarise(sdz = sd(value))

# get and keep the 1000 most variable probes (based on highest SD).
variable_genes_top1000 <- head(variable_genes[order(variable_genes$sdz, decreasing = T),], 1000)
df_pvclust <- df_beta[df_beta$probeID %in% variable_genes_top1000$probeID,]

# add probeIDs as rownames, then remove the column.
rownames(df_pvclust) <- df_pvclust$probeID
df_pvclust$probeID <- NULL

# make pvclust object, using the top 1000 most variable probes
results <- pvclust(df_pvclust, nboot = 10000, method = "euclidean", 
                   method.hclust = "complete", parallel = F)

# Plot
pdf("04_plots/Figure_6a.pdf", width = 10)
plot(results)
pvrect(results, alpha=0.95)
dev.off()

# export source data
df_pvclust_export <- df_pvclust
df_pvclust_export$probeID <- rownames(df_pvclust_export)
write.table(df_pvclust_export, file = "06_source_data/Fig_6a.tsv", quote = F, row.names = F, sep = "\t")



