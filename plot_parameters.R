# author of code: Antonio Lentini.
# modifies ggplots to look the same across plots.
require(ggplot2)
theme_AL_simple <- function (...) { 
  theme_bw(base_size=12, base_family="") %+replace% 
    theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(colour="black"),axis.line=element_line(colour="black"), ... )
}
theme_AL_simple_rotX <- function (...) { 
  theme_bw(base_size=12, base_family="") %+replace% 
    theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line=element_line(colour="black"),axis.text = element_text(colour="black"),axis.text.x=element_text(angle = 90,hjust=1), ... )
}
theme_AL_box <- function (...) { 
  theme_bw(base_size=12, base_family="") %+replace% 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(colour="black"),axis.line=element_line(colour="black"), ... )
}
theme_AL_box_rotX <- function (...) { 
  theme_bw(base_size=12, base_family="") %+replace% 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(colour="black"),axis.line=element_line(colour="black"), axis.text.x=element_text(angle = 90,hjust=1), ... )
}
theme_AL_box_rotX_45 <- function (...) { 
  theme_bw(base_size=12, base_family="") %+replace% 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text = element_text(colour="black"),axis.line=element_line(colour="black"), axis.text.x=element_text(angle = 45), ... )
}
theme_AL_simple_rotX_45 <- function (...) { 
  theme_bw(base_size=12, base_family="") %+replace% 
    theme(panel.border = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line=element_line(colour="black"),axis.text = element_text(colour="black"),axis.text.x=element_text(angle = 45), ... )
}

#Extracts estimated RNA-seq count data from a sleuth object and summarizes it to gene-level using the tximport algorithm
#Uses the data.table package for speed but can be changed to base::data.frame eg. using reshape2::dcast
summarizeSleuthToGene <- function(sleuth_object, tx2gene, which_df = c("obs_norm_filt", "obs_norm", "obs_raw"), ... ){
  if( class(sleuth_object) != 'sleuth' ) stop( '"sleuth_object" must be a sleuth object' )
  if( missing(sleuth_object) ) stop( 'argument "tx2gene" is missing with no default' )
  which_df <- match.arg(which_df, c("obs_norm_filt", "obs_norm", "obs_raw") )
  require(sleuth)
  require(tximport)
  require(data.table)
  dat <- data.table(sleuth_object[[which_df]])
  cnt_mat <- as.matrix(data.frame(dcast.data.table(data = dat, formula = target_id ~ sample, value.var = "est_counts"),row.names = 1))
  tpm_mat <- as.matrix(data.frame(dcast.data.table(data = dat, formula = target_id ~ sample, value.var = "tpm"),row.names = 1))
  len_mat <- as.matrix(data.frame(dcast.data.table(data = dat, formula = target_id ~ sample, value.var = "eff_len"),row.names = 1))
  txi <- list(abundance=tpm_mat,counts=cnt_mat,length=len_mat)
  summarizeToGene(txi,tx2gene=tx2gene, ... )
}

# Get distance from PCA
get.pca.dist <- function(txi,meta){
  library(ggfortify)
  pc <- prcomp(t(txi$counts))
  #autoplot(pc,data=meta,colour="Celltype")
  pc_avg <- aggregate(pc$x,list(meta$Pseudotime),mean)
  #plot(pc$x[,1:2],pch=19,col=rainbow(6)[as.numeric(factor(meta$Celltype))]);points(pc_avg[,2:3],pch=as.character(seq_along(pc_avg[,1])),col="red",cex=2)
  dst <- dist(pc_avg)
  idx <- which(!duplicated(combn(pc_avg[,1],2)[1,]))
  time_dist <- c(0,dst[idx-c(0,0,0,0,1)])
  names(time_dist) <- c("ETP","DN","DPearly","DPlate","CD4SP","CD8SP")
  
  time_dist_sum <- cumsum(time_dist)
  time_dist_sum[6] <- cumsum(time_dist[-5])[5]
  
  # Takes ~2-3 weeks from ETP --> SP
  nf <- max(time_dist_sum) / 420 # hours over 2.5 w
  time_h <- round(time_dist_sum/nf)
  time_h <- time_h[order(names(time_h))]
  return(time_h)
}


chromz <- c("chr1", "chr2", "chr3" ,"chr4", "chr5" , "chr6" ,"chr7" , "chr8" , "chr9" ,"chr10" , "chr11" , "chr12" , "chr13" ,"chr14" ,"chr15" ,"chr16" , "chr17", "chr18" , "chr19" , "chr20", "chr21" , "chr22" ,"chrX")

chrom_colors <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#7E6148FF","#B09C85FF","#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF","#7E6148FF","#B09C85FF","#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF", "#F39B7FFF")

chrom_colors_gg <- c("1" = "#E64B35FF","2" = "#4DBBD5FF","3" ="#00A087FF","4" ="#3C5488FF","5" ="#F39B7FFF","6" ="#8491B4FF","7" ="#91D1C2FF","8" ="#7E6148FF","9" ="#B09C85FF","10" ="#E64B35FF","11" ="#4DBBD5FF","12" ="#00A087FF","13" ="#3C5488FF","14" ="#F39B7FFF","15" ="#8491B4FF","16" ="#91D1C2FF","17" ="#7E6148FF","18" ="#B09C85FF","19" ="#E64B35FF","20" ="#4DBBD5FF","21" ="#00A087FF","22" ="#3C5488FF", "X" ="#F39B7FFF")


color_valzz <- c("brain-1" = "#0072B5FF",
                 "brain-2" ="#0072B5FF" ,
                 "brain-3" ="#0072B5FF",
                 "brain-4" = "#0072B5FF", 
                 "brain-5" = "#0072B5FF" ,
                 "brain-6" = "#0072B5FF" ,
                 "brain-7" = "#0072B5FF" ,
                 "brain-8" = "#0072B5FF" ,
                 "brain-9" = "#0072B5FF",
                 "brain-10" = "#0072B5FF",
                 "brain-11" = "#0072B5FF",
                 "brain-12" = "#0072B5FF",
                 "brain-13" ="#0072B5FF",
                 "esophagus-1" = "#4DBBD5FF" ,
                 "esophagus-2" ="#4DBBD5FF",
                 "esophagus-3" = "#4DBBD5FF",
                 "skin-1" ="#00A087FF",
                 "skin-2" = "#00A087FF",
                 "heart-1" = "#3C5488FF",
                 "heart-2" = "#3C5488FF",
                 "adipose-1" = "#F39B7FFF",
                 "adipose-2" ="#F39B7FFF",
                 "artery-1" = "#91D1C2FF",
                 "artery-2" = "#91D1C2FF" ,
                 "artery-3" = "#91D1C2FF",
                 "adrenal gland" = "#8491B4FF",
                 "kidney-1" = "#7E6148FF",
                 "kidney-2" = "#7E6148FF",
                 "bladder" = "#DC0000FF",
                 "nerve" = "#B09C85FF",
                 "salivary gland" = "#BC3C29FF",
                 "throid" = "#E64B35FF",
                 "pituitary" = "#E18727FF",
                 "pancreas" = "#20854EFF",
                 "spleen" = "#7876B1FF",
                 "liver" = "#6F99ADFF",
                 "lung" = "#FFDC91FF",
                 "breast" ="#EE4C97FF",
                 "muscle" = "#374E55FF",
                 "small intestine" = "#DF8F44FF",
                 "colon-1" = "#DF8F44FF",
                 "colon-2" = "#DF8F44FF",
                 "stomach" ="#DF8F44FF",
                 "lymphocytes" = "grey",
                 "fibroblasts" = "#B24745FF",
                 "whole blood" = "#79AF97FF",
                 "thymocytes" = "red")


# patient_order
patient_order <- c("M1", "M2", "M3", "F1", "F2", "F3")

# cell order
cell_order <- c("ETP","DN", "DPearly", "DPlate", "CD4SP", "CD8SP")
cell_order_alt <- c("ETP", "DN", "DPE", "DPL", "SP4", "SP8")


# cell colors
colors = c("ETP" = "#E31C1E", "DN" = "#EF7D19", "DPearly" = "#4EAF4A", "DPlate" = "#367EB9", "CD4SP" = "#EA82B3", "CD8SP" = "#944F99")

# patient colors
patient_colors <- c("F1" = "#000000", "F2" = "#E69F00", "F3" = "#56B4E9", "M1" = "#009E73",
                    "M2" =  "#0072B2", "M3" = "#CC79A7", "F3 (v6 + UTR)" = "#D55E00")

# Chromosome order
chrom_order <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22", "chrX")