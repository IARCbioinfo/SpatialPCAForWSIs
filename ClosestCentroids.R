print(.libPaths())
.libPaths("/home/mathiane/R/x86_64-pc-linux-gnu-library/4.0")

print("# --------------------------------------- Load libraies ------------------------------------------------------------------------")
library(stringr)
library("optparse")
library(tictoc)
library(ggplot2)
library(gridExtra)
library(readxl)
library(dplyr)



print("# --------------------------------- Command line arguments ------------------------------------------")

option_list = list(
  make_option(c("--sample_id"), type="character",
              help="Patient ID matching WSI filename", metavar="character"),
  make_option(c("--proj_tab_SPCA"), type="character",
              help="Path to the projectors table produced by the spatial PCA.", metavar="character"),
  make_option(c("--outdir"), type="character", 
              help="Folder where the results will be store", metavar="character"),
  make_option(c("--centroids_tab"), type="character",
              help="Path to the table of centroids of Leiden communities", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

sample_id = opt$sample_id 
exp_folder = opt$outdir 
dir.create(exp_folder, showWarnings = FALSE)

print("# --------------------------------- Load spatial PCA projections ------------------------------------------")
proj_path_review <- read.csv(opt$proj_tab_SPCA)


head(proj_path_review)
set.seed(123458)
df_enc_vector <- proj_path_review 
# Extract projections of the petient of nterest
df_enc_vector = df_enc_vector[df_enc_vector$sample_id == sample_id,]
head(df_enc_vector)

print("# ---------- Load Leiden community centroids ----------------------------------------------------------------------------")
df_centroids = read.csv(opt$centroids_tab)
head(df_centroids)

df_centroids <- df_centroids[order(df_centroids$cluster),] 
df_centroids$cluster
length(unique(df_centroids$cluster))


df_enc_vector <- df_enc_vector[,2:ncol(df_enc_vector)]

n_clusters_louvain = dim(df_centroids)[1]


closest_centroid<-  function(x){
  # Find the closest centroids to each spatial PCA projection
  dist <- c()
  for (row_centroid in 1:nrow(df_centroids)){
    df_centroids_c_cluster =  df_centroids[row_centroid,]
    sum_of_dist <-c()
    for(c in 2:21){
      sum_of_dist <-c(sum_of_dist, (df_enc_vector[x,c]-df_centroids_c_cluster[,c])**2)
    }
    d_to_clust <- sqrt(sum(sum_of_dist))
    dist <- c(dist, d_to_clust)
  }
  return(df_centroids[which.min(dist),]$cluster)
}


centroids_pred <-   lapply(seq(1,nrow(df_enc_vector)), function(x)  closest_centroid(x)) 
# add clusters assignation to the spatial PCA table
df_enc_vector$cluster = unlist(centroids_pred)
df_enc_vector$cluster = as.factor(df_enc_vector$cluster)


print("-------------------------------------- Save results ------------------------------------------------")
table_name = paste(exp_folder, paste( paste0( "SPCA_centroids_leiden_ntiles100000_KNN_6000_Res_01_", sample_id), paste( n_clusters_louvain, '.csv', sep = '') , sep= '_') , sep='/')
write.csv(df_enc_vector,table_name,  row.names = F)

