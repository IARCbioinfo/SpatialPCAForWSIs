###################################################################################################################################
# Authors: Emilie Mathian 
# This script searches for the Leiden community among the low-dimensional projections (20 dimensions) created by spatial PCA. 
# The number of clusters can be modulated by varying the `KNN` parameter, i.e. the number of nodes linked together as a function of their distance, when creating the graph; 
# and `Resolustion` when detecting the Leiden community. A smaller value for `Resolution` will tend to increase the number of communities and vice versa. 
# Time and memory costs increase with the number of points (vectors) considered to create the graph, see the ntiles parameter. 
# The centroids of the communities are saved in a csv file, so that the communities can then be assigned to other representations.
# based on the minimum distance of these vectors from the community centroids.

# Contacts: mathiane@iarc.who.int
#          Rare cancer genomics teams, International Agency of Research on Cancer - WHO 
####################################################################################################

print("#-------------------------------------- Install package ------------------------------------------------------------------------")
print(.libPaths())
.libPaths("/home/mathiane/R/x86_64-pc-linux-gnu-library/4.0")
packages <- c("stringr",  "cluster", "Matrix","SpatialPCA", "optparse" , "tictoc" , "ggplot2", "gridExtra",  "readxl", "dplyr", "igraph", "parallel")
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
print("# --------------------------------------- Load libraies ------------------------------------------------------------------------")
library(stringr)
library(cluster)
library(Matrix)
library("optparse")
library(tictoc)
library(SpatialPCA)
library(ggplot2)
library(gridExtra)
library(readxl)
library(dplyr)
library(igraph)
library(parallel)


print("# --------------------------------- Command line arguments ------------------------------------------")

option_list = list(
  make_option(c("--ntiles"), type="integer", default = 100000,
              help="Number of tiles"),
  make_option(c( "--KNN"), type="integer", default = 6000,
              help="Number of neighbors for the KNN graph"),
  make_option(c("--Resolution"), type="numeric", default = 0.1,
              help="Resolution for leiden"),
  make_option(c("--proj_tab_SPCA"), type="character",
              help="Path to the projectors table produced by the spatial PCA.", metavar="character"),
  make_option(c("--outdir"), type="character",
              help="Folder where the results will be store", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

n_tiles_louvain =opt$ntiles 
K =opt$KNN 
Resolution =opt$Resolution 
exp_folder = opt$outdir  
dir.create(exp_folder, showWarnings = FALSE)
print("# -------------------- Load data -------------------------")
leiden_clustering = function( latent_dat,knearest=100, Resolution=0.02){
  set.seed(1234)
  
  PCvalues = latent_dat
  dim(PCvalues)
  info.spatial = as.data.frame(t(PCvalues))
  colnames(info.spatial) =  paste0("factor", 1:nrow(PCvalues))
  knn.norm = FNN::get.knn(as.matrix(t(PCvalues)), k = knearest)
  print("KNN")
  knn.norm = data.frame(from = rep(1:nrow(knn.norm$nn.index),
                                   k=knearest), to = as.vector(knn.norm$nn.index), weight = 1/(1 + as.vector(knn.norm$nn.dist)))
  nw.norm = igraph::graph_from_data_frame(knn.norm, directed = FALSE)
  is_weighted(nw.norm)
  nw.norm = igraph::simplify(nw.norm)
  print("Computing Leiden clusters ")
  lc.norm = igraph::cluster_leiden(nw.norm, resolution = Resolution)
  print("Nb louvain community")
  clusterlabel <-  as.vector(igraph::membership(lc.norm))
  
  print("Nb community after merging")
  print(unique(clusterlabel))
  print(length(unique(clusterlabel)))
  
  return("cluster_label"=clusterlabel)
}
# 
correct_tiles_id <- function(df_enc_vector_c){
  tiles_id <- unlist(lapply(df_enc_vector_c$X, function(x)   paste(substr(x,1,7), 
                                                                        str_split(x, pattern ='_')[[1]][2], 
                                                                        str_split(x, pattern ='_')[[1]][3], sep= '_') ))
  return(tiles_id)
}


print("------------------------------------ # load the projections created by the spatial PCA -----------------------------")
proj_path_review <- read.csv(opt$proj_tab_SPCA)

head(proj_path_review)
set.seed(123458)
df_enc_vector <- proj_path_review 
print("------------------------------------ # Sample n tiles to compute Leiden Communities -----------------------------")

tiles_for_louvain <-sample(1:nrow(df_enc_vector), n_tiles_louvain , replace=TRUE)#
df_enc_vector_to_louvain <- df_enc_vector[tiles_for_louvain, ]
# Get vectors associated with the spatial PCA
# Warning here we are supposing that "axis_1" matches with the 3rd column and that "axis_20" matches with the 22nd column.
df_enc_vector_val <- df_enc_vector_to_louvain[,3:22]
# Extract tiles info in another table
df_enc_vector_info <- data.frame("img_id" = df_enc_vector_to_louvain[,2])


df_enc_vector_val <- t(df_enc_vector_val)
colnames(df_enc_vector_val) <- df_enc_vector_info$img_id


print("------------------------------------ # Compute Leiden communities -----------------------------")

clusterlabel_samples_proj = leiden_clustering(latent_dat=df_enc_vector_val, knearest=K, Resolution=Resolution) 

# Add clusters to the tables
df_enc_vector_to_louvain$cluster <- clusterlabel_samples_proj
n_clusters_leiden = length(unique(df_enc_vector_to_louvain$cluster))
df_enc_vector_to_louvain <- df_enc_vector_to_louvain[,2:ncol(df_enc_vector_to_louvain)]

print("# ------------------------------------ Compute Louvain Centroids -----------------------------")
centroid_by_axis <- function(axis, n_tot){
  return (sum(axis)/n_tot)
}


c = 0
for (clst in unique(df_enc_vector_to_louvain$cluster)){
  print(clst)
  df_enc_vector_to_louvain_c <- df_enc_vector_to_louvain[df_enc_vector_to_louvain$cluster == clst,]
  pos_by_axis = c()
  
  data_frame_centroids_c = data.frame(cluster=clst)
  for(i in 2:21){
    pos_by_axis_c <-  centroid_by_axis(df_enc_vector_to_louvain_c[,i], dim(df_enc_vector_to_louvain_c)[1] )
    data_frame_centroids_c = cbind(data_frame_centroids_c, pos_by_axis_c)
    
  }
  colnames(data_frame_centroids_c)[2:21] = unlist(lapply(2:21, function(x) paste('axis_', x-1, sep=''))) 
  if (c==0){
    df_centroids = data_frame_centroids_c
  }
  else{
    df_centroids = rbind(df_centroids, data_frame_centroids_c)
  }
  c= c+1
}

c
centroids_file_name = paste(exp_folder , paste0(
                                  paste0(
                                  paste0(
                                    paste0(
                                      paste0( 
                                        paste0(
                                          paste0(
                                            paste0('SPCA_centroids_leiden_ntiles'
                                                   ,n_tiles_louvain)
                                                    ,"_KNN_"),
                                                    K),
                                                    "_Res_"),
                                                    gsub("\\.", "", as.character(Resolution))), 
                                                    "_ncluster_"),
                                                     n_clusters_leiden),
                                                   '.csv'), sep = '/')
df_centroids
dim(df_centroids)

write.csv(df_centroids,centroids_file_name, row.names = F )
