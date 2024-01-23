print("#-------------------------------------- Install package ------------------------------------------------------------------------")
print(.libPaths())
.libPaths("/home/mathiane/R/x86_64-pc-linux-gnu-library/4.0")
packages <- c("stringr",  "cluster", "Matrix","SpatialPCA", "optparse" , "tictoc" )
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
print("# --------------------------------------- Load libraies ------------------------------------------------------------------------")
library(rlang)
library(stringr)
library(cluster)
library(Matrix)
library(SpatialPCA)
library("optparse")
library(tictoc)

source('ImgSpatialPCAMultipleSamles.R')
environment(SpatialPCA_Multiple_Sample) <- asNamespace('SpatialPCA')
assignInNamespace("SpatialPCA_Multiple_Sample", SpatialPCA_Multiple_Sample, ns = "SpatialPCA")


source('ImgSpatialPCA.R')
environment(CreateSpatialPCAObject) <- asNamespace('SpatialPCA')
assignInNamespace("CreateSpatialPCAObject", CreateSpatialPCAObject, ns = "SpatialPCA")

#setwd("~/script_spatial_pca")
tic("total")

print("# --------------------------------- Command line arguments ------------------------------------------------------------------------")
option_list = list(
  make_option(c("-n", "--n_tiles"), type="numeric", 
              help="Number of tiles selected for the PCA", metavar="number"),
  make_option(c("-o", "--output_folder"), type="character",
              help="Folder where to store the results", metavar="character"),
  make_option(c("-p", "--path2projectors"), type="character",
              help="Path to the csv file containing all the encoded vectors", metavar="character"),
  make_option(c("-d", "--dimension_enc_vectors"), type="integer", default=128, 
        help="Dimension of encoded vectors")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
n_tiles = opt$n_tiles #10000
output_folder = opt$output_folder
path2projectors = opt$path2projectors
d_enc_vectors = opt$dimension_enc_vectors
# -------------------- My function -----------------------
correct_tiles_id <- function(df_enc_vector_c){
  # Tiles are expected to have an id such as patientID_coordsX_coordsY.jpg
  # This function garantee that tile ID respect the pattern mentionned above.
  tiles_id <- unlist(lapply(df_enc_vector_c$img_id, function(x)   paste(str_split(x, pattern ='_')[[1]][1], 
                                                                        str_split(x, pattern ='_')[[1]][2], 
                                                                        str_split(x, pattern ='_')[[1]][3], sep= '_') ))
  return(tiles_id)
}

print("# --------------------------------- Read scales encoded vectors ------------------------------------------------------------------------")
df_enc_vector <- read.csv(path2projectors)
colnames(my_dataframe)[colnames(my_dataframe) == "tne_id_c"] ="sample_id"
df_enc_vector <- df_enc_vector[!duplicated(df_enc_vector$img_id),]
df_enc_vector <- df_enc_vector[,2:ncol(df_enc_vector)]

print("#--------------------- Samples n number of tiles to build the PCA ---------------------------------------------------------")
select_tiles_v <-  c(rep(0,nrow(df_enc_vector)-n_tiles), rep(1,n_tiles))
select_tiles_v <- sample(select_tiles_v)
df_enc_vector$TilesForPCA <- select_tiles_v

df_enc_vectors_to_pca <- df_enc_vector[df_enc_vector$TilesForPCA == 1,]
df_enc_vectors_to_proj <- df_enc_vector[df_enc_vector$TilesForPCA == 0,]
rm(df_enc_vector)

print("# -----------------------List of tables  --------------------------------------------")
# Get all sample name
Samples_from_PCA <- df_enc_vectors_to_pca$sample_id 
print(unique(Samples_from_PCA))
print(length(unique(Samples_from_PCA)))

# Remove Samples with less than 3 tiles  from spatial PCA
for (sampl in unique(Samples_from_PCA)){
  if(nrow(df_enc_vectors_to_pca[df_enc_vectors_to_pca$sample_id  == sampl,]) < 3){
    print(sampl)
    c_sampl <- df_enc_vectors_to_pca[df_enc_vectors_to_pca$sample_id  == sampl, ]
    df_enc_vectors_to_proj <- rbind(df_enc_vectors_to_proj, c_sampl)
    df_enc_vectors_to_pca <- df_enc_vectors_to_pca[df_enc_vectors_to_pca$sample_id  != sampl, ]
  }
}


# Create lists of data frames for coordinates and encoded vectors 
Samples_from_PCA <- df_enc_vectors_to_pca$sample_id
coords_df_list_for_pca <- list()
encoded_vectors_df_list_for_pca <- list()
c = 1
for (ele in unique(Samples_from_PCA)){
  # Coords data frame
  df_enc_vector_c_for_PCA <- df_enc_vectors_to_pca[df_enc_vectors_to_pca$sample_id == ele,]
  # Corrds matrices  
  coords_df_c_for_PCA  <- data.frame(x_coord=as.numeric(df_enc_vector_c_for_PCA$x))
  coords_df_c_for_PCA$y_coord <- as.numeric(df_enc_vector_c_for_PCA$y)
  rownames(coords_df_c_for_PCA) = df_enc_vector_c_for_PCA$img_id_c
  coords_df_c_for_PCA = as.matrix(coords_df_c_for_PCA)
  coords_df_list_for_pca[[c]] <- coords_df_c_for_PCA
  
  # Encoded vecotrs data frame
  enc_val_for_pca = df_enc_vector_c_for_PCA[,1:d_enc_vectors]
  enc_val_t_for_pca <- t(enc_val_for_pca)
  enc_val_t_for_pca <- Matrix(nrow = nrow(enc_val_t_for_pca), ncol = ncol(enc_val_t_for_pca), data = enc_val_t_for_pca, sparse = TRUE)
  colnames(enc_val_t_for_pca) <- df_enc_vector_c_for_PCA$img_id_c
  rownames(enc_val_t_for_pca) <- colnames(enc_val_for_pca)[1:d_enc_vectors]
  encoded_vectors_df_list_for_pca[[c]] <- enc_val_t_for_pca
  
  c = c + 1
}

print("# -------------------------------------- Compute spatial PCA -------------------------------------")

LIBD_multi = SpatialPCA_Multiple_Sample(encoded_vectors_df_list_for_pca, coords_df_list_for_pca ,gene.type="spatial",sparkversion="spark",
                                        numCores_spark=5, gene.number=d_enc_vectors, customGenelist=NULL, min.loctions = 20, min.features=20, bandwidth_common=0.1)

# Save spatial PCA in a R obj
save(LIBD_multi, file =  paste(output_folder, "SpatialPCAObject_50k_V4.RData", sep='/'))

multi_spatial_pcs <- LIBD_multi$MultipleSample_SpatialPCA_object@SpatialPCs
multi_spatial_pcs <- t(multi_spatial_pcs)
colnames(multi_spatial_pcs)<-unlist(lapply(1:20, function(x) paste('axis', x, sep='_')))
multi_spatial_pcs <- data.frame(multi_spatial_pcs)
img_id <- unlist(lapply( rownames(multi_spatial_pcs), function(x)  substr( x, str_locate(x, 'TNE')[[1]], nchar(x)))) 
rownames(multi_spatial_pcs) <- img_id

samples_id <- unlist(lapply(rownames(multi_spatial_pcs), function(x) str_split(x, '_')[[1]][1] ))
unique(samples_id)
multi_spatial_pcs$sample_id <- samples_id

multi_spatial_pcs <- merge(multi_spatial_pcs, df_enc_vectors_to_pca[,129:ncol(df_enc_vectors_to_pca)], by.x="row.names", by.y="img_id_c")

write.csv(multi_spatial_pcs, paste(output_folder, "spatial_pca_coords_50k_V4.csv", sep='/'))  

#output_pca_filename <- paste(paste("MultiSpatialPCA_3D", TNE_name, sep  =''), '.csv', sep='')
#write_csv(multi_spatial_pcs, output_pca_filename)

#  ------------------------ See projection by sample ID ----------------------------

# ggplot(multi_spatial_pcs, aes(axis_1, axis_2)) +
#   geom_point(aes(colour = samples_id), size = 0.5, alpha = 0.8)+ theme(legend.position = "none")
# 
# ggplot(multi_spatial_pcs, aes(axis_1, axis_3)) +
#   geom_point(aes(colour = samples_id), size = 0.5, alpha = 0.8)+ theme(legend.position = "none")
# 
# 
# ggplot(multi_spatial_pcs, aes(axis_2, axis_3)) +
#   geom_point(aes(colour = samples_id), size = 0.5, alpha = 0.8)+ theme(legend.position = "none")

# ------------------------ Computation louvain clusters  -----------------------
# clusterlabel_all = louvain_clustering(clusternum=10,latent_dat=LIBD_multi$MultipleSample_SpatialPCA_object@SpatialPCs,knearest=70 ) 
# multi_spatial_pcs$cluster <- clusterlabel_all
# 
# ggplot(multi_spatial_pcs, aes(axis_1, axis_2)) +
#   geom_point(aes(colour = cluster), size = 0.5, alpha = 0.8) 
# 
# ggplot(multi_spatial_pcs, aes(axis_1, axis_3)) +
#   geom_point(aes(colour = cluster), size = 0.5, alpha = 0.8) 
# 
# ggplot(multi_spatial_pcs, aes(axis_2, axis_3)) +
#   geom_point(aes(colour = cluster), size = 0.5, alpha = 0.8) 


# ------------------------ Individual representation of the PCA ------------------------------------------------------
# cbp=c("#9C9EDE" ,"#5CB85C" ,"#E377C2", "#4DBBD5" ,"#FED439" ,"#FF9896", "#FFDC91","#00AFBB", "#E7B800", "#FC4E07")
# sample_int_id = 247
# unique(multi_spatial_pcs$samples_id)
# multi_spatial_pcs$x <- unlist(lapply(rownames(multi_spatial_pcs), function(x) as.numeric(str_split(x, '_')[[sample_int_id]][2])))
# multi_spatial_pcs$y <- unlist(lapply(rownames(multi_spatial_pcs), function(x) as.numeric(str_split(x, '_')[[sample_int_id]][3])))
# 
# coords_table_sample  <- multi_spatial_pcs[multi_spatial_pcs$samples_id == unique(multi_spatial_pcs$samples_id)[sample_int_id], names(multi_spatial_pcs) %in% c('x', 'y')]
# dim(coords_table_sample)
# coords_label_sample <-  multi_spatial_pcs[multi_spatial_pcs$samples_id == unique(multi_spatial_pcs$samples_id)[sample_int_id],]$cluster
# plot_cluster(location=coords_table_sample,clusterlabel=coords_label_sample,pointsize=1.5,title_in=paste0(unique(multi_spatial_pcs$samples_id)[sample_int_id]),color_in=cbp)
# 

rm(multi_spatial_pcs)
rm(df_enc_vectors_to_pca)
#print("# -------------------------------- Projection --------------------------------------------------------")

## TO REMOVE
# df_enc_vectors_to_proj <- df_enc_vectors_to_proj[sample(1:nrow(df_enc_vectors_to_proj), 1000000), ]
##
# enc_vectors_to_proj <- t(df_enc_vectors_to_proj[,1:128])
# dim(enc_vectors_to_proj)
# enc_vectors_to_proj_info <- df_enc_vectors_to_proj[,129:ncol(df_enc_vectors_to_proj)]
# head(enc_vectors_to_proj_info)
# 
# Z_PCA <- t(LIBD_multi$MultipleSample_SpatialPCA_object@W) %*% enc_vectors_to_proj
# Z_PCA_c <- Z_PCA
# Z_PCA <- t(Z_PCA)
# colnames(Z_PCA) <-unlist(lapply(1:3, function(x) paste('axis', x, sep='_')))
# 
# Z_PCA <- cbind(Z_PCA,enc_vectors_to_proj_info)
# 
# write.csv(Z_PCA, paste(output_folder, "spatial_pca_coords_projection.csv", sep='/'))  


# ------------------------------------ Compute Louvain clusters on a random partition --------------
## Warning: Add tiles include to build the spatial PCA
# tiles_for_louvain <-sample(1:1000000, 40000, replace=TRUE)
# Z_PCA_sampling_louvain_com  <- Z_PCA_c[,tiles_for_louvain]
# clusterlabel_samples_proj = louvain_clustering(clusternum=10,latent_dat=Z_PCA_sampling_louvain_com,knearest=70 ) 
# Z_PCA_sampling_louvain_com_clst <- data.frame(t(Z_PCA_sampling_louvain_com))
# Z_PCA_sampling_louvain_com_clst$cluster <- clusterlabel_samples_proj
# 
# head(Z_PCA_sampling_louvain_com_clst)
# enc_vectors_to_proj_info_for_louvain <- enc_vectors_to_proj_info[tiles_for_louvain,]
# 
# Z_PCA_sampling_louvain_com_clst <- cbind(Z_PCA_sampling_louvain_com_clst, enc_vectors_to_proj_info_for_louvain)
# 
# ggplot(Z_PCA_sampling_louvain_com_clst, aes(X1,X2)) +
#   geom_point(aes(colour = as.factor(cluster)), size = 1, alpha = 0.8) +  theme(legend.position = "none")
# 
# ggplot(Z_PCA_sampling_louvain_com_clst, aes(X1,X3)) +
#   geom_point(aes(colour = as.factor(cluster)), size = 1, alpha = 0.8) +  theme(legend.position = "none")
# 
# ggplot(Z_PCA_sampling_louvain_com_clst, aes(X2,X3)) +
#   geom_point(aes(colour = as.factor(cluster)), size = 1, alpha = 0.8) +  theme(legend.position = "none")
# 
# 
# ## Conpute centroids by cluster 
# list_cluster <- c()
# list_centroids_axis1 <- c()
# list_centroids_axis2 <- c()
# list_centroids_axis3 <- c()
# for (clst in unique(Z_PCA_sampling_louvain_com_clst$cluster)){
#   Z_PCA_sampling_louvain_com_clst_c <- Z_PCA_sampling_louvain_com_clst[Z_PCA_sampling_louvain_com_clst$cluster == clst,]
#   axis_1_cent <- sum(Z_PCA_sampling_louvain_com_clst_c$X1)/dim(Z_PCA_sampling_louvain_com_clst_c)[1]
#   axis_2_cent <- sum(Z_PCA_sampling_louvain_com_clst_c$X2)/dim(Z_PCA_sampling_louvain_com_clst_c)[1]
#   axis_3_cent <- sum(Z_PCA_sampling_louvain_com_clst_c$X3)/dim(Z_PCA_sampling_louvain_com_clst_c)[1]
#   list_cluster <- c(list_cluster, clst)
#   list_centroids_axis1 <- c(list_centroids_axis1 , axis_1_cent)
#   list_centroids_axis2 <- c(list_centroids_axis2 ,axis_2_cent)
#   list_centroids_axis3 <- c(list_centroids_axis3, axis_3_cent)
# }
# 
# df_centroids <- data.frame(cluster= list_cluster,
#                            axis1_centrois = list_centroids_axis1,
#                            axis2_centrois = list_centroids_axis2,
#                            axis3_centrois = list_centroids_axis3)
# 
# Z_PCA_sampling_louvain_centroids  <- t(Z_PCA_c[,-c(tiles_for_louvain)])
# enc_vectors_to_proj_info_cent <- enc_vectors_to_proj_info[-c(tiles_for_louvain),]
# Z_PCA_sampling_louvain_centroids <- cbind(Z_PCA_sampling_louvain_centroids, enc_vectors_to_proj_info_cent)
# head(Z_PCA_sampling_louvain_centroids)
# colnames(Z_PCA_sampling_louvain_centroids)[1] <- "axis1"
# colnames(Z_PCA_sampling_louvain_centroids)[2] <- "axis2"
# colnames(Z_PCA_sampling_louvain_centroids)[3] <- "axis3"
# ## To remove
# 
# Z_PCA_sampling_louvain_centroids <- Z_PCA_sampling_louvain_centroids[sample(1:nrow(Z_PCA_sampling_louvain_centroids), 500000),]
# ## Dist to cluster
# 
# 
# # clst 2
# closest_centroid<-  function(x,df_centroids, Z_PCA_sampling_louvain_centroids){
#   for (row_centroid in 1:nrow(df_centroids)){
#   Z_PCA_sampling_louvain_centroids[x,]
#   df_centroids_c_cluster =  df_centroids[row_centroid,]
#   d_to_clust <- sqrt((Z_PCA_sampling_louvain_centroids[x,1]-df_centroids_c_cluster$axis1_centrois)**2 +
#                        (Z_PCA_sampling_louvain_centroids[x,2]-df_centroids_c_cluster$axis2_centrois)**2 +
#                        (Z_PCA_sampling_louvain_centroids[x,3]-df_centroids_c_cluster$axis3_centrois)**2)
#   dist <- c(dist, d_to_clust)
#   return(df_centroids[which.min(dist),]$cluster)
#   }
# }
# 
# centroids_pred <- c()
# for(i in 1:nrow(Z_PCA_sampling_louvain_centroids)){
#   c_clst <- closest_centroid(i,df_centroids, Z_PCA_sampling_louvain_centroids)
#   centroids_pred <- c(centroids_pred, c_clst)
# }
toc()
