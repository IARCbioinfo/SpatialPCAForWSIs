########################################################################################################################
# Package: SpatialPCA
# Version: 1.2.0
# Date   : 2022-8-12
# Title : Spatially Aware Dimension Reduction for Spatial Transcriptomics
# Authors: L. Shang and X. Zhou
# Modifications and adaptation for WSI : Emilie Mathian
# Contacts: shanglu@umich.edu
#          University of Michigan, Department of Biostatistics
#          For WSIs: mathiane@iarc.who.int
#          Rare cancer genomics teams, International Agency of Research on Cancer - WHO 
####################################################################################################


#' Multiple sample SpatialPCA.
#' In this extension, we construct the covariance matrix for the latent factors in the form of a block diagonal matrix:
#' it consists of the kernel matrices constructed using the spatial location information within each dataset, with zero correlation for pairs of locations across datasets.
#' This way, the latent factors within each dataset are correlated a priori across spatial locations, while the latent factors across datasets are not correlated a priori.
#' Certainly, if one wants to model the a priori correlation between latent factors across datasets, due to, for example, their similarity in the features extracted from histology images, then one can also modify the kernel matrices by constructing them using features other than spatial location information.
#' @param count_list A list of g by n count matrix, g is gene number, n is location number.
#' @param location_list A list of n by d location matrix, n is location number, d is location dimension. The rownames of each location matrix should match with the colnames of its corresponding count matrix.
#' @return Returns SpatialPCA object with estimated Spatial PCs on locations.
#' @export

SpatialPCA_Multiple_Sample = function(count_list,location_list,gene.type="spatial",sparkversion="spark",numCores_spark=5,gene.number=3000, customGenelist=NULL,min.loctions = 20, min.features=20,bandwidth_common=0.1){
  require(Seurat)
  if(length(count_list) != length(location_list)) {
    stop("The number of count matrix should match with the number of location matrix. ")
  }# end fi
  
  # create seurat object for each dataset
  seurat_list = list()
  for(count in 1:length(count_list)){
    print("************************************** HERE 0 ***********************************************************************")
    print("************************************** Colnames *********************************************************************")
    print(colnames(count_list[[count]]))

    colnames(count_list[[count]]) = paste0("Sample",count,"_",colnames(count_list[[count]]))
    rownames(location_list[[count]]) = colnames(count_list[[count]])
    seurat_list[[count]] <- count_list[[count]] 
  }

  features <-  rownames(count_list[[1]])

  count_df <- data.frame(count_list[1])
  for (i in 2:length(count_list)){
    print("************************************** HERE 1 ***********************************************************************")

    count_df <- cbind(count_df, data.frame(count_list[i]))
  }
  data_m <- data.matrix(count_df)
  integrated_data <- Matrix(nrow = nrow(data_m), ncol =ncol(data_m) , data = data_m)

  # create each SpatialPCA object, this step is essentially for spatial gene selection
  spatialpca_list = list()
  for(count in 1:length(seurat_list)){
    print(paste0("Creating SpatialPCA object for dataset ",count))
    spatialpca_list[[count]] <- CreateSpatialPCAObject(counts=count_list[[count]], location=location_list[[count]], project = "SpatialPCA",gene.type=gene.type,sparkversion=sparkversion,
                                                       numCores_spark=numCores_spark,gene.number=gene.number, customGenelist=customGenelist,min.loctions = min.loctions, min.features=min.features)
  }
  matched_data = list()
  matched_location = list()
  Kernal_mat = list()
  for(count in 1:length(count_list)){
    print(paste0("Creating SpatialPCA kernel for dataset ",count))
    match_gene = na.omit(match(rownames( spatialpca_list[[count]]@normalized_expr), rownames(integrated_data)))
    match_spot = grep(paste0("1_",count),colnames(integrated_data))
    matched_data[[count]] = integrated_data#[match_gene,match_spot]
    matched_location[[count]] = spatialpca_list[[count]]@location

    spatialpca_list[[count]]@normalized_expr = t(scale(t(data.matrix(matched_data[[count]]))))
    spatialpca_list[[count]]@location = scale(matched_location[[count]])
    spatialpca_list[[count]] = SpatialPCA_buildKernel(spatialpca_list[[count]], kerneltype="gaussian", bandwidth.set.by.user=bandwidth_common)
    Kernal_mat[[count]] = spatialpca_list[[count]]@kernelmat

  }
  common_genes = unique(unlist(lapply(matched_data, function(x){ row.names(x)})))
  MultipleSample_merge = spatialpca_list[[1]] # initialize using the first data, will replace with integrated data later
  MultipleSample_merge@normalized_expr =data.matrix(integrated_data)#  t(scale(t(data.matrix(integrated_data))))
  print(str(Kernal_mat))
  options(spam.force64 = TRUE)
  MultipleSample_merge@kernelmat = Matrix::bdiag(Kernal_mat) 
  MultipleSample_merge@kernelmat = as.matrix(MultipleSample_merge@kernelmat)#as.matrix(MultipleSample_merge@kernelmat)

  MultipleSample_merge@params$expr = MultipleSample_merge@normalized_expr
  MultipleSample_merge@location = do.call("rbind",lapply(spatialpca_list, function(x){ x@location}))
  # perform SpatialPCA on integrated expression and integrated block diagonal kernel
  MultipleSample_merge = SpatialPCA_EstimateLoading(MultipleSample_merge,fast=TRUE,SpatialPCnum=20)
  MultipleSample_merge = SpatialPCA_SpatialPCs(MultipleSample_merge, fast=TRUE)
  colnames(MultipleSample_merge@SpatialPCs) = rownames(MultipleSample_merge@location)
  id_list = list()
  SpatialPC_list = list()
  Location_spatialpc_list = list()
  for(sample in 1:length(count_list)){

    id_list[[sample]] = grep(paste0("Sample",sample), colnames(MultipleSample_merge@SpatialPCs))

    SpatialPC_list[[sample]] = MultipleSample_merge@SpatialPCs[,id_list[[sample]]]

    Location_spatialpc_list[[sample]] = spatialpca_list[[sample]]@location

  }
  
  names(SpatialPC_list) = paste0("Sample",1:length(SpatialPC_list))
  names(Location_spatialpc_list) = paste0("Sample",1:length(Location_spatialpc_list))
  return(list("MultipleSample_SpatialPCA_object"=MultipleSample_merge,
              "SpatialPC_list" = SpatialPC_list,
              "Location_spatialpc_list"=Location_spatialpc_list))
}
















