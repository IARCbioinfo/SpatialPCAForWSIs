#!/bin/bash
#SBATCH --job-name=ClosestCentroids
#SBATCH --output=ClosestCentroids-%j.out
#SBATCH --error=ClosestCentroids-%j.error
#SBATCH --mem=20G
#SBATCH --partition=high_p
#SBATCH --account=gem

Rscript ClosestCentroids.R  --sample_id TNE0001 \
                            --proj_tab_SPCA ~/LNENWork/LNEN_Molecular_groups_barlow_twin/spatial_pca/ProjBySampleID_v4/ProjSpatialPCA_allSamples_50k_v4.csv \
                            --outdir LeidenCom_Spc50k_n100000_knn6000_res01_CC \
                            --centroids_tab LeidenCom_Spc50k_n100000_knn6000_res01/SPCA_centroids_leiden_ntiles100000_KNN_6000_Res_01_ncluster_164.csv

