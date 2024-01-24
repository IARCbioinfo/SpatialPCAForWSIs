#!/bin/bash
#SBATCH --job-name=LeidenSPCA_50K_res01
#SBATCH --output=LeidenSPCA_50K_res01-%j.out
#SBATCH --error=LeidenSPCA_50K_res001-%j.error
#SBATCH --mem=300G
#SBATCH --partition=high_p
#SBATCH --account=gem

Rscript LeidenCommunitySpatialPCA.R  --Resolution 0.1 \
                                     --ntiles 100000 \
                                     --KNN 6000 \
                                     --proj_tab_SPCA ~/LNENWork/LNEN_Molecular_groups_barlow_twin/spatial_pca/ProjBySampleID_v4/ProjSpatialPCA_allSamples_50k_v4.csv \
                                     --outdir LeidenCom_Spc50k_n100000_knn6000_res01

