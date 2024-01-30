#!/bin/bash
#SBATCH --job-name=SpatialPCA50K
#SBATCH --output=SpatialPCA50k-%j.out
#SBATCH --error=SpatialPCA50k-%j.error
#SBATCH --mem=300G
#SBATCH --partition=high_p
#SBATCH --account=gcs

Rscript RunMultiSPCARandomSampling.R --n_tiles 50000 \
                                    --output_folder spatial_pca_test \
                                    --path2projectors /home/mathiane/LNENWork/LNEN_Molecular_groups_barlow_twin/projector/projector_concat_path_review.csv \
                                    --dimension_enc_vectors 128 