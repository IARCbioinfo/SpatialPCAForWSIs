#!/bin/bash
#SBATCH --job-name=ProjSPCA50K
#SBATCH --output=ProjSPCA50K-%j.out
#SBATCH --error=ProjSPCA50K-%j.error
#SBATCH --mem=180G
#SBATCH --partition=high_p
#SBATCH --account=gcs

Rscript ProjectionElaboratedByPatient.R --sample_id TNE0001  \
                                    --spca_obj ~/LNENWork/LNEN_Molecular_groups_barlow_twin/spatial_pca/SpatialPCAObject_50k_V4.RData \
                                    --dimension_enc_vectors 128 \
                                    --proj_tab_norm /home/mathiane/LNENWork/LNEN_Molecular_groups_barlow_twin/projector/barlow_twins_encoded_vectors_for_spatial_PCA_scale_path_review.csv \
                                    --outdir spatial_PCA_proj_by_sampleLeidenCommu


#