![PresentationImg](SPCA_ImgPresentation.png)
# Spatial PCA for WSI: 
Spatial Principal Component Analysis (PCA), proposed by [L. Shang and X. Zhou, NAT COM 2022](https://www.nature.com/articles/s41467-022-34879-1), has been developed to project single cell data into a lower dimensional space while integrating the spatial information into the modelling. Here, we proposed an adaptation of the method for whole slide images (WSIs). To get a low-dimensional representation of these huge images (~20,000 x 20,000 pixels), they are sliced into patches called tiles. For each tile, a vector of features is computed by training a deep learning model; see our [Barlow Twins implementation for WSIs](https://github.com/IARCbioinfo/LNENBarlowTwins/tree/master). These encoded vectors are independent of the tile positions within a WSI. However, we can assume that tiles that are close to each other are more likely to have a similar representation in feature space than distant tiles, as they are more likely to share common morphological features. To model this assumption, we adapted [spatial PCA](https://www.nature.com/articles/s41467-022-34879-1) by removing variable selection and using a multi-sample strategy. Given the quadratic memory and time cost of the algorithm, a random set of vectors must be selected for each patient (~185 tiles per patient), experimentally 50,000 encoded vectors are sufficient to produce a consistent latent space. Intermediate matrices extracted from the SpatialPCA R object created are then used to project new vectors into the low-dimensional space created by the spatial PCA (see supplementary method equation 13  of  [L. Shang and X. Zhou, NAT COM 2022](https://www.nature.com/articles/s41467-022-34879-1)).

- Original article: [L. Shang and X. Zhou, NAT COM 2022](https://www.nature.com/articles/s41467-022-34879-1).
- Original code: [https://github.com/shangll123/SpatialPCA_analysis_codes](https://github.com/shangll123/SpatialPCA_analysis_codes)
- Method used to obtain a spatially informed low-diension represention in "Assessment of the current and emerging criteria for the histopathological classification of lung neuroendocrine tumours in the lungNENomics project." ESMO Open 2023 (under review)

## Installation
- Clone this repository: tested on R 4.1.2
- All needed packages will be install automatocally when the script will be launches
- Please note that the original function of [Spatial PCA package](https://github.com/shangll123/SpatialPCA_analysis_codes) are override by the ones in `ImgSpatialPCA.R` and `ImgSpatialPCAMultipleSamles.R`. 

## Organization of the repository
- `RunMultiSPCARandomSampling.R` allows the Spatial PCA to be run.
- `ImgSpatialPCA.R` contains the based function to create the Spatial PCA and overrides `CreateSpatialPCAObject` function of the original package.
- `ImgSpatialPCAMultipleSamles.R` adapted the Spatial PCA to several samples and overrides `SpatialPCA_Multiple_Sample` function of the original package.

## Step 1: Creation of the Spatial PCA latent space  
- To create a spatial PCA R object run `RunMultiSPCARandomSampling.R` an example of configuration file is given in `RunSpatialPCA50K.sh`
- Command line for cluster running with slurm
```
sbatch RunSpatialPCA50K.sh
```
### Description of the process

1. Load encoded vectors, those one has to be concatenated in a single csv file such as (see: `path2projectors`):

|   | X0           | X1            | X2          | X3            | ... | X124         | X125         | X126         | X127         | img_id                  | sample_id   | img_id_c               | x           | y           |
| - | ------------ | ------------- | ----------- | ------------- | --- | ------------ | ------------ | ------------ | ------------ | ------------------------ | -------- | ---------------------- | ----------- | ----------- |
| 1 | 0.010731053  | -0.017491885  | -0.05379057 | 0.0060576447  | ... | -0.021526879 | 0.038895514 | 0.021861676 | -0.0008289963 | TNE1019_30721_19585     | TNE1019  | TNE1019_30721_19585   | 30721       | 19585       |
| 2 | 0.0031735892 | -0.0024470983 | -0.04042089 | 7.895916e-05  | ... | -0.01900657  | -0.0067212125| 0.0070669674| -0.015635846  | TNE1019_33409_28801     | TNE1019  | TNE1019_33409_28801   | 33409       | 28801       |

2. Extraction of n random row of in the data frame (n = `n_tiles`)
3. Creation of list of tables of features and coordinates per samples
4. Creation of the Spatial PCA
5. Save the SpatialPCA R object and coordinates in `output_folder`

- :warning: WARNING :warning:
    - For a representation containing 100 encoded vectors, a machine with 300 GB of RAM is required, and the R object that is created has a size of 6 GB. 
    - The encoded vectors must not be normalised, this step is included in the pipeline.

## TO DO LIST
- :construction: Projections
- :construction: Leiden community 
- :construction: Random forest