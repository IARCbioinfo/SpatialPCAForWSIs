print("#-------------------------------------- Install package ------------------------------------------------------------------------")
#libPaths(c(tempdir(), .libPaths()))
print(.libPaths())
.libPaths("/home/mathiane/R/x86_64-pc-linux-gnu-library/4.0")
packages <- c( "Matrix","optparse", "ggplot2" )
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
print("# --------------------------------------- Load libraies ------------------------------------------------------------------------")
library(stringr)
library(Matrix)
library(RSpectra)
library("optparse")
# To comment
library(ggplot2)
library(gridExtra)
setwd("/home/mathiane/LNENWork/SpatialPCAForWSIs")


source('ImgSpatialPCAMultipleSamles.R')
environment(SpatialPCA_Multiple_Sample) <- asNamespace('SpatialPCA')
assignInNamespace("SpatialPCA_Multiple_Sample", SpatialPCA_Multiple_Sample, ns = "SpatialPCA")


source('ImgSpatialPCA.R')
environment(CreateSpatialPCAObject) <- asNamespace('SpatialPCA')
assignInNamespace("CreateSpatialPCAObject", CreateSpatialPCAObject, ns = "SpatialPCA")

print("# --------------------------------- Command line arguments ------------------------------------------")

option_list = list(
  make_option(c("--sample_id"), type="character",
              help="Patient ID matching WSI filename", metavar="character"),
  make_option(c("--spca_obj"), type="character",
              help="Path to the SpatialPCA R object", metavar="character"),
  make_option(c("--dimension_enc_vectors"), type="integer", default=128, 
        help="Dimension of encoded vectors"),
  make_option(c("--proj_tab_norm"), type="character",
              help="Access path to the projectors table. Projectors must be centred and standardised.", metavar="character"),
  make_option(c("--outdir"), type="character",
              help="Folder where the results will be store", metavar="character")
  
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

sampleID =opt$sample_id 
SPCA_r_obj =opt$spca_obj 
d_enc_vectors = opt$dimension_enc_vectors
path_projectors_norm = opt$proj_tab_norm
outdir = opt$outdir

dir.create(outdir, showWarnings = FALSE)
exp_folder = outdir
# Output filename
nfilename = paste(exp_folder, paste('Proj', paste( sampleID , '.csv', sep = '') , sep = '_'), sep = '/')
if ( (file.exists(nfilename) == FALSE) & (sampleID != "TNE233")){


print("# -------------------- Load data -------------------------")
load(SPCA_r_obj)


print("------------------------------------Extract information form the original spatial PCA ---------------------")
W <- LIBD_multi$MultipleSample_SpatialPCA_object@W
tau  <- LIBD_multi$MultipleSample_SpatialPCA_object@tau
X <- LIBD_multi$MultipleSample_SpatialPCA_object@params$X

print("Spatial PCA training object loaded")
SpatialPCs_ori <- LIBD_multi$MultipleSample_SpatialPCA_object@SpatialPCs
SpatialPCs_ori_t <- t(SpatialPCs_ori)
head(SpatialPCs_ori_t)
rownames(SpatialPCs_ori_t) <- unlist(lapply(rownames(SpatialPCs_ori_t), function(x)  substr(x, str_locate(x, '_')+1, nchar(x))  ))
# ecoded vector
df_enc_vector <- read.csv(path_projectors_norm)

df_enc_vector$x <- as.integer(df_enc_vector$x)
df_enc_vector$y <- as.integer(df_enc_vector$y)


df_enc_vector_c_sampleid = df_enc_vector[df_enc_vector$tne_id == sampleID,]
head(df_enc_vector_c_sampleid)


df_enc_vector_c_sampleid = df_enc_vector_c_sampleid[!duplicated(df_enc_vector_c_sampleid$img_id_c ),]


df_coords_csampleid <- data.frame('x' = df_enc_vector_c_sampleid$x)
df_coords_csampleid$y <-  df_enc_vector_c_sampleid$y
df_coords_csampleid <- scale(df_coords_csampleid)
bandwidth = LIBD_multi$MultipleSample_SpatialPCA_object@bandwidth

rm(LIBD_multi)
gc(verbose=T)
print("Spatial PCA training object deleted")

K_csampleid = exp(-1*as.matrix(dist(df_coords_csampleid)^2)/bandwidth)
print(dim(K_csampleid))

# SVD of K
n = dim(K_csampleid)[1]
fast_eigen_num = ceiling(n*0.1)
EIGEN = eigs_sym(K_csampleid, k=fast_eigen_num, which = "LM")
U=EIGEN$vectors
delta=EIGEN$values

rm(EIGEN)
gc(verbose=T)
print("SVD computed")


enc_val_for_pca = df_enc_vector_c_sampleid[,2:(d_enc_vectors+1)]
enc_val_t_for_pca <- t(enc_val_for_pca)
enc_val_t_for_pca <- Matrix(nrow = nrow(enc_val_t_for_pca), ncol = ncol(enc_val_t_for_pca), data = enc_val_t_for_pca, sparse = TRUE)
colnames(enc_val_t_for_pca) <- df_enc_vector_c_sampleid$img_id_c
rownames(enc_val_t_for_pca) <- colnames(enc_val_for_pca)[2:(d_enc_vectors+1)]
Y <- enc_val_t_for_pca
rownames(df_enc_vector_c_sampleid) = df_enc_vector_c_sampleid$img_id_c
coords_df_c_for_PCA = as.matrix(df_enc_vector_c_sampleid)


print("--- Compute M ---")
H = matrix(1, dim(Y)[2],1)
HH_inv=solve(t(H)%*%H,tol = 1e-40)
HH = H%*%HH_inv%*%t(H)
M=diag(dim(Y)[2])-HH
# Y=expr
q=1

rm(H)
rm(HH)
rm(HH_inv)
gc(verbose=T)
print("M computed")


### New projection

W_hat_t = t(W)
YM <- Y %*% M
rm(Y)
WtYM = W_hat_t%*%YM
WtYMK = WtYM %*% K_csampleid
WtYMU = WtYM %*% U
rm(WtYM)
Ut=t(U)
UtM = Ut %*% M
rm(YM)
rm(M)
gc(verbose=T)
print("WtYMK and WtYMU computed")


UtMU = UtM %*% U
middle_inv = solve(1/tau * diag(1/delta) + UtMU, tol = 1e-40)
rm(UtMU)
UtMK = UtM %*% K_csampleid
rm(UtM)
gc(verbose=T)
print("middle_inv computed")

SpatialPCs_proj = tau*WtYMK - tau*WtYMU %*% middle_inv %*% UtMK
SpatialPCs_proj_t <- t(SpatialPCs_proj)


colnames(SpatialPCs_proj_t)<-unlist(lapply(1:20, function(x) paste('axis', x, sep='_')))
SpatialPCs_proj_t <- data.frame(SpatialPCs_proj_t)
img_id <- df_enc_vector_c_sampleid$img_id_c
rownames(SpatialPCs_proj_t) <- img_id
write.csv(SpatialPCs_proj_t, nfilename)

}
