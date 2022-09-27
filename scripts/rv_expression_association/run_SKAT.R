install.packages("SKAT")

library(googleCloudStorageR)
library(gargle)
library(SKAT)

# token authorisation
scope <- c("https://www.googleapis.com/auth/cloud-platform")
token <- token_fetch(scopes = scope)
gcs_auth(token = token)

# set bucket
googleCloudStorageR::gcs_global_bucket("gs://cpg-tob-wgs-test")

# # get genotype files (plink format)
# File.Bed <- googleCloudStorageR::gcs_get_object("v0/plink_files/igll5_rare_regulatory.bed")
# File.Bim <- googleCloudStorageR::gcs_get_object("v0/plink_files/igll5_rare_regulatory.bim")
# File.Fam <- googleCloudStorageR::gcs_get_object("plink_files/igll5_rare_regulatory.fam")

# File.SetID <- googleCloudStorageR::gcs_get_object("v0/skat/igll5_rare_regulatory.SetID")

# File.SSD <- "IGLL5_50K_RV_VEP.SSD"
# File.Info <- "IGLL5_50K_RV_VEP.Info"

# # generate SSD file
# ssd <- Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info)

# # get expression file (pseudobulk)
# E <- 

# # get covariates
# X <- 1  # intercept only

# # set up SKAT
# Z <- Get_Genotypes_SSD(ssd)
# y.c <- 


# get input files generated elsewhere
# Z <- 
# y.c <- 
# K <- 
X <- 1  # intercept only

# run null model (no Kinship)
obj <- SKAT_Null_Model(y.c ~ X, out_type = "C")

# get all three p-values
pv_skat <- SKAT(Z, obj)$p.value  # SKAT
pv_burden <- SKAT(Z, obj, r.corr = 1)$p.value  # burden
pv_skat_o <- SKAT(Z, obj, method="SKATO")$p.value  # SKAT-O

print(paste0("no Kinship p-values, SKAT: ", pv_skat, ", burden: ", pv_burden, ", SKAT-O: ", pv_skat_o))

# run null model (with Kinship matrix)
obj_K <- SKAT_NULL_emmaX(y ~ X, K = K, out_type = "C")

# get all three p-values
pv_skat <- SKAT(Z, obj_K)$p.value  # SKAT
pv_burden <- SKAT(Z, obj_K, r.corr = 1)$p.value  # burden
pv_skat_o <- SKAT(Z, obj_K, method="SKATO")$p.value  # SKAT-O

print(paste0("with Kinship p-values, SKAT: ", pv_skat, ", burden: ", pv_burden, ", SKAT-O: ", pv_skat_o))