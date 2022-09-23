library(googleCloudStorageR)
library(gargle)
library(SKAT)

# token authorisation
scope <- c("https://www.googleapis.com/auth/cloud-platform")
token <- token_fetch(scopes = scope)
gcs_auth(token = token)

# set bucket
googleCloudStorageR::gcs_global_bucket("gs://cpg-tob-wgs-test")

# get genotype files (plink format)
G <- 

# get expression file (pseudobulk)
E <- 

# get covariates
X <- 

# set up SKAT
Z <-  G
y.c <- 

# run null model
obj <- SKAT_Null_Model(y.c ~ X, out_type = "C")


# get all three p-values
pv_skat <- SKAT(Z, obj)$p.value  # SKAT
pv_burden <- SKAT(Z, obj, r.corr = 1)$p.value  # burden
pv_skat_o <- SKAT(Z, obj, method="SKATO")$p.value  # SKAT-O

print(paste0("SKAT: ", pv_skat, ", burden: ", pv_burden, ", SKAT-O: ", pv_skat_o))