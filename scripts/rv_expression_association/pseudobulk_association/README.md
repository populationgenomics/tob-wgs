## Pseudobulk rare variants association analysis

Scripts to find associations between rare variants (from WGS data) and single-cell expression (from scRNA-seq data) aggregated as "pseudobulk".

As a first step to evaluate our ability to test for association between rare variation (freq<5%) and single-cell expression, we consider aggregated expression across donors, separately for each cell type.
"Pseudobulk" expression (in this case mean expression for each gene and donor) resembles more closely "bulk" expression data, which is less sparse than single-cell data and thus is easier to model.
Additionally, it results in a single observation per individual (as opposed to multiple cells) which also makes it easier to deal with.
Thus, we perform this analysis as a benchmark method.

### Step 1: obtain variants

Script to get genotype data for variants that are:
* in and around a given gene (at a given window size)
* in promoter region (as annotated by VEP)
* biallelic SNVs
* rare (freq<5%)

We save the object as a hail table for downstream analyses, and export to plink for association testing.

To run:
```
analysis-runner --dataset "tob-wgs" \
    --description "get set of promoter variants for a gene, convert to plink" \
    --access-level "standard" \
    --output-dir "tob_wgs_rv/plink" \
    get_promoter_variants.py
```

### Step 2: prepare input files
Next, we prepare the input files to perform the association (Step 3)