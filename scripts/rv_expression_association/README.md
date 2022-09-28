## Analysis plan

This runs a gene-set association analysis to test for an association between i) a set of rare genetic variants located in and around (now, up to 50kb up and downstream of) a gene and ii) the expression level of the gene itself (in a given cell type).

Testing is done on one gene within chromosome 22, _VPREB3_, in naive B cells, but the analysis will be later expanded to multiple genes and other cell types.

### Overview
This folder contains three scripts:
* [get_vep_variants.py](get_vep_variants.py) is a Python script which takes a hail MT + HT object, filters for relevant variants (detailed below), and exports to plink files
* [prepare_inputs.py](prepare_inputs.py) is a Python script which takes in genotype and expression files and prepares the input files to run SKAT
* [run_SKAT.R](run_SKAT.R) is an R script that runs SKAT using the inputs generated from the previous scripts and returns association p-values

### Step 1 - select genetic variants

At the moment, I have subset both the MT object and the [VEP](https://asia.ensembl.org/info/docs/tools/vep/index.html)-annotated object manually, using [this script](https://github.com/populationgenomics/analysis-runner/blob/main/scripts/subset_matrix_table.py) to a specific genomic region around this one specific gene (_VPREB3_), _e.g._, for the VEP-annotated hail table (from local version of populationgenomics/analysis-runner/scripts/):
```
analysis-runner --dataset tob-wgs \
    --description "subset vep annotated ht" \
    --output-dir "v0" --access-level standard \
    python3 subset_hail_table.py -i gs://cpg-tob-wgs-main/tob_wgs_vep/104/vep104.3_GRCh38.ht \
    --chr chr22 --pos 23702743-23804425 --out VPREB3_50K_window_vep
```
To get the interval, for now in a notebook:
```
import hail as hl
import pandas as pd

gene_file = 'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/gene_location_files/GRCh38_geneloc_chr22.tsv'
gene_df = pd.read_csv(gene_file, sep='\t', index_col=0)

chrom = 22
gene_name = 'VPREB3'
window_size = 50000

interval_start = int(gene_df[gene_df['gene_name'] == gene_name]['start']) - int(window_size)
interval_end = int(gene_df[gene_df['gene_name'] == gene_name]['end']) + int(window_size)

# clip to chromosome boundaries
left_boundary = max(1, interval_start)
right_boundary = min(interval_end, hl.get_reference('GRCh38').lengths[f'chr{chrom}'])

# get gene-specific genomic interval
gene_interval = f'chr{chrom}:{left_boundary}-{right_boundary}'
gene_interval
```
In the future, either i) add a previous step subsetting the MT + HT objects to the right genomic region taking a gene name as input, or ii) add that step to this script.

This step selects QC-passing, biallelic single nucleotide variants (SNVs) that are rare (alternative allele frequency < 5%) and that are predicted by VEP to have regulatory consequences (_e.g._, in promoters, enhancers, transcription factor binding sites..).
Then, it creates [PLINK](https://zzz.bwh.harvard.edu/plink/) files (.bed, .bim, .fam) for those variants only.
To run this, use:
```
analysis-runner --dataset "tob-wgs" \
    --description "get set of variants for a gene, convert to plink" \
    --access-level "test" \
    --output-dir "v0" \
    get_vep_variants.py
```

### Step 2 - prepare input files for SKAT

[SKAT](https://github.com/leelabsg/SKAT) is an R package to run gene-set association tests.
It requires as inputs:
* a genotype matrix (Z), samples X genetic variants
* a phenotype vector (y.c), samples X 1 (in my case, this will be the expression level of a gene across samples, see below)
* (optionally,) a kinship matrix (K), samples X samples
* (optionally,) a covariance matrix (X), samples X covariates

I use this (Python) script to prepare these input files, using similar steps as in the [CellRegMap pipeline](https://github.com/populationgenomics/cellregmap-pipeline).
To run:
```
analysis-runner --dataset "tob-wgs" \
    --description "prepare SKAT input files" \
    --access-level "test" \
    --output-dir "v0" \
    prepare_inputs.py
```

#### current steps

For now, considering all regulatory variants (or, _e.g.,_ selecting promoter variants only) and considering pseudo-bulk (mean) expression for each individual, borrowing the expression file (here, donor-level aggregated expression in naive B cells) from Kat's eQTL files.

### Step 3 - run SKAT

This is an R script that run SKAT-O in multiple modes (after matching sample IDs across objects):
* with and without the kinship matrix K
* reporting the p-values for each of the burden, SKAT and optimised SKAT-O tests. (print only atm, change to saving p-values to file)

To run:
```
analysis-runner --dataset "tob-wgs" \
    --description "run SKAT" \
    --access-level "test" \
    --output-dir "v0" \
    --image australia-southeast1-docker.pkg.dev/analysis-runner/images/driver-r:1.2 \
    run_SKAT.R
``` 
