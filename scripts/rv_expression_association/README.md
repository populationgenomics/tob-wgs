## Analysis plan

The overall goal of these scripts is to run a gene-set association test to test for an association between i) a number of rare genetic variants located around a gene and ii) the expression level of the gene itself (in a given cell type).

At first, I am testing specifically one chromosome 22 gene (_IGLL5_), in naive B cells.
I consider variants located at most 50kb up or down-stream of the gene body (including those within the gene itself).
Later, I'll want to run this for many genes.

### Overview
This folder contains three scripts:
* [igll5_vep.py](igll5_vep.py) is a Python script which takes a MT + HT object, filters relevant variants (detail below) and exports to plink files
* [prepare_inputs.py](prepare_inputs.py) is a Python script which takes in genotype and expression files and prepares the input files to run SKAT
* [run_SKAT.R](run_SKAT.R) is an R script that runs SKAT using the inputs generated from the previous script

### Step 1 - select genetic variants

At the moment, I have subsetted both the MT object and the [VEP](https://asia.ensembl.org/info/docs/tools/vep/index.html)-annotated object manually, using [this script](https://github.com/populationgenomics/analysis-runner/blob/main/scripts/subset_matrix_table.py) to a specific genomic region around this one specific gene (IGLL5), _e.g._, for the VEP-annotated hail table (from local version of populationgenomics/analysis-runner/scripts/):
```
analysis-runner --dataset tob-wgs \
    --description "subset vep annotated ht" \
    --output-dir "v0" --access-level standard \
    python3 subset_hail_table.py -i gs://cpg-tob-wgs-main/tob_wgs_vep/104/vep104.3_GRCh38.ht \
    --chr chr22 --pos 22837780-22946111 --out IGLL5_50K_window_vep
```
To get the interval, for now in a notebook:
```
import hail as hl
import pandas as pd

gene_file = 'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/gene_location_files/GRCh38_geneloc_chr22.tsv'
gene_df = pd.read_csv(gene_file, sep='\t', index_col=0)

chrom = 22
gene_name = 'IGLL5'
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
In the future, either 1) add a previous step subsetting the MT + HT objects to the right genomic region taking a gene name as input, or add that to this script.

This step selects QC-passing, biallelic SNP that are rare (alternative allele frequency < 5%) and that are predicted by VEP to have regulatory consequences.
Then, it creates plink files (.bed, .bim, .fam) for those variants only.
To run this, use:
```
analysis-runner --dataset "tob-wgs" \
    --description "get set of variants for a gene, convert to plink" \
    --access-level "test" \
    --output-dir "v0" \
    igll5_vep.p
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
    --description "open adata using scanpy" \
    --access-level "test" \
    --output-dir "v0" \
    prepare_inputs.py
```

### Step 3 - run SKAT

This is an R script that run SKAT-O in multiple modes:
* with and without the kinship matrix
* reporting the p-values for each of the burden, SKAT and optimised SKAT-O tests.

To run (??):
```Rscript run_SKAT.R``` 
