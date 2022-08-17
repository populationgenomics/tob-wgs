## Outlier analysis

### Generate positive association examples
The first step for this analysis is to construct a set of "true associations" by identifying expression outliers (_e.g.,_ individuals that have extreme - outlier - expression for a specific gene) and then identifying genetic variants within the vicinity of that gene for which the outlier individual(s) have alternative alleles (either 0/1 or 1/1).

[This script](../rv_expression_association/get_rv_outliers.py) takes a gene and an individual IDs as inputs and returns the following file:

```donor ID | gene ID | variant ID | position | CADD | MAF (OneK1K) | MAF (gnomad) | regulatory consequences (VEP) | ...```

Variants are only included if they are also annotated to have regulatory consequences (from VEP).
CADD scores and MAF are annotated but not filtered for.

To run this script using the analysis runner, run:
```
analysis-runner --dataset "tob-wgs" \
    --description "get set of variants for given outlier" \
    --access-level "test" \
    --output-dir "tob_wgs_rv/expression_outliers" \
    get_rv_outliers.py --onek1k-id "943_944" --gene-name "IGLL5" --chrom "22" 
```

### Test them for association

#### Convert genotypes to plink format
In order to test for associations between a gene's expression and a set of variants (using an updated [CellRegMap](https://github.com/limix/CellRegMap/), or [SAIGE-GENE](https://saigegit.github.io//SAIGE-doc/docs/set.html), or other), the specific sets of variants need (for now?) to be converted back to plink format files (bed, bim, fam).

[This future script]() does that for a specific gene using the file generated above.
TO DO

#### Run association 
TO DO

## Systematic analysis

For all genes (with some criterion, _e.g.,_ min expression) create a similar file as above, this time without donor annotations (and therefore without filter for variants with alt allele for that one specific individual:

```gene ID | variant ID | position | CADD | MAF (OneK1K) | MAF (gnomad) | regulatory consequences (VEP) | ...```

Then generate input files (including converting to plink etc) for whatever association tool we have established as optimal.
