## Outlier analysis

### Generate positive association examples
The first step for this analysis is to construct a set of "true associations" by identifying expression outliers (e.g., individuals that have extreme - outlier - expression for a specific gene) and then identifying genetic variants within the vicinity of that gene for which the outlier individual(s) have alternative alleles (both 0/1 and 1/1).

[This script](../rv_expression_association/get_rv_outliers.py) takes a gene and an individual IDs as inputs and returns the following file:

```donor ID | gene ID | variant ID | position | CADD | MAF (OneK1K) | MAF (gnomad) | regulatory consequences (VEP) | ...```

Variants are only included if they are also annotated to have regulatory consequences (from VEP).
CADD scores and MAF are annotated but not filtered for.

### Test them for association

#### Convert genotypes to plink format
In order to test for associations between a gene's expression and a set of variants (using CellRegMap, or SAIGE-GENE), the specific sets of variants need (for now?) to be converted back to plink format files (bed, bim, fam).

[This script]() does that for a specific gene using the file generated above.

#### Run association 
TO DO

## Systematic analysis

For all genes (with some criterion, e.g. min expression) create a similar file as above, this time without donor annotations (and therefore without filter for variants with alt allele for that one specific individual:

```gene ID | variant ID | position | CADD | MAF (OneK1K) | MAF (gnomad) | regulatory consequences (VEP) | ...```

Then generate input files (including converting to plink etc) for whatever association tool we have established as optimal.
