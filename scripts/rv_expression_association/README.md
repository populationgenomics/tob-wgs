## Outlier analysis

### Generate positive association examples
The first step for this analysis is to construct a set of "true associations" by identifying expression outliers (e.g., individuals that have extreme - outlier - expression for a specific gene) and then identifying genetic variants within the vicinity of that gene for which the outlier individual(s) have alternative alleles (both 0/1 and 1/1).

[This script](../rv_expression_association/get_rv_outliers.py) takes a gene and an individual IDs as inputs and returns the following file:

```donor ID | gene ID | variant ID | position | CADD | MAF (OneK1K) | MAF (gnomad) | regulatory consequences (VEP) | ...```

Variants are only included if they are also annotated to have regulatory consequences (from VEP).
CADD scores and MAF are annotated but not filtered for.

### Convert to plink format
In order to test for associations between a gene's expression and a set of variants (using CellRegMap, or SAIGE-GENE), the specific set of variants needs to be converted back to plink format files (bed, bim, fam).
