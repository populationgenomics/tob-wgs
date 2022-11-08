This script counts the number of rare variants (freq<5%) present in a Hail Matrix object provided as input.
In addition to filtering for variants that are rare, the following additional filters are applied:
* Variants that pass GATK QC filters
* biallelic SNVs (single nucleotide changes, with only two possible alleles)
* not ref-ref (i.e. variable within the given dataset)

To run:
```
analysis-runner --dataset "tob-wgs" \
    --description "count TOB rare variants" \
    --access-level "standard" \
    --output-dir "v0" \
    count_variants.py --mt-path mt/v7.mt
```
