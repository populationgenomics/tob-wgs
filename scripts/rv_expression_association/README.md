This script counts how many rare variants (freq<5%) in the TOB-WGS dataset there are.
The filters applied are:
* QC
* biallelic SNVs
* not HomRef
* freq<5%

To run:
```
analysis-runner --dataset "tob-wgs" \
    --description "count TOB rare variants" \
    --access-level "standard" \
    --output-dir "v0" \
    count_variants.py
```
