## Parameters and settings affecting gene-set association tests power

These scripts aim to assess what parameters most influence power to detect effects using a gene-set association test.
In a first instance, using the three tests implemented in the [SKAT R package](https://cran.r-project.org/web/packages/SKAT/SKAT.pdf) + [ACAT](https://github.com/yaowuliu/ACAT):
* SKAT: SNP-Set (Sequence) Kernel Association Test ([Wu et al, AJHG 2011](https://www.sciencedirect.com/science/article/pii/S0002929711002229))
* burden: standard burden test
* SKAT-O: "omnibus" optimised combined test ([Lee et al, AJHG 2012](https://www.sciencedirect.com/science/article/pii/S0002929712003163)).
* ACAT-V: gene-set test using the aggregated Cauchy association test (ACAT; [Liu et al, AJHG 2019](https://www.sciencedirect.com/science/article/pii/S0002929719300023))

Later, this will be extended to other methods, including my own [CellRegMap-RV (new name coming soon)](https://github.com/annacuomo/CellRegMap/blob/main/cellregmap/_cellregmap.py#L657-L687). 

### Approach
* Genotypes matrix (```genotypes```): real variants (SNVs from a genomic region, around the _VPREB3_ gene) from the TOB-WGS dataset
  * samples: subsetting to 100 and then 1,000 individuals
  * variants: singletons only (so freq=0.005 and 0.0005, respectively)
* Phenotype vector pheno as ```pheno = genotypes*beta + noise``` 
  * where noise is random and
    * normally distributed with mean=0, sd=1, or
    * following a Poisson distribution with lambda=1
  * genotypes and beta vary as below

#### Scenarios
* only testing 10 variants, all causal (beta !=0), same magnitude and direction of effect (beta=1)
* testing more variants (20, 50), but where only 10 are causal, same magnitude and direction of effect (beta=1)
* varying direction of effect: only testing the 10 causal variants, but for 2 (or 5) variants beta=-1 (beta=1 for the remaining 8 (or 5) variants)
* varying magnitude of effect: only testing the 10 causal variants, but varying magnitude of effect (beta=[0.1, 0.2, ..., 1])

I am also recording how "normal" y is in each case (using the [Shapiro test](https://en.wikipedia.org/wiki/Shapiro%E2%80%93Wilk_test) p-value as measure).

#### Future
Other aspects I would like to include are:
* frequency (not just singletons)
* getting y (without noise) through a g-link function, _e.g._ to look more Poisson (for now, I have added Poisson noise)
* more cells per donor (while I am not explicitly saying so, here I am modelling only one cell/observation per donor) - SKAT should be able to handle these repeats if I add the repeat structure as background GRM, I hope, but this may be where CellRegMap-RV should look better?

Even further in the future is to build more complex simulations drawing from real scRNA-seq data.

### Run scripts
To run (after authenticating: ```gcloud auth application-default login```).

R scripts (SKAT)
```
analysis-runner --dataset "tob-wgs" \
    --description "test SKAT" \
    --access-level "test" \
    --output-dir "v0" \
    --image australia-southeast1-docker.pkg.dev/analysis-runner/images/driver-r:1.3 \
    test_SKAT.R
```

python scripts (CellRegMap)
```
analysis-runner --dataset "tob-wgs" \
    --description "test CRM-RV" \
    --access-level "test" \
    --output-dir "v0" \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/cellregmap:0.0.3 \
    test_CellRegMap_RV.py
```
