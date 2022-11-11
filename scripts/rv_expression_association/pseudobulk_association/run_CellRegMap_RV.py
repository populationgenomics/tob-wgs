#!/usr/bin/env python3
# pylint: disable=wrong-import-position,missing-module-docstring,no-value-for-parameter,import-error

# This script runs our CellRegMap-RV method to test for an association between a set of 
# rare genetic variants (within the promoter of a given gene)
# and the expression of the gene itself (in one cell type, and aggregated across donors). 
# The tests performed are a variance-component test, a burden test (max) and an omnibus test
# The resulting association p-values are recorded.

# import python modules
import sys
import subprocess
import click
import numpy as np
import pandas as pd

# install CellRegMap (new version) from github
subprocess.run(
    [
        sys.executable,
        '-m',
        'pip',
        'install',
        'git+https://github.com/annacuomo/CellRegMap',
        '--force-reinstall',  # install github version and overwrite current
        '--no-dependencies',  # same dependencies, no need to uninstall and reinstall those
    ],
    check=True,
)

from cellregmap import (
    run_gene_set_association,
    run_burden_association,
    omnibus_set_association,
)

def get_crm_pvs(pheno, covs, genotypes, contexts=None):
    """
    CellRegMap-RV tests
    * score test (variance)
    * burden test (max)
    * omnibus (Cauchy) test
    """
    pv0 = run_gene_set_association(y=pheno, G=genotypes, W=covs, E=contexts)[0]
    pv1 = run_burden_association(y=pheno, G=genotypes, W=covs, E=contexts, mask="mask.max")[0]
    pv2 = omnibus_set_association(np.array([pv0, pv1]))
    return [pv0, pv1, pv2]


@click.command()
@click.option('--pheno-file', required=True) 
@click.option('--geno-file', required=True) 
def main(
    pheno_file: str,
    geno_file: str,
):  
    pheno = pd.read_csv(pheno_file)
    geno = pd.read_csv(geno_file)


if __name__ == '__main__':
    main()

