#!/usr/bin/env python3
# pylint: disable=wrong-import-position

# This script runs our CellRegMap-RV method to test for an association between a set of 
# rare genetic variants (within the promoter of a given gene)
# and the expression of the gene itself (in one cell type, and aggregated across donors). 
# The tests performed are a variance-component test, a burden test and an omnibus test
# The resulting association p-values are recorded.

# import python modules
import sys
import subprocess

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
