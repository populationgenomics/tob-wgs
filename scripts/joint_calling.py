#!/usr/bin/env python

"""
Drive the joint calling workflow
"""

import subprocess
import sys

subprocess.run(
    'wget https://raw.githubusercontent.com/populationgenomics/joint-calling/develop/workflows/batch_workflow.py -O batch_workflow.py',
    shell=True,
    check=False,
)
subprocess.run('mkdir joint_calling', shell=True, check=False)
subprocess.run(
    'wget https://raw.githubusercontent.com/populationgenomics/joint-calling/develop/joint_calling/utils.py -O joint_calling/utils.py',
    shell=True,
    check=False,
)
is_test = len(sys.argv) > 1 and sys.argv[1] == 'test'
subprocess.run(
    'python batch_workflow.py '
    + '--callset tob-wgs '
    + '--version v0 '
    + ('--batch 0 ' if is_test else '--batch 0 --batch 1 ')
    + '--from '
    + ('test ' if is_test else ' main')
    + '--to '
    + ('temporary ' if is_test else ' analysis')
    + '--keep-scratch '
    + '--billing-project tob-wgs ',
    shell=True,
    check=False,
)
