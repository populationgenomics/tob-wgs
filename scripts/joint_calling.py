#!/usr/bin/env python

"""
Drive the joint calling workflow
"""

import subprocess

subprocess.run(
    'wget https://raw.githubusercontent.com/populationgenomics/joint-calling/develop/workflows/batch_workflow.py',
    shell=True,
    check=False,
)
subprocess.run('mkdir joint_calling', shell=True, check=False)
subprocess.run(
    'wget https://raw.githubusercontent.com/populationgenomics/joint-calling/develop/joint_calling/utils.py joint_calling',
    shell=True,
    check=False,
)
subprocess.run(
    'python batch_workflow.py '
    '--callset tog-wgs '
    '--version v1 '
    '--batch 0 --batch 1 '
    '--from main '
    '--to temporary '
    '--keep_scratch '
    '--billing-project tob-wgs-standard ',
    shell=True,
    check=False,
)
