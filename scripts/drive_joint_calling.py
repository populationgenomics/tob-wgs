#!/usr/bin/env python

"""
Drive the joint calling workflow
"""

import subprocess
import sys


def run_cmd(cmd):
    """Print the command and run"""
    print(cmd)
    subprocess.run(
        cmd,
        shell=True,
        check=False,
    )


is_test = len(sys.argv) > 1 and sys.argv[1] == 'test'
version = sys.argv[2] if len(sys.argv) > 2 else 'v0'

run_cmd(
    'PYTHONPATH=$PWD/../joint-calling python ../joint-calling/workflows/batch_workflow.py '
    + '--callset tob-wgs '
    + f'--version {version} '
    + ('--batch 0 ' if is_test else '--batch 0 --batch 1 ')
    + '--from '
    + ('test ' if is_test else 'main ')
    + '--to '
    + ('temporary ' if is_test else 'analysis ')
    + '--keep-scratch '
    + '--billing-project tob-wgs '
)
