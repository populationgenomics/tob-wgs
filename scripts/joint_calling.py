#!/usr/bin/env python

"""
Drive the joint calling workflow
"""

import subprocess

subprocess.run(
    'batch_workflow.py '
    '--callset tog-wgs '
    '--version v0 '
    '--output_bucket gs://cpg-tob-wgs-temporary/v0/ '
    '--keep_scratch',
    shell=True,
    check=False,
)
