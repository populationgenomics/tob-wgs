#!/usr/bin/env python

"""
Helper to submit scripts to dataproc
"""

import sys
import subprocess

subprocess.run(sys.argv[1:], check=False)
