### Analysis plan

This folder contains three scripts:
* igll5_vep.py is a Python script which takes a MT + HT object, filters relevant variants (detail below) and exports to plink files
* prepare_inputs.py is a Python script which takes in genotype and expression files and prepares the input files to run SKAT
* run_SKAT.R is an R script that runs SKAT using the inputs generated from the previous script