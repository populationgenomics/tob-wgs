#!/usr/bin/env python3

import hail as hl
from hail.methods import export_plink
# import numpy as np
# import pandas as pd
from cpg_utils.hail_batch import dataset_path, init_batch, reference_patch
# from cloudpathlib import AnyPath

MT = dataset_path('v0/IGLL5_50K_window.mt')

init_batch()
mt = hl.read_matrix_table(MT)
mt = hl.experimental.densify(mt)

vep_ht = reference_path('tob_wgs_vep/104/vep104.3_GRCh38.ht')
mt = mt.annotate_rows(vep = vep_ht[mt.row_key].vep)

export_plink(mt, 'plink_files', ind_id = mt.s)