#!/usr/bin/env python3

import hail as hl
import numpy as np
# import pandas as pd
from cpg_utils.hail_batch import dataset_path, init_batch
from cloudpathlib import AnyPath

# MT = dataset_path('mt/v7.mt')
MT = dataset_path('v0/IGLL5_50K_window.mt')

init_batch()
mt = hl.read_matrix_table(MT)
# mt = mt.filter_rows(mt.locus.contig == 'chr22')
mt = hl.experimental.densify(mt)
print(mt.count())

mt = mt.filter_rows(hl.len(hl.or_else(mt.filters, hl.empty_set(hl.tstr))) == 0)
mt = hl.variant_qc(mt)

mt = mt.filter_rows(hl.len(mt.alleles) == 2)
mt = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1]))

# gene_file = (
#     'scrna-seq/grch38_association_files/gene_location_files/GRCh38_geneloc_chr22.tsv'
# )
# gene_df = pd.read_csv(AnyPath(dataset_path(gene_file)), sep='\t', index_col=0)

# gene_name = 'IGLL5'
# window_size = 50000

# interval_start = int(gene_df[gene_df['gene_name'] == gene_name]['start']) - int(
#     window_size
# )
# interval_end = int(gene_df[gene_df['gene_name'] == gene_name]['end']) + int(window_size)

# chrom = 22

# # clip to chromosome boundaries
# left_boundary = max(1, interval_start)
# right_boundary = min(interval_end, hl.get_reference('GRCh38').lengths[f'chr{chrom}'])
# gene_interval = f'chr{chrom}:{left_boundary}-{right_boundary}'

# mt = hl.filter_intervals(
#     mt,
#     [hl.parse_locus_interval(gene_interval, reference_genome='GRCh38')],
# )

print(mt.count())

sample = 'CPG9951'
donor_mt = mt.filter_cols(mt.s == sample)
donor_mt = donor_mt.filter_rows(hl.agg.any(donor_mt.GT.is_non_ref()))

print(np.nanmin(donor_mt.variant_qc.AF[1].collect()))

donor_mt05 = donor_mt.filter_rows(donor_mt.variant_qc.AF[1] < 0.05)
donor_mt01 = donor_mt.filter_rows(donor_mt.variant_qc.AF[1] < 0.01)

print(donor_mt05.count())
print(donor_mt01.count())
print(donor_mt05.variant_qc.AF.show())