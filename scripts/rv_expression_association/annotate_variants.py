#!/usr/bin/env python3

import click
import hail as hl
import pandas as pd
from cpg_utils.hail_batch import dataset_path, init_batch, output_path


@click.command()
@click.option('--input-mt-path', required=True)  # 'mt/v7.mt'
@click.option('--annotation-df=path', required=True)  # 'tob_wgs_rv/open_chromatin_annotation/predicted_l1_celltypes_avg_peaks_chr21.csv'
@click.option('--output-mt-path', required=True)  # 'tob_wgs_rv/open_chromatin_annotation/v7_open_chromatin.mt'
def annotate_variants(
    input_mt_path: str,
    annotation_df_path: str,
    output_mt_prefix: str,
):
    init_batch()

    # read hail matrix table object (WGS data)
    mt = hl.read_matrix_table(dataset_path(input_mt_path))
    mt = hl.experimental.densify(mt)

    # read in annotation data frame
    openchr_df = pd.read_csv(dataset_path(annotation_df_path), index_col=0)
    # import as hail Table
    openchr_ht = hl.Table.from_pandas(openchr_df) 
    # annotation
    annotated_mt = mt.annotate_rows(peak_avg_count = openchr_ht)  
    # save mt
    annotated_mt.write(output_mt_path, overwrite=True)


if __name__ == '__main__':
    annotate_variants()
