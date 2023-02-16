#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,missing-module-docstring,no-value-for-parameter

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
    output_mt_path: str,
):
    init_batch()

    # read hail matrix table object (WGS data)
    mt = hl.read_matrix_table(dataset_path(input_mt_path))
    mt = hl.experimental.densify(mt)
    
    # filter out low quality variants and consider only variable loci (no ref-only)
    mt = mt.filter_rows(
        (hl.len(hl.or_else(mt.filters, hl.empty_set(hl.tstr))) == 0)  # QC
        & (hl.len(mt.alleles) == 2)  # remove hom-ref
    )

    # read in annotation data frame
    # this data frame has columns indicating the genomic region and average openness in each cell type
    #                        |       B      |     CD4 T    | ...
    # chr21-5065291-5066183  |   3.643175   |   1.791078   | ...
    openchr_df = pd.read_csv(dataset_path(annotation_df_path), index_col=0)
    
    # turn into (there may be a better way):
    #         peak           |       B      |     CD4 T    | ...
    # chr21:5065291-5066183  |   3.643175   |   1.791078   | ...
    openchr_df["peak"] = str(openchr_df.index[0]).replace('-', ':', 1)
    
    # import as hail Table
    openchr_ht = hl.Table.from_pandas(openchr_df) 
    
    # annotation
    # not sure how to do this, does it make sense to loop over peaks / genomic region
    # and annotate all loci in the same way, then move on to the next peak?
    # sth like: loci[loci in interval].open_chromatin.{cell_type} = df[interval, {celltype}]
    annotated_mt = mt.annotate_rows(peak_avg_count=openchr_ht)  
    
    # save mt
    annotated_mt.write(output_path(output_mt_path), overwrite=True)


if __name__ == '__main__':
    annotate_variants()
