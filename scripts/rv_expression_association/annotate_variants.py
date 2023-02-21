#!/usr/bin/env python3
# pylint: disable=missing-function-docstring,missing-module-docstring,no-value-for-parameter

import click
import hail as hl
import pandas as pd
from cpg_utils.hail_batch import dataset_path, init_batch, output_path


@click.command()
@click.option('--input-mt-path', required=True)  # 'mt/v7.mt'
@click.option(
    '--annotation-df-path', required=True
)  # 'tob_wgs_rv/open_chromatin_annotation/predicted_l1_celltypes_avg_peaks_chr21.csv'
@click.option(
    '--output-ht-path', required=True
)  # 'tob_wgs_rv/open_chromatin_annotation/open_chromatin_annotated.ht'
def annotate_variants(
    input_mt_path: str,
    annotation_df_path: str,
    output_ht_path: str,
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
    # this data frame has, for the genomic region, columns indicating average chromatin openness per cell type
    #                        |       B      |     CD4 T    | ...
    # chr21-5065291-5066183  |   3.643175   |   1.791078   | ...
    openchr_df = pd.read_csv(dataset_path(annotation_df_path), index_col=0)

    # turn into (there may be a better way):
    #         interval       |       B      |     CD4 T    | ...
    # chr21:5065291-5066183  |   3.643175   |   1.791078   | ...
    openchr_df['interval'] = [index.replace('-', ':', 1) for index in openchr_df.index]

    # import as hail Table
    openchr_ht = hl.Table.from_pandas(openchr_df)
    # parse as interval, set as key
    openchr_ht = openchr_ht.annotate(
        interval=hl.parse_locus_interval(openchr_ht.interval, reference_genome='GRCh38')
    ).key_by('interval')

    # annotate
    mt = mt.filter_rows(mt.locus.contig == 'chr21')  # object in test is chr21 only
    variants_ht = mt.rows()
    variants_ht = variants_ht.annotate(
        open_chromatin=openchr_ht.index(variants_ht.locus)
    )

    # save ht
    variants_ht.write(output_path(output_ht_path), overwrite=True)


if __name__ == '__main__':
    annotate_variants()
