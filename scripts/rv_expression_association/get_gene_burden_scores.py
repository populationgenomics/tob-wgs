#!/usr/bin/env python3

import click
import logging
import hail as hl
import pandas as pd
from cpg_utils.hail_batch import dataset_path, init_batch, output_path
from cloudpathlib import AnyPath

# use logging to print statements, display at info level
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

MT = dataset_path('mt/v7.mt')


# get gene name as argument using click
@click.command()
@click.option('--gene-name', required=True)
@click.option('--gene-file', required=True)
@click.option(
    '--window-size',
    required=True,
    help='the size of the flanking region, to be added twice, on either side of the gene body',
)
def main(
    gene_name: str,
    gene_file: str,  # e.g., 'scrna-seq/grch38_association_files/gene_location_files/GRCh38_geneloc_chr1.tsv'
    window_size: int,
):
    """for a given gene
    get all (biallelic SNV) variants within a given window,

    obtain individual-level burden scores as
    - Number of low-frequency variants (1-5% MAF) within the interval
    - Number of rare variants (<1% MAF) within the interval

    create the following table:
    individual | gene score low freq | gene score rare variants 
    """

    # define output filenames and check if they already exist
    scores_filename = AnyPath(output_path(f'{gene_name}_{window_size}.csv'))
    logging.info(f'Output file: {scores_filename}')

    # skip if file already exists
    if scores_filename.exists():
        raise Exception(f'File {scores_filename} already exists, exiting')

    init_batch()
    # get WGS object, (hail matrix table)
    mt = hl.read_matrix_table(MT)
    logging.info(f'Number of total variants: {mt.count()[0]}')

    # get relevant chromosome
    gene_df = pd.read_csv(AnyPath(dataset_path(gene_file)), sep='\t', index_col=0)
    chrom = gene_df[gene_df['gene_name'] == gene_name]['chr']

    # filter to relevant chromosome to speed densification up
    mt = mt.filter_rows(mt.locus.contig == ('chr' + chrom))
    logging.info(f'Number of variants on chromosome {chrom}: {mt.count()[0]}')

    # densify
    mt = hl.experimental.densify(mt)

    # filter out low QC variants
    mt = mt.filter_rows(hl.len(hl.or_else(mt.filters, hl.empty_set(hl.tstr))) == 0)
    mt_path = output_path('densified_chrom_and_qc_filtered.mt', 'tmp')
    mt = mt.checkpoint(mt_path, overwrite=True)
    logging.info(f'Number of QC-passing variants: {mt.count()[0]}')

    # get gene body position (start and end) and build interval
    # include variants up to {window size} up- and downstream
    interval_start = int(gene_df[gene_df['gene_name'] == gene_name]['start']) - int(
        window_size
    )
    interval_end = int(gene_df[gene_df['gene_name'] == gene_name]['end']) + int(
        window_size
    )
    # clip to chromosome boundaries
    left_boundary = max(1, interval_start)
    right_boundary = min(
        interval_end, hl.get_reference('GRCh38').lengths[f'chr{chrom}']
    )
    # get gene-specific genomic interval
    gene_interval = f'chr{chrom}:{left_boundary}-{right_boundary}'
    logging.info(f'Interval considered: {gene_interval}')  # 'chr22:23219960-23348287'

    mt = hl.filter_intervals(
        mt, [hl.parse_locus_interval(gene_interval, reference_genome='GRCh38')]
    )
    logging.info(f'Number of variants within interval: {mt.count()[0]}')

    # focus on SNVs for now
    mt = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1]))
    # mt = mt.filter_rows(mt.vep.variant_class == 'SNV')

    # filter for biallelic only (this filters both multi-allelic and hom-ref)
    mt = mt.filter_rows(hl.len(mt.alleles) == 2)  # biallelic already
    logging.info(
        f'Number of variants after filtering for biallelic SNVs: {mt.count()[0]}'
    )

    logging.info('Get AF within OneK1K sample')
    # get MAF within sample by calculating their variant QC
    mt = hl.variant_qc(mt)
    alt_af_sample_list = mt.variant_qc.AF[1].collect()

    logging.info('Get within sample MAF rather than alt allele counts')
    maf_sample_list = [1-i if i >= 0.5 else i for i in alt_af_sample_list]

    low_frequency_variants = maf_sample_list[maf_sample_list>0.01 and maf_sample_list<0.05]
    rare_variants = maf_sample_list[maf_sample_list<0.01]

    logging.info('Count variants at different MAF windows')
    lf_scores: List[int] = [] # low frequency variants (1-5%)
    rv_scores: List[int] = [] # rare variants (<1%)

    samples = mt.s.collect()
    
    for sample in samples:
        donor_mt = mt.filter_cols(mt.s == sample)
        ## low frequency variants
        donor_mt_lq = donor_mt.filter_rows(hl.set(low_frequency_variants).contains(donor_mt.row_key))  # consider relevant variants
        donor_mt_lq = donor_mt_lq.filter_rows(hl.agg.any(donor_mt_lq.GT.is_non_ref())) # this will get non ref instead of MAF in sample
        lf_scores[sample] = donor_mt_lq.count()[0]
        ## rare variants
        donor_mt_rv = donor_mt.filter_rows(hl.set(rare_variants).contains(donor_mt.row_key))  # consider relevant variants
        donor_mt_rv = donor_mt_rv.filter_rows(hl.agg.any(donor_mt_rv.GT.is_non_ref())) # this will get non ref instead of MAF in sample
        rv_scores[sample] = donor_mt_rv.count()[0]
   

    logging.info('Preparing results data')
    results_data = {
        'individual': samples,
        'score_lowfreq': lf_scores,
        'score_rare': rv_scores,
    }
    df = pd.DataFrame(results_data)
    with scores_filename.open('w') as of:
        df.to_csv(of, index=False)


if __name__ == '__main__':
    main()
