#!/usr/bin/env python3

import click
import logging
import hail as hl
import pandas as pd
from cpg_utils.hail_batch import reference_path
from cpg_utils.hail_batch import output_path, init_batch
from cloudpathlib import AnyPath

# use logging to print statements, display at info level
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)


# get OneK1K sample ID and gene name as arguments using click
@click.command()
@click.option('--onek1k-id', required=True)
@click.option('--gene-name', required=True)
@click.option('--chrom', required=True)
@click.option(
    '--window-size',
    required=True,
    help='the size of the flanking region, to be added twice, on either side of the gene body',
)
def main(
    onek1k_id: str,
    gene_name: str,
    chrom: str,
    window_size: int,
):
    """for a given individual,
    and gene for which that individual is an expression outlier
    get all variants within a given window,

    then filter for variants that:
    - are alt (0/1 or 1/1) for that individual
    - are biallelic SNVs
    - have regulatory consequences based on VEP annotations

    also annotate variants with the following:
    - CADD score
    - MAF within the OneK1K sample
    - MAF in gnomad

    create the following table:
    donor ID | gene ID | variant ID | position | CADD | MAF (OneK1K) | MAF (gnomad) | ..
    """

    # file matching OneK1K IDs to CPG (internal) and TOB (external) IDs
    sample_key_filename = (
        'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/OneK1K_CPG_IDs.tsv'
    )
    sample_key_df = pd.read_csv(AnyPath(sample_key_filename), sep='\t', index_col=0)

    cpg_id = sample_key_df[sample_key_df.index == onek1k_id]['InternalID'].values[0]
    logging.info(f'CPG ID: {cpg_id}')  # e.g., 'CPG9951'

    # define output filename and check if it already exists
    output_filename = AnyPath(output_path(f'{cpg_id}_{gene_name}_{window_size}.csv'))
    logging.info(f'Output file: {output_filename}')

    # skip if file already exists
    if output_filename.exists():
        raise Exception(f'File {output_filename} already exists, exiting')

    init_batch()
    # get VEP-annotated WGS object (hail matrix table)
    mt = hl.read_matrix_table('gs://cpg-tob-wgs-test/tob_wgs_vep/v1/vep105_GRCh38.mt')
    logging.info(f'Number of total variants: {mt.count()[0]}')

    # filter out low QC variants
    mt = mt.filter_rows(hl.len(mt.filters) == 0)

    # select matrix down to that one donor
    donor_mt = mt.filter_cols(mt.s == cpg_id)

    logging.info(f'Number of QC-passing variants: {donor_mt.count()[0]}')

    # check this file in the future (GENCODE??)
    # get gene body position (start and end) and build interval
    # include variants up to 10kb up- and downstream
    gene_file = (
        'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/gene_location_files/GRCh38_geneloc_chr'
        + chrom
        + '.tsv'
    )
    gene_df = pd.read_csv(AnyPath(gene_file), sep='\t', index_col=0)
    interval_start = int(gene_df[gene_df['gene_name'] == gene_name]['start']) - int(
        window_size
    )
    interval_end = int(gene_df[gene_df['gene_name'] == gene_name]['end']) + int(
        window_size
    )

    # get gene-specific genomic interval
    gene_interval = 'chr' + chrom + ':' + str(interval_start) + '-' + str(interval_end)
    logging.info(
        f'Interval considered: {gene_interval}'
    )  # 'chr22:23219960-23348287'

    donor_mt = hl.filter_intervals(
        donor_mt, [hl.parse_locus_interval(gene_interval, reference_genome='GRCh38')]
    )
    logging.info(f'Number of variants within interval: {donor_mt.count()[0]}')

    # remove variants for which this individual is 0/0
    donor_mt = hl.variant_qc(donor_mt)
    donor_mt = donor_mt.filter_rows(donor_mt.variant_qc.n_non_ref > 0)
    logging.info(
        f'Number of non-ref variants for this indvidual: {donor_mt.count()[0]}'
    )

    # focus on SNVs for now
    donor_mt = donor_mt.filter_rows(donor_mt.vep.variant_class == 'SNV')
    # filter for biallelic only
    donor_mt = donor_mt.filter_rows(hl.len(donor_mt.alleles) == 2)
    logging.info(
        f'Number of variants after filtering for biallelic SNVs: {donor_mt.count()[0]}'
    )

    # filter to only variants with some regulatory consequences
    donor_mt = donor_mt.filter_rows(
        hl.len(donor_mt.vep.regulatory_feature_consequences) > 0
    )
    logging.info(
        f'Number of variants after filtering for variants with regulatory consequences: {donor_mt.count()[0]}'
        )
    )

    # retain info about het variants
    is_het_list = donor_mt.variant_qc.n_het.collect()

    # annotate variants with CADD scores, gnomad etc
    logging.info('Annotate variants with CADD scores and gnomad AF')
    # ref_ht = hl.read_table(
    #     'gs://cpg-reference/seqr/v0-1/combined_reference_data_grch38-2.0.4.ht'
    # )
    ref_ht = reference_path('seqr/v0-1/combined_reference_data_grch38-2.0.4.ht')
    donor_mt = donor_mt.annotate_rows(
        cadd=ref_ht[donor_mt.row_key].cadd,
        gnomad_genomes=ref_ht[donor_mt.row_key].gnomad_genomes,
    )

    # get CADD scores
    cadd_list = donor_mt.cadd.PHRED.collect()

    # get gnomad MAF (popmax now, need NFE too?)
    maf_popmax = donor_mt.gnomad_genomes.AF_POPMAX_OR_GLOBAL.collect()

    logging.info('Get AF within OneK1K sample')
    # get MAF within sample
    # by selecting the appropriate rows in the full object
    relevant_loci = donor_mt.row_key.collect()
    mt = mt.filter_rows(hl.set(relevant_loci).contains(mt.row_key))
    # calculating their variant QC
    mt = hl.variant_qc(mt)
    maf_sample_list = mt.variant_qc.AF[1].collect()

    # get relevant regulatory info from VEP (biotype)
    regulatory_consequences = donor_mt.vep.regulatory_feature_consequences[
        'biotype'
    ].collect()

    # get gene info
    gene_start = gene_df[gene_df['gene_name'] == gene_name]['start']
    gene_end = gene_df[gene_df['gene_name'] == gene_name]['end']
    ensembl_gene_id = gene_df[gene_df['gene_name'] == gene_name].index[0]

    logging.info('Preparing results data')
    results_data = {
        'onek1k_id': onek1k_id,
        'cpg_id': cpg_id,
        'gene_name': gene_name,
        'variant_id': donor_mt.row_key[0].collect(),
        'is_het': is_het_list,
        'cadd': cadd_list,
        'maf_onek1k': maf_sample_list,
        'maf_gnomad_popmax': maf_popmax,
        'regulatory_consequences': regulatory_consequences,
        'ensembl_gene_id': ensembl_gene_id,
        'chrom': chrom,
        'gene_start': int(gene_start),
        'gene_end': int(gene_end),
    }
    df = pd.DataFrame(results_data)
    with output_filename.open('w') as of:
        df.to_csv(of, index=False)


if __name__ == '__main__':
    main()
