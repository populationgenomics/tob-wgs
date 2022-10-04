#!/usr/bin/env python3

import click
import logging
import hail as hl
import pandas as pd
from cpg_utils.hail_batch import dataset_path, init_batch, output_path, reference_path
from cloudpathlib import AnyPath

# use logging to print statements, display at info level
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

MT = dataset_path('mt/v7.mt')


# get gene name as argument using click
@click.command()
@click.option('--gene-name', required=True)
# @click.option('--sample-mapping-file', required=True)  # this is not needed anymore for a specific individuals, but we should make sure is filtered for individuals we have scRNA-seq data for
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
    get all variants within a given window,

    then filter for variants that:
    - are biallelic SNVs
    - have regulatory consequences based on VEP annotations

    also annotate variants with the following:
    - CADD score
    - MAF within the OneK1K sample
    - MAF in gnomad

    create the following table:
    gene ID | variant ID | position | CADD | MAF (OneK1K) | MAF (gnomad) | ..

    and convert mt into plink files for testing (using CellRegMap-RV)
    """

    # define output filename and check if it already exists
    output_filename = AnyPath(output_path(f'{gene_name}_{window_size}.csv'))
    logging.info(f'Output file: {output_filename}')

    # skip if file already exists
    if output_filename.exists():
        raise Exception(f'File {output_filename} already exists, exiting')

    init_batch()
    # get VEP-annotated WGS object, *swap path* (hail matrix table)
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
    # logging.info(f'Number of variants after filtering for SNVs: {mt.count()[0]}')

    # SKIP FOR NOW NO VEP ANNOTATIONS
    # # filter to only variants with some regulatory consequences
    # mt = mt.filter_rows(
    #     hl.len(mt.vep.regulatory_feature_consequences) > 0
    # )
    # logging.info(
    #     f'Number of variants after filtering for variants with regulatory consequences: {mt.count()[0]}'
    # )

    logging.info('Export genotypes to PLINK')
    plink_out_path = AnyPath(
        output_path(f'{gene_name}_{window_size}')
    )  # figure out syntax
    hl.export_plink(mt, plink_out_path, ind_id=mt.s)

    # annotate variants with CADD scores, gnomad etc
    logging.info('Annotate variants with CADD scores and gnomad AF')
    # ref_ht = hl.read_table(
    #     'gs://cpg-reference/seqr/v0-1/combined_reference_data_grch38-2.0.4.ht'
    # )
    ref_ht = reference_path('seqr/v0-1/combined_reference_data_grch38-2.0.4.ht')
    mt = mt.annotate_rows(
        cadd=ref_ht[mt.row_key].cadd,
        gnomad_genomes=ref_ht[mt.row_key].gnomad_genomes,
    )

    # get CADD scores
    cadd_list = mt.cadd.PHRED.collect()

    # get gnomad MAF (popmax now, need NFE too?)
    maf_popmax = mt.gnomad_genomes.AF_POPMAX_OR_GLOBAL.collect()

    logging.info('Get AF within OneK1K sample')
    # get MAF within sample by calculating their variant QC
    mt = hl.variant_qc(mt)
    maf_sample_list = mt.variant_qc.AF[1].collect()

    # # get relevant regulatory info from VEP (biotype)
    # regulatory_consequences = mt.vep.regulatory_feature_consequences[
    #     'biotype'
    # ].collect()

    # get gene info
    gene_start = gene_df[gene_df['gene_name'] == gene_name]['start']
    gene_end = gene_df[gene_df['gene_name'] == gene_name]['end']
    ensembl_gene_id = gene_df[gene_df['gene_name'] == gene_name].index[0]

    logging.info('Preparing results data')
    results_data = {
        'gene_name': gene_name,
        'variant_id': mt.row_key[0].collect(),
        'cadd': cadd_list,
        'maf_onek1k': maf_sample_list,
        'maf_gnomad_popmax': maf_popmax,
        # 'regulatory_consequences': regulatory_consequences,
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
