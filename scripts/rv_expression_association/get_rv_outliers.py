import click
import hail as hl
import pandas as pd

# get OneK1K sample ID as an argument using click
@click.command()
@click.option('--onek1k-id', required=True)
@click.option('--gene-name', required=True)
@click.option('--chrom', required=True)

# for a given individual,
# and gene(s?) for which that individual is an expression outlier
# get all variants within a window, and annotate with the following:
# CADD score
# MAF within OneK1K
# MAF in gnomad
# regulatory information from VEP
# other?

# create the following table:
# donor ID | gene ID | variant ID | position | CADD | MAF (OneK1K) | MAF (gnomad) | ..


def main(
    onek1k_id: str,
    gene_name: str,
    chrom: str,
):

    # 'CPG9951' (943_944) is an expression outlier for gene: 'IGLL5'

    # file matching OneK1K IDs to CPG (internal) and TOB (external) IDs
    sample_key_filename = (
        'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/OneK1K_CPG_IDs.tsv'
    )
    sample_key_df = pd.read_csv(sample_key_filename, sep='\t', index_col=0)

    cpg_id = sample_key_df[sample_key_df['OneK1K_ID'] == onek1k_id]['InternalID']
    # switch print to log
    print(cpg_id)  # 'CPG9951'

    # define output filename and check if it already exists
    output_filename = cpg_id + '-' + gene_name + '.csv'

    # get VEP-annotated WGS object (hail matrix table)
    mt = hl.read_matrix_table('gs://cpg-tob-wgs-test/tob_wgs_vep/v1/vep105_GRCh38.mt')

    # select matrix down to that one donor
    donor_mt = mt.filter_cols(mt.s == cpg_id)

    # check this file in the future (GENCODE??)
    # get gene body position (start and end) and build interval
    # include variants up to 10kb up- and downstream
    gene_file = (
        'gs://cpg-tob-wgs-main/scrna-seq/grch38_association_files/gene_location_files/GRCh38_geneloc_chr'
        + chrom
        + '.tsv'
    )
    gene_df = pd.read_csv(gene_file, sep='\t', index_col=0)
    interval_start = float(gene_df[gene_df['gene_name'] == gene_name]['start']) - 10000
    interval_end = float(gene_df[gene_df['gene_name'] == gene_name]['end']) + 10000

    # get gene-specific genomic interval
    gene_interval = (
        'chr' + chrom + ':' + str(interval_start) + '-' + str(interval_end)
    )  # 'chr22:23219960-23348287'
    print(gene_interval)  # switch print to log
    donor_mt = hl.filter_intervals(
        donor_mt, [hl.parse_locus_interval(gene_interval, reference_genome='GRCh38')]
    )

    # remove variants for which this individual is 0/0
    donor_mt = hl.variant_qc(donor_mt)
    donor_mt = donor_mt.filter_rows(donor_mt.variant_qc.n_non_ref > 0)

    # focus on SNVs for now
    donor_mt = donor_mt.filter_rows(donor_mt.vep.variant_class == 'SNV')
    # filter for biallelic only
    donor_mt = donor_mt.filter_rows(hl.len(donor_mt.alleles) == 2)

    # annotate variants with CADD scores, gnomad etc
    ref_ht = hl.read_table(
        'gs://cpg-reference/seqr/v0-1/combined_reference_data_grch38-2.0.4.ht'
    )
    donor_mt = donor_mt.annotate_rows(
        cadd=ref_ht[donor_mt.row_key].cadd,
        gnomad_genomes=ref_ht[donor_mt.row_key].gnomad_genomes,
    )

    # get CADD scores
    cadd_list = donor_mt.cadd.PHRED.collect()

    # get gnomad MAF (both popmax and nfe)

    # get MAF within sample 
    # by selecting the appropriate rows in the full object
    relevant_loci = donor_mt.row_key.collect()
    mt = mt.filter_rows(hl.set(relevant_loci).contains(mt.row_key))
    # calculating their variant QC
    mt = hl.variant_qc(mt)
    maf_sample_list = mt.AF[0].collect()

    # get relevant regulatory info from VEP (enough?)

    results_data = {
        'onek1k_id': onek1k_id,
        'cpg_id': cpg_id,
        'gene_name': gene_name,
        'variant_id': donor_mt.row_key[0].collect(),
        'cadd': cadd_list,
        'maf (OneK1K)': maf_sample_list,
    }
    df = pd.DataFrame(results_data)
    df.to_csv(output_filename)


if __name__ == '__main__':
    main()
