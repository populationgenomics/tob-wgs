import hail as hl
import pandas as pd
from cloudpathlib import AnyPath
from cpg_utils.hail_batch import (
    # output_path,
    init_batch,
)


def filter_joint_call_mt(
    *,
    keys_path: str,
    joint_mt_path: str,
    frequency_table_path: str,
    vep_annotation_path: str,
    output_path: str,
    force: bool = False
):
    """Filter hail matrix table

    Input:
    keys_path: path to a tsv file with information on
    OneK1K amd CPG IDs (columns) for each sample (rows).

    Returns:
    Path to a hail matrix table, with rows (alleles) filtered on the following requirements:
    1) biallelic, 2) meets VQSR filters, 3) gene quality score higher than 20,
    4) call rate of 0.8, and 5) variants with MAF <= 0.01.
    """

    init_batch()
    if AnyPath(output_path).exists() and not force:
        return output_path

    mt = hl.read_matrix_table(joint_mt_path)
    samples = mt.s.collect()
    n_samples = len(samples)
    mt = mt.naive_coalesce(10000)
    mt = hl.experimental.densify(mt)
    # filter to biallelic loci only
    mt = mt.filter_rows(hl.len(mt.alleles) == 2)
    # filter out variants that didn't pass the VQSR filter
    mt = mt.filter_rows(hl.len(hl.or_else(mt.filters, hl.empty_set(hl.tstr))) == 0)
    # VQSR does not filter out low quality genotypes. Filter these out
    mt = mt.filter_entries(mt.GQ <= 20, keep=False)
    # filter out samples with a genotype call rate > 0.8 (as in the gnomAD supplementary paper)
    call_rate = 0.8
    mt = mt.filter_rows(
        hl.agg.sum(hl.is_missing(mt.GT)) > (n_samples * call_rate), keep=False
    )
    # filter out variants with MAF <= 0.01
    ht = hl.read_table(frequency_table_path)
    mt = mt.annotate_rows(freq=ht[mt.row_key].freq)
    mt = mt.filter_rows(mt.freq.AF[1] > 0.01)
    # add in VEP annotation
    vep = hl.read_table(vep_annotation_path)
    mt = mt.annotate_rows(
        vep_functional_anno=vep[mt.row_key].vep.regulatory_feature_consequences.biotype
    )
    # add OneK1K IDs to genotype mt
    sampleid_keys = pd.read_csv(AnyPath(keys_path), sep='\t')
    genotype_samples = pd.DataFrame(samples, columns=['sampleid'])
    sampleid_keys = pd.merge(
        genotype_samples,
        sampleid_keys,
        how='left',
        left_on='sampleid',
        right_on='InternalID',
    )
    sampleid_keys.fillna('', inplace=True)
    sampleid_keys = hl.Table.from_pandas(sampleid_keys)
    sampleid_keys = sampleid_keys.key_by('sampleid')
    mt = mt.annotate_cols(onek1k_id=sampleid_keys[mt.s].OneK1K_ID)
    # repartition to save overhead cost
    mt = mt.naive_coalesce(1000)
    mt.write(output_path)
    return output_path
