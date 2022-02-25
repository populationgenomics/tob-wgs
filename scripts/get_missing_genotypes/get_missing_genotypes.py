"""Get GQ and GT for missing variants"""

import hail as hl

NAGIM = 'gs://cpg-nagim-test/mt/v1-3.mt/'


def get_gt_and_gq(mt, locus, sample):
    """Get GQ and GT scores for individual samples"""
    na_locus = hl.parse_locus(locus, 'GRCh38')
    parsed_locus = mt.filter_rows(mt.locus == na_locus, keep=True)
    gt = parsed_locus.filter_cols(parsed_locus.s == sample).GT.collect()[0].alleles
    gq = parsed_locus.filter_cols(parsed_locus.s == sample).GQ.collect()
    print(f'The GT for sample {sample} is {gt} and the GQ score is {gq}')


def query():
    """Query script entry point"""

    hl.init(default_reference='GRCh38')

    nagim = hl.read_matrix_table(NAGIM)
    get_gt_and_gq(mt=nagim, locus='chr1:28663', sample='CPG1107')
    get_gt_and_gq(mt=nagim, locus='chr1:94725', sample='CPG1198')
    get_gt_and_gq(mt=nagim, locus='chr1:588091', sample='CPG1081')


if __name__ == '__main__':
    query()
