#!/usr/bin/env python3

# from bokeh.io.export import get_screenshot_as_png

# from cpg_utils.config import get_config
from cpg_utils.hail_batch import dataset_path, init_batch  #, output_path

# , reference_path
import hail as hl

HT = dataset_path('tob_wgs_vep/104/vep104.3_GRCh38.ht')


def main():
    init_batch()

    ht = hl.read_table(VEP_HT)
    ht = ht.filter_rows(ht.locus.contig == 'chr22')  # check syntax
    print(ht.count())
    # ht = hl.experimental.densify(ht)
    # ht = hl.filter_intervals(
    #     ht,
    #     [hl.parse_locus_interval('chr22:23219960-23348287', reference_genome='GRCh38')],
    # )
    print(ht.count())
    # ht = hl.variant_qc(ht)
    # ht = ht.filter_rows(hl.len(hl.or_else(ht.filters, hl.empty_set(hl.tstr))) == 0)
    print(ht.count())
    # filter for biallelic
    ht = ht.filter_rows(hl.len(ht.alleles) == 2)
    print(ht.count())
    print(ht.freq.AF.show())


    # p1 = hl.plot.histogram(ht.variant_qc.AF[1])
    # p1_filename = output_path('histogram_alt_af_all_gene_variants.png', 'web')
    # with hl.hadoop_open(p1_filename, 'wb') as f:
    #     get_screenshot_as_png(p1).save(f, format='PNG')



if __name__ == '__main__':
    main()
