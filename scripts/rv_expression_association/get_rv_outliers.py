import hail as hl

mt = hl.read_matrix_table('gs://cpg-tob-wgs-test/tob_wgs_vep/v1/vep105_GRCh38.mt')

# 'CPG9951' (943_944) is an expression outlier for gene: 'IGLL5'
# select matrix down to that one donor 
donor_mt = mt.filter_cols(mt.s == 'CPG9951')

# from this file on Garvan HPC: /share/ScratchGeneral/anncuo/OneK1K/GeneLocations.tsv
# the row corresposnding to IGLL5 looks like this:
# gene_name  gene_id          seqid  start     end       strand
# IGLL5      ENSG00000254709  22     23229960  23238287  +
# so adding 10kb upstream and downstream, 
# i get the interval: 23,229,960 - 10,000 : 23,238,287 + 10,000: 22:23219960-23348287
donor_mt = hl.filter_intervals(donor_mt, [hl.parse_locus_interval('chr22:23219960-23348287', reference_genome='GRCh38')])

# remove variants for which this individual is 0/0
donor_mt = hl.variant_qc(donor_mt)
donor_mt = donor_mt.filter_rows(donor_mt.variant_qc.n_non_ref > 0)

# annotate variants with CADD scores etc
ref_ht = hl.read_table('gs://cpg-reference/seqr/v0-1/combined_reference_data_grch38-2.0.4.ht')
donor_mt = donor_mt.annotate_rows(cadd=ref_ht[donor_mt.row_key].cadd)

# get CADD scores
cadd_list = donor_mt.cadd.PHRED.collect()

# plot histogram of CADD scores
dp_hist = donor_mt.aggregate_entries(hl.expr.aggregators.hist(donor_mt.cadd.PHRED, 0, 30, 30))
p = hl.plot.histogram(dp_hist, legend='CADD score', title='CADD score Histogram')

# figure out how / and if I can save that list
# how to save the plot