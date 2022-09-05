from cpg_utils.hail_batch import dataset_path

DEFAULT_JOINT_CALL_TABLE_PATH = dataset_path('mt/v7.mt/')
DEFAULT_FREQUENCY_TABLE_PATH = dataset_path(
    'joint-calling/v7/variant_qc/frequencies.ht/', 'analysis'
)
DEFAULT_VEP_ANNOTATION_TABLE_PATH = dataset_path('tob_wgs_vep/104/vep104.3_GRCh38.ht/')
DEFAULT_GENCODE_GTF_PATH = 'gs://cpg-reference/gencode/gencode.v38.annotation.gtf.bgz'  # reference_path('gencode_gtf'),

MULTIPY_IMAGE = 'australia-southeast1-docker.pkg.dev/cpg-common/images/multipy:0.16'
