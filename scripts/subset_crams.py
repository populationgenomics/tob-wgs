"""
Righto then, what do?

- we have a directory of BED files, each with the external sample ID as a name
- we want to find the internal IDs for each sample
- for each sample, create a job to create a mini-CRAM using its relevant BED
- write the resulting files to the main bucket

Uses GCS_OAUTH_TOKEN within the samtools container to read GCS CRAMs directly

Requires manual transfer later to the release bucket
"""


import click

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import authenticate_cloud_credentials_in_job
from cpg_workflows.batch import get_batch, dataset_path

from sample_metadata.model.analysis_type import AnalysisType
from sample_metadata.apis import AnalysisApi, SampleApi


def get_samples_from_metamist(samples: list[str]) -> dict[str, str]:
    """
    return a mapping of external to internal IDs
    """
    return SampleApi().get_sample_id_map_by_external(
        project='tob-wgs', request_body=samples
    )


def get_cram_files(samples: list[str]) -> dict[str, str]:
    """
    find the most recent registered CRAM file for each sample
    """
    results_list = AnalysisApi().get_latest_analysis_for_samples_and_type(
        analysis_type=AnalysisType('cram'), project='tob-wgs', request_body=samples
    )

    # return a mapping of CPG ID to CRAM file
    return {data['sample_ids'][0]: data['output'] for data in results_list}


@click.command
@click.option('--beds', help='directory to find BED files')
def main(beds: str):
    """
    get all relevant file paths for these samples
    produce a single job per sample, generating a subset CRAM
    write that cram + index to the release bucket
    (release bucket not present in config...)

    :param beds: str, path to cloud folder containing BED files
        BED file naming is "<EXT_ID>.bed"
    """

    # find all BED files in the input directory
    bed_files = to_path(beds).glob('*.bed')

    # from each of the BED files, get the ext ID
    ext_ids = {file.name.split('.', maxsplit=1)[0]: str(file) for file in bed_files}

    # get lookup of all CPG: Ext IDs for these samples
    sample_map = get_samples_from_metamist(list(ext_ids.keys()))

    # get all (latest) CRAM files for these samples
    cram_map = get_cram_files(list(sample_map.values()))

    # set the image and reference to use
    samtools_image = get_config()['images']['samtools']

    # let's start up a hail batch!
    batch = get_batch('Generate CRAM subsets - tob-wgs')

    # read reference in once per batch
    batch_reference = batch.read_input('gs://cpg-common-main/references/hg38/v0/Homo_sapiens_assembly38.fasta')

    # iterate over  all samples & BED files
    for ext_id, bed_file in ext_ids.items():
        cpg_id = sample_map[ext_id]
        cram_file = cram_map[cpg_id]
        cram_job = batch.new_job(name=f'subset CRAM {ext_id}/{cpg_id}')

        # authenticate credentials in job
        authenticate_cloud_credentials_in_job(cram_job)

        # run the job inside the samtools image
        cram_job.image(samtools_image)

        # big storage, probably not huge on compute (small regions)
        cram_job.storage('60G')
        cram_job.memory('16G')

        # read in the sample BED file
        sample_bed = batch.read_input(bed_file)

        # set a group for the output file (cram and index)
        cram_job.declare_resource_group(
            output_cram={'cram': '{root}.cram', 'cram.crai': '{root}.cram.crai'}
        )

        # samtools view
        # -T: CRAM reference fasta
        # -C: output CRAM
        # -L: target regions file, specific to sample
        # --write-index: ...
        # sets GCS_OAUTH_TOKEN to access data directly
        cram_job.command(
            'GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token) '
            'samtools view '
            f'-T {batch_reference} '
            f'-L {sample_bed} '
            '-C --write-index '
            f'-o {cram_job.output_cram["cram"]} '
            f'{cram_file}'
        )

        # write the CRAM and relevant index
        get_batch().write_output(
            cram_job.output_cram,
            dataset_path(f'MRR_cram_extracts/2023_01_30/{ext_id}.mini'),
        )

    get_batch().run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
