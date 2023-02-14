"""
Righto then, what do?

- we have a directory of BED files, each with the external sample ID as a name
- we want to find the internal IDs for each sample

"""


import os
import click

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_workflows.batch import get_batch

from sample_metadata.model.analysis_type import AnalysisType
from sample_metadata.apis import AnalysisApi, SampleApi

REFERENCE = 'ffrom config?'


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
    """

    # find all BED files in the input directory
    bed_files = to_path(beds).glob('*.bed')

    # for each of the BED files, get the ext ID
    ext_ids = {file.name.split('.', maxsplit=1)[0]: str(file) for file in bed_files}

    # get lookup of all CPG: Ext IDs for these samples
    sample_map = get_samples_from_metamist(list(ext_ids.keys()))

    # get all (latest) CRAM files for these samples
    cram_map = get_cram_files(list(sample_map.values()))

    # set an output path to write files to
    release_cram = 'gs://cpg-tob-wgs-release/cram'

    # set the image and reference to use
    samtools_image = get_config()['images']['samtools']
    cram_reference = get_config()['references']['broad']['ref_fasta']

    # let's start up a hail batch!
    batch = get_batch('Generate CRAM subsets - tob-wgs')

    # read reference in once per batch
    batch_reference = batch.read_input(cram_reference)

    # iterate over  all samples & BED files
    for ext_id, bed_file in ext_ids.items():
        cpg_id = sample_map[ext_id]
        cram_file = cram_map[cpg_id]
        cram_job = batch.new_job(name=f'subset CRAM {ext_id}/{cpg_id}')
        cram_job.image(samtools_image)

        # big storage, probably not huge on compute (small regions)
        cram_job.storage('60G')

        # read in the cram data for this sample
        cram = batch.read_input_group(
            **{'cram': cram_file, 'cram.crai': f'{cram_file}.crai'}
        ).cram

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
        cram_job.command(
            f'samtools view '
            f'-T {batch_reference} '
            f'-L {sample_bed} '
            f'-C --write-index '
            f'-o {cram_job.output_cram["cram"]} {cram}'
        )

        cram_out_path = os.path.join(release_cram, ext_id)

        # write the CRAM and relevant index
        get_batch().write_output(cram_job.output_cram, cram_out_path)


if __name__ == '__main__':
    main()
