import logging
import subprocess
from os.path import basename

import click

import hail as hl

from cpg_utils import to_path
from cpg_utils.config import get_config
from cpg_utils.hail_batch import init_batch
from metamist.graphql import gql, query

ACTIVE_INACTIVE_QUERY = gql(
    """
    query getActiveInactiveSGs($project: String!) {
        project(name: $project) {
            id
            samples {
                externalId
                id
                inactive_sgs: sequencingGroups(activeOnly: {eq: false}) {
                    id
                }
                active_sgs: sequencingGroups(activeOnly: {eq: true}) {
                    id
                }
            }
        }
    }
    """,
)

config = get_config()


def get_dict_of_gvcf_directories(project: str, nagim: bool = False) -> dict[str, str]:
    sgid_gvcf_path_map = {}
    prefix = f'gs://cpg-{project}/gvcf/{"nagim" if nagim else ""}'

    # Get GCP file paths
    cmd = ['gsutil', 'ls', prefix]
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        shell=True,  # noqa: S602
        check=False,
    )
    if result.returncode == 0:
        gvcf_paths: list[str] = result.stdout.splitlines()
        for gvcf_path in gvcf_paths:
            if gvcf_path.endswith('.gz'):
                sgid = basename(gvcf_path).split('.')[0]
                sgid_gvcf_path_map[sgid] = gvcf_path
    return sgid_gvcf_path_map


def create_combiner(
    project: str,
    gvcf_paths: list[str],
    sample_names: list[str],
    external_header: str,
    nagim: bool = False,
) -> hl.vds.combiner.VariantDatasetCombiner:
    out_path = (
        f'gs://cpg-{project}/dragmap_parity/{"nagim_" if nagim else "new_"}vds.vds'
    )
    # make Combiner objects
    combiner = hl.vds.new_combiner(
        output_path=out_path,
        temp_path=f'gs://cpg-{project}-tmp/dragmap_parity/',
        gvcf_paths=gvcf_paths,
        gvcf_sample_names=sample_names,
        gvcf_external_header=external_header,
        reference_genome='GRCh38',
        use_genome_default_intervals=True,
    )
    return combiner, out_path


def create_keyed_hail_tables(active_inactive_sg_map: dict[str, dict[str, str]]):
    tob_ids = []
    active_sgids = []
    inactive_sgids = []
    for tobid, sgids in active_inactive_sg_map.items():
        tob_ids.append(tobid)
        active_sgids.append(sgids['active'])
        inactive_sgids.append(sgids['inactive'])

    data = [
        {'tobid': tob_ids[i], 'active': active_sgids[i], 'inactive': inactive_sgids[i]}
        for i in range(len(tob_ids))
    ]
    ht_active_key = hl.Table.parallelize(
        data,
        hl.tstruct(tobid=hl.tstr, active=hl.tstr, inactive=hl.tstr),
        key='active',
    )
    ht_inactive_key = hl.Table.parallelize(
        data,
        hl.tstruct(tobid=hl.tstr, active=hl.tstr, inactive=hl.tstr),
        key='inactive',
    )
    return ht_active_key, ht_inactive_key


def check_active_inactive_gvcf_paths(
    project: str,
    active_inactive_sg_map: dict[str, dict[str, str]],
):
    nagim_gvcf_paths_dict: dict[str, str] = get_dict_of_gvcf_directories(
        project,
        nagim=True,
    )
    gvcf_paths_dict: dict[str, str] = get_dict_of_gvcf_directories(project)

    nagim_sgids = set(nagim_gvcf_paths_dict.keys())
    new_sgids = set(gvcf_paths_dict.keys())

    checked_nagim_gvcf_paths = []
    checked_new_gvcf_paths = []
    expids = []
    for expid, sgids in active_inactive_sg_map.items():
        inactive_sgids = sgids['inactive']
        active_sgid = sgids['active']
        # check inactive sgid has nagim gvcf
        if inactive_sgids in nagim_sgids and active_sgid in new_sgids:
            checked_nagim_gvcf_paths.append(nagim_gvcf_paths_dict[inactive_sgids])
            checked_new_gvcf_paths.append(gvcf_paths_dict[active_sgid])
            expids.append(expid)

    return checked_nagim_gvcf_paths, checked_new_gvcf_paths, expids


def get_active_inactive_sg_map(
    project: str,
    samples_to_skip: list[str],
) -> dict[str, dict[str, str]]:
    if '-main' in project:
        project = project.replace('-main', '')

    active_inactive_response = query(
        ACTIVE_INACTIVE_QUERY,
        variables={'project': project},
    )

    active_inactive_sg_map = {}
    for sample in active_inactive_response['project']['samples']:
        if len(sample['inactive_sgs']) >= 1:
            if sample['externalId'] in samples_to_skip:
                logging.info(f'Skipping {sample["externalId"]}')
                continue
            assert len(sample['active_sgs']) == len(sample['inactive_sgs']) == 1
            active_inactive_sg_map[sample['externalId']] = {
                'active': sample['active_sgs'][0]['id'],
                'inactive': sample['inactive_sgs'][0]['id'],
            }
    return active_inactive_sg_map


def create_vds(
    project: str,
    checked_nagim_gvcf_paths: list[str],
    checked_new_gvcf_paths: list[str],
    expids: list[str],
) -> hl.vds.VariantDataset:
    # make Combiner objects
    nagim_combiner, nagim_vds_path = create_combiner(
        project=project,
        gvcf_paths=checked_nagim_gvcf_paths,
        sample_names=expids,
        external_header=checked_nagim_gvcf_paths[0],
        nagim=True,
    )
    new_combiner, new_vds_path = create_combiner(
        project=project,
        gvcf_paths=checked_new_gvcf_paths,
        sample_names=expids,
        external_header=checked_new_gvcf_paths[0],
    )

    # run the combiners
    nagim_combiner.run()
    new_combiner.run()

    return hl.vds.read_vds(nagim_vds_path), hl.vds.read_vds(new_vds_path)


def rekey_matrix_table(mt: hl.MatrixTable, keyed_ref_table: hl.Table) -> hl.MatrixTable:
    """
    `keyed_ref_table` must have a column `tobid`
    """

    # Annotate the MatrixTable with the tobid from keyed_ref_table
    mt = mt.annotate_cols(tobid=keyed_ref_table[mt.s].tobid)

    # Capture the columns where tobid is None
    none_tobid_samples = mt.filter_cols(hl.is_missing(mt.tobid)).s.collect()

    # Print or log the columns that will be set to None (optional)
    logging.info(
        f"""Columns that will be set to None and dropped because of could not find corresponding 'newer' CPGID:" \
        {none_tobid_samples}""",
    )

    # Filter out columns where tobid is None
    mt = mt.filter_cols(hl.is_defined(mt.tobid))

    # Re-key the MatrixTable by the tobid
    return mt.key_cols_by('tobid')


@click.command()
@click.option('--project', required=True)
@click.option('--output-version', required=True)
@click.option('--new-vds-path', required=False, default=None)
@click.option('--nagim-vds-path', required=False, default=None)
@click.option('--nagim-mt-path', required=False, default=None)
@click.option('--samples-to-skip', required=False, multiple=True, default=None)
@click.option('--test/--no-test', default=True)
def main(
    project: str,
    output_version: str,
    new_vds_path: str | None,
    nagim_vds_path: str | None,
    nagim_mt_path: str | None,
    samples_to_skip: tuple[str] | None = None,
    test: bool = True,
):
    # 2 samples were duplicates and have no new CPG IDs
    # https://centrepopgen.slack.com/archives/C018KFBCR1C/p1709006419612279?thread_ts=1708911201.676049&cid=C018KFBCR1C
    if samples_to_skip:
        print(f"Skipping samples: {', '.join(samples_to_skip)}")
    else:
        print("No samples to skip.")

    project = project + '-test' if test else project + '-main'

    logging.info(f'Running dragmap_parity_check for project {project}')

    active_inactive_sg_map = get_active_inactive_sg_map(project, list(samples_to_skip))

    ht_active_key, ht_inactive_key = create_keyed_hail_tables(active_inactive_sg_map)

    new_vds = hl.vds.read_vds(new_vds_path)

    # read in vds/matrixtables
    if nagim_vds_path:
        nagim_vds = hl.vds.read_vds(nagim_vds_path)
    elif nagim_mt_path:
        nagim_mt = hl.read_matrix_table(nagim_mt_path)
    else:
        # create vds' from gvcf paths
        # TODO: Test the below still works (need to delete currently saved vds')
        (
            checked_nagim_gvcf_paths,
            checked_new_gvcf_paths,
            expids,
        ) = check_active_inactive_gvcf_paths(
            project,
            active_inactive_sg_map,
        )
        nagim_vds, new_vds = create_vds(
            project,
            checked_nagim_gvcf_paths,
            checked_new_gvcf_paths,
            expids,
        )

    logging.info(f'nagim path: {nagim_vds_path if nagim_vds_path else nagim_mt_path}')
    logging.info(f'new path: {new_vds_path}')
    # prepare vds' for comparison
    # As per documentation, hl.methods.concordance() requires the dataset to contain no multiallelic variants.
    # as well as the entry field to be 'GT', also expects MatrixTable, not a VariantDataset
    if nagim_vds_path:
        nagim_vds = hl.vds.split_multi(nagim_vds, filter_changed_loci=True)
        nagim_mt = nagim_vds.variant_data
    else:
        nagim_mt = hl.split_multi_hts(nagim_mt)

    new_vds = hl.vds.split_multi(new_vds, filter_changed_loci=True)
    new_mt = new_vds.variant_data

    # checkpoint the split MatrixTables
    output_prefix = f'gs://cpg-{project}-analysis/dragmap_parity/{output_version}'
    logging.info(f'Output prefix {output_prefix}')

    # rekey the MatrixTables
    nagim_mt = rekey_matrix_table(nagim_mt, ht_inactive_key)
    new_mt = rekey_matrix_table(new_mt, ht_active_key)

    # compare the two VDS'
    summary, samples, variants = hl.concordance(nagim_mt, new_mt)

    # write the concordance results
    logging.info(
        f'Writing concordance results to {output_prefix}',
    )
    samples.checkpoint(f'{output_prefix}/cols_concordance.ht')
    variants.checkpoint(f'{output_prefix}/variants_concordance.ht')
    with to_path(f'{output_prefix}/summary_concordance.txt').open('w') as f:
        f.write(str(summary))


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    init_batch()
    main()
