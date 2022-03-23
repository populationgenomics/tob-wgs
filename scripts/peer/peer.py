"""
Test run this with:

    cd scripts/peer
    analysis-runner \
        --access-level test \
        --description 'Test run peer analysis' \
        --output-dir 'peer/2022-03-23' \
        --dataset tob-wgs \
        python3 peer.py \
        --scores-path 'gs://cpg-tob-wgs-test/tob_wgs_pca/nfe_no_outliers/v0/scores.json' \
        --covariates-path 'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/covariates_files/covariates.tsv' \
        --sample-id-keys-path 'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/metadata/keys_metadata_sheet.csv' \
        --path-to-cell-files 'gs://cpg-tob-wgs-test/scrna-seq/grch38_association_files/expression_files/'
"""


import os
import logging
import click
import hailtop.batch as hb
import pandas as pd
from typing import Tuple
from cpg_utils.hail import remote_tmpdir, output_path
from google.cloud import storage

PEER_DOCKER = 'australia-southeast1-docker.pkg.dev/cpg-common/images/peer:1.3.2'


def get_covariates(
    scores_path, covariates_path, expression_file, sample_id_keys_path
) -> Tuple[str, str]:
    """
    Get covariate data by merging PCA scores with age and sex info.
    Only needs to be run once. This returns a TSV (as a string)
    """
    # load in scores which have outliers removed
    scores_df = pd.read_json(scores_path)
    sampleid = scores_df.s
    # only get the first 4 PCs, as indicated by scree plot
    scores = pd.DataFrame(scores_df['scores'].to_list()).iloc[:, 0:4]
    scores.columns = ['PC1', 'PC2', 'PC3', 'PC4']
    scores.insert(loc=0, column='sampleid', value=sampleid)
    # filter to only CPG samples
    scores = scores[scores.sampleid.str.contains('CPG')]
    # change CPG IDs to One1K1K IDs
    sampleid_keys = pd.read_csv(sample_id_keys_path, sep=',')
    scores = pd.merge(
        scores, sampleid_keys, how='left', left_on='sampleid', right_on='CPG_ID'
    )
    # remove samples which don't have a OneK1K ID and only keep PCA scores
    scores = scores[~scores.OneK1K_ID.isna()][
        ['PC1', 'PC2', 'PC3', 'PC4', 'OneK1K_ID']
    ]
    # Merge PCA scores with Seyhan's covariate data and drop Seyhan's RNA-seq pca scores
    # NB: these are the ones in lowercase (e.g., 'pc1')
    covariates = pd.read_csv(covariates_path, sep='\t')
    covariates = pd.merge(
        scores, covariates, how='left', left_on='OneK1K_ID', right_on='sampleid'
    ).drop(['pc1', 'pc2', 'pc3', 'pc4'], axis=1)
    # check whether there's any missing data in the covariate data
    missing_cov_data = list(covariates[covariates.sampleid.isna()].OneK1K_ID)
    logging.info(f'{missing_cov_data} is missing covariate data')
    # Remove samples with no covariate data and drop
    # the OneK1K_ID, since it's redundant with sampleid
    covariates = covariates[~covariates.sampleid.isna()].drop(
        ['OneK1K_ID'], axis=1
    )
    # Match expression data to covariate data, since PEER needs expr and covs
    # dfs to be the same length
    expression = pd.read_csv(expression_file, sep='\t')
    merged_expr_covs = pd.merge(
        expression, covariates, how='right', left_on='sampleid', right_on='sampleid'
    )
    samples_missing_expression = merged_expr_covs.shape[0] - merged_expr_covs.dropna().shape[0]
    logging.info(f'{samples_missing_expression} samples are missing expression data')
    # remove any rows with NA values
    merged_expr_covs = merged_expr_covs.dropna()
    expression = merged_expr_covs.drop(
        ['PC1', 'PC2', 'PC3', 'PC4', 'sex', 'age'], axis=1
    )
    covariates = merged_expr_covs[
        ['sampleid', 'PC1', 'PC2', 'PC3', 'PC4', 'sex', 'age']
    ]
    # make sure samples in covariates and expression data are equal
    assert expression.sampleid.equals(covariates.sampleid)
    # remove sampleid from both covariates and expression dfs
    expression.drop('sampleid', axis=1, inplace=True)
    covariates.drop('sampleid', axis=1, inplace=True)
    # finally, make sure age and sex in covariate df are int
    covariates[['sex', 'age']] = covariates[['sex', 'age']].astype(int)

    # return expression data and covariates
    return covariates.to_csv(index=False), expression.to_csv(index=False)


def get_at_index(obj, idx):
    """
    Get the object at index, required to break up the result
    from the get_covariates data within Hail Batch    
    """

    return obj[idx]


def run_peer_job(job: hb.batch.job.Job, expression_file, covariates_file):
    """
    Run peer analysis, except because peer is in Python2, we can't take
    advantage of the python_job concept in Batch. Instead, we'll rely on
    the two inputs `expression_file` and `covariates_file` to be TSV files.
    """

    job.image(PEER_DOCKER)

    # write python script to container
    job.command(
        """
cat <<EOT >> run_peer.py

import peer
import numpy as np
import pandas as pd

def run_peer(expression_file, covariates_file, factors_output_path, weights_output_path, precision_output_path, residuals_output_path):
    \"""
    Get covariate data for each cell type
    \"""

    # load in data
    expr = pd.read_csv(expression_file, header=None, skiprows=1)
    covs = pd.read_csv(covariates_file, header=None, skiprows=1)

    # Set PEER paramaters as per the PEER website
    model = peer.PEER()
    model.setPhenoMean(expr)
    model.setCovariates(covs)
    model.setNk(10)
    model.update()

    # Calculate and save the PEER factors
    factors = model.getX()
    np.savetxt(factors_output_path, factors, delimiter=',')
    # Calculate and save the weights for each factor
    weights = model.getW()
    np.savetxt(weights_output_path, weights, delimiter=',')
    # Calculate and save the precision values
    precision = model.getAlpha()
    np.savetxt(precision_output_path, precision, delimiter=',')
    # Calculate and save the residuals
    residuals = model.getResiduals()
    np.savetxt(residuals_output_path, residuals, delimiter=',')


if __name__ == "__main__":
    import sys
    run_peer(*sys.argv[1:])

EOT"""
    )

    job.command(
        f"""python run_peer.py \
            {expression_file} \
            {covariates_file} \
            {job.factors_output_path} \
            {job.weights_output_path} \
            {job.precision_output_path} \
            {job.residuals_output_path}"""
    )

    return job


def process_cell_type_on_batch(
    batch: hb.Batch,
    cell_type_name: str,
    expression_file,
    scores_path,
    covariates_path,
    sample_id_keys_path,
):
    """
    Run PEER calculation in hail batch
    """
    expression_f = batch.read_input(expression_file)
    job_prefix = f'{cell_type_name}: '
    load_data = batch.new_python_job(job_prefix + 'load covariates')
    intermediate_tuple = load_data.call(
        get_covariates, scores_path, covariates_path, expression_f, sample_id_keys_path
    )
    covariates_csv = load_data.call(get_at_index, intermediate_tuple, 0).as_str()
    expression_csv = load_data.call(get_at_index, intermediate_tuple, 1).as_str()
    peer_job = batch.new_job(job_prefix + 'peer')
    run_peer_job(peer_job, expression_csv, covariates_csv)

    batch.write_output(
        peer_job.factors_output_path,
        output_path(f'{cell_type_name}_peer_factors_file.txt'),
    )
    batch.write_output(
        peer_job.factors_output_path, output_path(f'{cell_type_name}_weights_file.txt')
    )
    batch.write_output(
        peer_job.factors_output_path,
        output_path(f'{cell_type_name}_precision_file.txt'),
    )
    batch.write_output(
        peer_job.factors_output_path,
        output_path(f'{cell_type_name}_residuals_file.txt'),
    )


def find_cell_types_from_path(path_to_cell_files):
    """
    Get cell types from expression data path
    """

    assert path_to_cell_files.startswith('gs://')
    cs_client = storage.Client()

    bucket_name = path_to_cell_files.split('gs://')[1].split('/')[0]
    bucket = cs_client.get_bucket(bucket_name)
    bucket_path = path_to_cell_files.split(f'gs://{bucket_name}/')[-1]

    logging.info(f'Going to fetch cell types from {bucket_path}')
    blobs = bucket.list_blobs(prefix=bucket_path, delimiter='/')

    cell_types = [
        os.path.basename(b.name)[:-15]
        for b in blobs
        if b.name.endswith('_expression.tsv')
    ]
    logging.info(f'Found {len(cell_types)} cell types: {cell_types}')

    return cell_types


@click.command()
@click.option('--scores-path', help='A path prefix of where PCA scores are located, eg: gs://MyBucket/output-folder/')
@click.option('--covariates-path', help='Path prefix of where covariate files are located, eg: gs://MyBucket/output-folder/',)
@click.option('--sample-id-keys-path', help='Path prefix to convert internal to external, eg: gs://MyBucket/output-folder/')
@click.option('--path-to-cell-files', help='Path prefix of where cell type files are located, eg: gs://MyBucket/output-folder/')
def main(scores_path, covariates_path, sample_id_keys_path, path_to_cell_files):
    """
    This function finds all the cell types, and runs the inner
    workflow for each cell_type, by constructing the path to the
    expression files and the path to the covariates files.


    """

    dataset = os.getenv('CPG_DATASET')
    driver_image = os.getenv('CPG_DRIVER_IMAGE')
    assert driver_image
    backend = hb.ServiceBackend(billing_project=dataset, remote_tmpdir=remote_tmpdir())
    batch = hb.Batch(name='PEER', backend=backend, default_python_image=driver_image)
    cell_types: list = find_cell_types_from_path(path_to_cell_files)

    for cell_type in cell_types:
        expression_file = f'{path_to_cell_files}/{cell_type}_expression.tsv'
        process_cell_type_on_batch(
            batch, 
            cell_type, 
            scores_path=scores_path, 
            covariates_path=covariates_path, 
            sample_id_keys_path=sample_id_keys_path, 
            expression_file=expression_file
        )

    batch.run(wait=False)


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
