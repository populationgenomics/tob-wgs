#!/usr/bin/env python3

"""Reproduction case for the Hail team."""

import hail as hl
from cpg_utils.hail_batch import dataset_path, init_batch

MT_PATH = dataset_path('mt/v7.mt')
NUM_PARTITIONS = 10000


def main():
    """Main entrypoint."""
    init_batch()
    mt = hl.read_matrix_table(MT_PATH)
    mt = mt.naive_coalesce(NUM_PARTITIONS)
    print(mt.rows()._force_count())  # pylint: disable=protected-access


if __name__ == '__main__':
    main()
