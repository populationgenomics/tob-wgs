#!/usr/bin/env python3

"""Reproduction case for the Hail team."""

from cpg_utils.hail_batch import (
    dataset_path,
    init_batch,
    output_path
)

MT_PATH = dataset_path('mt/v7.mt')
NUM_PARTITIONS = 10000

def main():
    init_batch()
    mt = hl.read_matrix_table(MT_PATH)
    intervals = mt._calculate_new_partitions(NUM_PARTITIONS)
    intervals_path = output_path(f'eqtl/debug/intervals_{NUM_PARTITIONS}', 'analysis')
    hl.experimental.write_expression(intervals, intervals_path)
    print(hl.read_matrix_table(MT_PATH, _intervals=intervals).rows()._force_count())

if __name__ == '__main__':
    main()
