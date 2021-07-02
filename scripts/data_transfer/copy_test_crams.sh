#!/usr/bin/env bash

set -ex

gsutil -m cp "gs://cpg-tob-wgs-archive/cram/batch1/TOB152[0-2].cram*" "gs://cpg-tob-wgs-test/$OUTPUT"
