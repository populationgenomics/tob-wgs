#!/bin/bash

set -ex

gsutil mv gs://cpg-tob-wgs-upload/validation.csv $OUTPUT

for sample in ERR1341796 HG01513 HG02238 NA12248 NA12878 NA20502 NA20826
do
    gsutil -m mv "gs://cpg-tob-wgs-upload/$sample.*" $OUTPUT
done
