#!/usr/bin/env bash

set -ex

gsutil cp "gs://cpg-tob-wgs-main-web/pdiakumis/snpchip_id_report/snpchip_ids.html" "gs://cpg-tob-wgs-test-tmp/$OUTPUT"
