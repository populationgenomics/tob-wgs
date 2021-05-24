# QC report

Generate a report and put it into the web bucket. Will make a test report into the temporary bucket first.

```sh
analysis-runner --dataset tob-wgs --access-level standard --output-dir "gs://cpg-tob-wgs-temporary/qc_report" --description "QC report" \
script-for-analysis-runner.sh --prod --batch batch2
```
