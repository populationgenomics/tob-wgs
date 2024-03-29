# QC report

Generate a test report and put it into the `test-web` bucket. The report will be available under [test-web.populationgenomics.org.au/tob-wgs/qc/qc-v1.html](https://test-web.populationgenomics.org.au/tob-wgs/qc/qc-v1.html):

```sh
analysis-runner \
   --dataset tob-wgs \
   --access-level test \
   --output-dir "qc-report" \
   --description "QC report test" \
script-for-analysis-runner.sh \
   --batch batch4 \
   --version v2 \
   --analysis-runner
```

A prod run that would put the report into the `main-web` bucket, to be available under [main-web.populationgenomics.org.au/tob-wgs/qc/qc-v1.html](https://main-web.populationgenomics.org.au/tob-wgs/qc/qc-v1.html):

```sh
analysis-runner \
   --dataset tob-wgs \
   --access-level standard \
   --output-dir "qc-report" \
   --description "QC report" \
script-for-analysis-runner.sh \
   --batch batch4 \
   --version v2 \
   --analysis-runner \
   --prod
```

`batch4` means that it will look for the metadata in the `batch4` folder, which includes metadata for all previous batches (1 through 4).
