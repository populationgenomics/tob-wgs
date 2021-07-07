# Exploration of TOB-WGS SNPchip IDs

Here we use the analysis-runner to render an RMarkdown notebook
into an HTML file. The command to generate the report is:

```bash
analysis-runner \
   --dataset "tob-wgs" \
   --access-level "test" \
   --output-dir "pdiakumis/snpchip_id_report" \
   --description "snpchip id report exploration" \
   run.R
```

Example output can be found at:
<https://test-web.populationgenomics.org.au/tob-wgs/pdiakumis/snpchip_id_report/snpchip_ids.html>

To run using inputs from the `main` bucket, use `--access-level "full"` in the above command.
