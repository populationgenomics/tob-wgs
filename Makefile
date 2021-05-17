VERSION := v0

.PHONY: submodule
submodule:
	-(cd ../joint-calling && git add --all && git commit -m 'WIP' --no-verify && git push)
	(cd joint-calling && git pull --rebase)
	(git add joint-calling && git commit -m 'Update joint-calling submodule' && git push)

.PHONY: run_test
run_test:
	analysis-runner \
	--dataset tob-wgs \
	--output-dir "gs://cpg-tob-wgs-hail/joint-vcf/test" \
	--description "Joint calling" \
	--access-level test \
	joint-calling/driver_for_analysis_runner.sh workflows/batch_workflow.py\
    --access-level test \
    --callset tob-wgs \
    --version test-$(VERSION) \
    --keep-scratch \
    --reuse

.PHONY: run_full
run_full:
	analysis-runner \
	--dataset tob-wgs \
	--output-dir "gs://cpg-tob-wgs-hail/joint-vcf" \
	--description "Joint calling" \
	--access-level full \
	joint-calling/driver_for_analysis_runner.sh workflows/batch_workflow.py \
	--access-level full \
	--batch 0 --batch 1 --batch 2 \
	--callset tob-wgs \
	--version $(VERSION)
