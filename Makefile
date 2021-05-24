VERSION := v2

.PHONY: pkg
pkg:
	-(cd ../joint-calling && git add --all && git commit -m 'WIP' --no-verify && git push)
	cd ../joint-calling && make
	sleep 60

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
	--description "joint calling" \
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
	--description "joint calling" \
	--access-level test \
	joint-calling/driver_for_analysis_runner.sh workflows/batch_workflow.py \
	--access-level full \
	--batch 0 --batch 1 \
	--callset tob-wgs \
	--version $(VERSION) \
	--keep-scratch --reuse

.PHONY: run_unit_test
run_unit_test:
	analysis-runner \
	--dataset tob-wgs \
	--access-level test \
	--output-dir "gs://cpg-tob-wgs-hail/joint-vcf/unit-test" \
	--description "test evaluation" \
	joint-calling/workflows/joint_calling.sh test/test_rf_evaluation.py
