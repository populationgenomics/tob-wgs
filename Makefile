VERSION := v1

.PHONY: update_joint_calling_submodule
update_joint_calling_submodule:
	-(cd ../joint-calling && git add --all && git commit -m 'WIP' --no-verify && git push)
	(cd joint-calling && git pull --rebase)
	(git add joint-calling && git commit -m 'Update joint-calling submodule' && git push)

.PHONY: joint_calling_test_to_temporary
joint_calling_test_to_temporary:
	analysis-runner \
	--dataset tob-wgs \
	--output-dir "gs://cpg-tob-wgs-hail/joint-calling/test" \
	--description "Joint calling" \
	--access-level test \
	joint-calling/driver_for_analysis_runner.sh workflows/batch_workflow.py \
	--from test \
	--to temporary \
	--callset tob-wgs \
	--version test-$(VERSION) \
	--keep-scratch \
	--reuse

.PHONY: joint_calling_test_to_test
joint_calling_test_to_test:
	analysis-runner \
	--dataset tob-wgs \
	--output-dir "gs://cpg-tob-wgs-hail/joint-calling/test" \
	--description "Joint calling" \
	--access-level standard \
	joint-calling/driver_for_analysis_runner.sh workflows/batch_workflow.py \
	--from test \
	--to test \
	--callset tob-wgs \
	--version $(VERSION) \
	--keep-scratch \
	--reuse

.PHONY: joint_calling_main_to_main
joint_calling_main_to_main:
	analysis-runner \
	--dataset tob-wgs \
	--output-dir "gs://cpg-tob-wgs-hail/joint-calling" \
	--description "Joint calling full" \
	--access-level full \
	joint-calling/driver_for_analysis_runner.sh workflows/batch_workflow.py \
	--batch 1 --batch 2 --batch 3 --batch 4 \
	--from main \
	--to main \
	--callset tob-wgs \
	--version $(VERSION) \
	--reuse
