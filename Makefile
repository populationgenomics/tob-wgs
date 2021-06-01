VERSION := v2

.PHONY: update_joint_calling_submodule
update_joint_calling_submodule:
	-(cd ../joint-calling && git add --all && git commit -m 'WIP' --no-verify && git push)
	(cd joint-calling && git pull --rebase)
	(git add joint-calling && git commit -m 'Update joint-calling submodule' && git push)

.PHONY: joint_calling_test_to_temporary
joint_calling_test_to_temporary:
	analysis-runner \
	--dataset tob-wgs \
	--output-dir "gs://cpg-tob-wgs-hail/joint-calling/test-to-temporary" \
	--description "Joint calling test-to-temporary" \
	--access-level test \
	joint-calling/driver_for_analysis_runner.sh workflows/batch_workflow.py \
	--from test \
	--to test-tmp \
	--callset tob-wgs \
	--version $(VERSION) \
	--keep-scratch \
	--reuse

.PHONY: joint_calling_test_to_test
joint_calling_test_to_test:
	analysis-runner \
	--dataset tob-wgs \
	--output-dir "gs://cpg-tob-wgs-hail/joint-calling/test-to-test" \
	--description "Joint calling test-to-test" \
	--access-level full \
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
	--output-dir "gs://cpg-tob-wgs-hail/joint-calling/main-to-main" \
	--description "Joint calling main-to-main" \
	--access-level full \
	joint-calling/driver_for_analysis_runner.sh workflows/batch_workflow.py \
	--batch batch1 --batch batch2 --batch batch3 --batch batch4 \
	--from main \
	--to main \
	--callset tob-wgs \
	--version $(VERSION) \
	--reuse
