VERSION := v5
TEST_VERSION := v6-4
SCATTER_COUNT_TEST := 10
SCATTER_COUNT_PROD := 100
ANALYSIS_PROJECT := tob-wgs
REUSE_ARG := --reuse

.PHONY: joint_calling_update_submodule
joint_calling_update_submodule:
	-(cd ../joint-calling && git add --all && git commit -m 'WIP' --no-verify && git push)
	(cd joint-calling && git pull --rebase)
	(git add joint-calling && git commit -m 'Update joint-calling submodule' && git push)

.PHONY: joint_calling_test_to_tmp
joint_calling_test_to_tmp:
	analysis-runner \
	--dataset $(ANALYSIS_PROJECT) \
	--output-dir "joint-calling/test-to-tmp" \
	--description "Joint calling test-to-temporary" \
	--access-level test \
	joint-calling/driver_for_analysis_runner.sh workflows/batch_workflow.py \
	--scatter-count $(SCATTER_COUNT_TEST) \
	--from test \
	--to tmp \
	--batch batch1 \
	--analysis-project $(ANALYSIS_PROJECT) \
	--input-project tob-wgs \
	--dataset-version ${TEST_VERSION} \
	--keep-scratch \
	$(REUSE_ARG)

.PHONY: joint_calling_test_to_test
joint_calling_test_to_test:
	analysis-runner \
	--dataset $(ANALYSIS_PROJECT) \
	--output-dir "joint-calling/test-to-test" \
	--description "Joint calling test-to-test" \
	--access-level test \
	joint-calling/driver_for_analysis_runner.sh workflows/batch_workflow.py \
	--scatter-count $(SCATTER_COUNT_TEST) \
	--from test \
	--to test \
	--batch batch1 \
	--analysis-project $(ANALYSIS_PROJECT) \
	--input-project tob-wgs \
	--dataset-version $(VERSION) \
	--keep-scratch \
	$(REUSE_ARG)

.PHONY: joint_calling_main_to_main
joint_calling_main_to_main:
	analysis-runner \
	--dataset $(ANALYSIS_PROJECT) \
	--output-dir "joint-calling/main-to-main" \
	--description "Joint calling main-to-main" \
	--access-level full \
	joint-calling/driver_for_analysis_runner.sh workflows/batch_workflow.py \
	--scatter-count $(SCATTER_COUNT_PROD) \
	--batch batch1 --batch batch2 --batch batch3 --batch batch4 --batch batch5 --batch batch6 --batch batch7 \
	--from main \
	--to main \
	--analysis-project $(ANALYSIS_PROJECT) \
	--input-project tob-wgs \
	--dataset-version $(VERSION) \
	$(REUSE_ARG)
