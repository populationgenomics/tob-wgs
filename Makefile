VERSION := v2
SCATTER_COUNT_TEST := 10
SCATTER_COUNT_PROD := 100
CALLSET := tob-wgs
REUSE_ARG := --reuse

.PHONY: joint_calling_update_submodule
joint_calling_update_submodule:
	-(cd ../joint-calling && git add --all && git commit -m 'WIP' --no-verify && git push)
	(cd joint-calling && git pull --rebase)
	(git add joint-calling && git commit -m 'Update joint-calling submodule' && git push)

.PHONY: joint_calling_clean_tmp
joint_calling_clean_tmp:
	gsutil -q rm -rf gs://cpg-tob-wgs-test-tmp/joint-calling/$(VERSION)

.PHONY: joint_calling_clean_test
joint_calling_clean_test:
	gsutil -q rm -rf gs://cpg-tob-wgs-test/joint-calling/$(VERSION)

.PHONY: joint_calling_test_to_tmp
joint_calling_test_to_tmp:
	analysis-runner \
	--dataset $(CALLSET) \
	--output-dir "gs://cpg-tob-wgs-hail/joint-calling/test-to-tmp" \
	--description "Joint calling test-to-temporary" \
	--access-level test \
	joint-calling/driver_for_analysis_runner.sh workflows/batch_workflow.py \
	--scatter-count $(SCATTER_COUNT_TEST) \
	--from test \
	--to tmp \
	--batch batch1 \
	--callset $(CALLSET) \
	--version v2 \
	--keep-scratch \
	$(REUSE_ARG)

.PHONY: joint_calling_test_to_tmp_extend
joint_calling_test_to_tmp_extend:
	analysis-runner \
	--dataset $(CALLSET) \
	--output-dir "gs://cpg-tob-wgs-hail/joint-calling/test-to-tmp_extend" \
	--description "Joint calling test-to-temporary" \
	--access-level test \
	joint-calling/driver_for_analysis_runner.sh workflows/batch_workflow.py \
	--scatter-count $(SCATTER_COUNT_TEST) \
	--existing-mt "gs://cpg-tob-wgs-test/mt/v2-raw.mt" \
	--from test \
	--to tmp \
	--batch batch2 \
	--callset $(CALLSET) \
	--version v3 \
	--keep-scratch \
	$(REUSE_ARG)

.PHONY: joint_calling_test_to_test
joint_calling_test_to_test:
	analysis-runner \
	--dataset $(CALLSET) \
	--output-dir "gs://cpg-tob-wgs-hail/joint-calling/test-to-test" \
	--description "Joint calling test-to-test" \
	--access-level test \
	joint-calling/driver_for_analysis_runner.sh workflows/batch_workflow.py \
	--scatter-count $(SCATTER_COUNT_TEST) \
	--from test \
	--to test \
	--callset $(CALLSET) \
	--version $(VERSION) \
	--keep-scratch \
	$(REUSE_ARG)

.PHONY: joint_calling_main_to_main
joint_calling_main_to_main:
	analysis-runner \
	--dataset $(CALLSET) \
	--output-dir "gs://cpg-tob-wgs-hail/joint-calling/main-to-main" \
	--description "Joint calling main-to-main" \
	--access-level full \
	joint-calling/driver_for_analysis_runner.sh workflows/batch_workflow.py \
	--scatter-count $(SCATTER_COUNT_PROD) \
	--batch batch1 --batch batch2 --batch batch3 --batch batch4 \
	--from main \
	--to main \
	--callset $(CALLSET) \
	--version $(VERSION) \
	$(REUSE_ARG)
