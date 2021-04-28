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
	analysis-runner --dataset tob-wgs --output-dir "gs://cpg-tob-wgs-temporary/joint-calling" --description "joint calling" --access-level test joint-calling/workflows/joint_calling.sh --access-level test --callset tob-wgs --version v2 --keep-scratch --reuse

.PHONY: run_full
run_full:
	analysis-runner --dataset tob-wgs --output-dir "gs://cpg-tob-wgs-temporary/joint-calling-full" --description "joint calling" --access-level test joint-calling/workflows/joint_calling.sh --access-level full --batch 0 --batch 1 --callset tob-wgs --version v2 --keep-scratch
