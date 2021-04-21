.PHONY: submodule_test
submodule_test: submodule pkg run_test

.PHONY: submodule_full
submodule_full: submodule pkg run_full

.PHONY: pkg_test
pkg_test: submodule pkg run_test

.PHONY: pkg_full
pkg_full: submodule pkg run_full

.PHONY: pkg
jc_pkg:
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
	analysis-runner --dataset tob-wgs --output-dir "gs://cpg-tob-wgs-temporary/joint-calling" --description "joint calling" --access-level test scripts/drive_joint_calling.py test

.PHONY: run_full
run_full:
	analysis-runner --dataset tob-wgs --output-dir "gs://cpg-tob-wgs-temporary/joint-calling-full" --description "joint calling" --access-level test scripts/drive_joint_calling.py full v1

.PHONY: sleep
sleep:
	sleep 60
