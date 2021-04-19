.default: test

.PHONY: test
test:
	jc_pkg jc_submodule sleep jc_test

.PHONY: full
full:
	jc_pkg jc_submodule sleep jc_full

.PHONY: jc_pkg
jc_pkg:
	-(cd ../joint-calling && git add --all && git commit -m 'WIP' && git push)
	cd ../joint-calling && make

.PHONY: jc_submodule
jc_submodule:
	-(cd ../joint-calling && git add --all && git commit -m 'WIP' && git push)
	(cd joint-calling && git pull --rebase)
	(git add joint-calling && git commit -m 'Update joint-calling submodule' && git push)

.PHONY: sleep
sleep:
	sleep 60

.PHONY: jc_test
jc_test:
	analysis-runner --dataset tob-wgs --output-dir "gs://cpg-tob-wgs-temporary/joint-calling" --description "joint calling" --access-level test scripts/drive_joint_calling.py test

.PHONY: jc_full
jc_full:
	analysis-runner --dataset tob-wgs --output-dir "gs://cpg-tob-wgs-temporary/joint-calling-full" --description "joint calling" --access-level test scripts/drive_joint_calling.py full v1
