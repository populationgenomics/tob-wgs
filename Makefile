.default: test

.PHONY: test
test: jc_submodule jc_test

.PHONY: full
full: jc_submodule jc_full

.PHONY: pkg
jc_pkg:
	-(cd ../joint-calling && git add --all && git commit -m 'WIP' && git push)
	cd ../joint-calling && make
	sleep 60

.PHONY: jc_submodule
jc_submodule:
	-(cd ../joint-calling && git add --all && git commit -m 'WIP' && git push)
	(cd joint-calling && git pull --rebase)
	(git add joint-calling && git commit -m 'Update joint-calling submodule' && git push)

.PHONY: jc_test
jc_test:
	analysis-runner --dataset tob-wgs --output-dir "gs://cpg-tob-wgs-temporary/joint-calling" --description "joint calling" --access-level test scripts/drive_joint_calling.py test

.PHONY: jc_full
jc_full:
	analysis-runner --dataset tob-wgs --output-dir "gs://cpg-tob-wgs-temporary/joint-calling-full" --description "joint calling" --access-level test scripts/drive_joint_calling.py full v1

.PHONY: sleep
sleep:
	sleep 60
