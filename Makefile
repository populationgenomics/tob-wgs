.default: joint-calling

.PHONY: pkg
pkg:
	-(cd ../joint-calling && git add --all && git commit -m 'WIP' && git push)
	cd ../joint-calling && make

.PHONY: joint-calling
joint-calling:
	-(cd ../joint-calling && git add --all && git commit -m 'WIP' && git push)
	(cd joint-calling && git pull --rebase)
	(git add joint-calling && git commit -m 'Update joint-calling submodule' && git push)
