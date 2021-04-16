.PHONY: joint-calling
joint-calling:
	-(cd ../joint-calling && git add --all && git commit -m 'WIP' && git push && make)
	(cd joint-calling && git pull --rebase)
	(git add joint-calling && git commit -m 'Update joint-calling submodule' && git push)
