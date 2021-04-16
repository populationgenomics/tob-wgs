.PHONY: joint-calling
joint-calling:
	cd ../joint-calling
	git add --all
	git commit -m 'WIP'
	git push
	cd ../tob-wgs/joint-calling
	git pull --rebase
	git push
	cd ..
	git add joint-calling
	git commit -m 'Update joint-calling submodule'
	git push
