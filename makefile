################################################################
## Management of some basic tasks

usage:
	@echo "Targets"
	@echo "Web    open the github pages for the Web site"
	@echo "repo   open the github repository to access the code"

WEB_URL=
web:
	@echo "Web site: ${WEB_URL}"

REPO=https://github.com/IFB-ElixirFr
repo:
	@echo "Github repo: ${REPO}"
