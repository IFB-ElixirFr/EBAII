################################################################
## Management of some basic tasks

usage:
	@echo "Targets"
	@echo "web    open the github pages for the Web site"
	@echo "repo   open the github repository to access the code"

WEB_URL=https://ifb-elixirfr.github.io/EBAII/2020/
web:
	@echo "Web site"
	@echo "	${WEB_URL}"

REPO=https://github.com/IFB-ElixirFr
repo:
	@echo "Github repo"
	@echo "	${REPO}"
