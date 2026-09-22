
SHELL=/usr/bin/bash
.ONESHELL:
.SHELLFLAGS := -euic
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

.PHONY:
all: SingleCell_GSEA.html SingleCell_DEA.html

SingleCell_GSEA.html: ebaii SingleCell_DEA.html
	conda activate ./ebaii && \
	R -e 'rmarkdown::render("SingleCell_GSEA.Rmd")'


SingleCell_DEA.html: ebaii 
	conda activate ./ebaii && \
	R -e 'rmarkdown::render("SingleCell_DEA.Rmd")'


ebaii:
	conda create -p ./ebaii --file environment.linux-64.pin.txt

.PHONY:
clean:
	rm -rf ./ebaii *.RDS
