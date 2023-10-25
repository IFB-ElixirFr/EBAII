
SHELL=/usr/bin/bash
.ONESHELL:
.SHELLFLAGS := -euic
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

.PHONY:
all: SingleCell_GSEA.html SingleCell_DEA.html

SingleCell_GSEA.html: ebaii SingleCell_DEA.html
	mamba activate ./ebaii && \
	R -e 'rmarkdown::render("SingleCell_GSEA.Rmd")'


SingleCell_DEA.html: ebaii Scaled_Normalized_Seurat_Object.RDS
	mamba activate ./ebaii && \
	R -e 'rmarkdown::render("SingleCell_DEA.Rmd")'

Scaled_Normalized_Seurat_Object.RDS:
	wget https://nextcloud.gustaveroussy.fr/s/p8Ab8Be43xFogYN/download/Scaled_Normalized_Seurat_Object.RDS

ebaii:
	mamba env create -p ./ebaii -f environment.yaml

.PHONY:
clean:
	rm -rf ./ebaii *.html
