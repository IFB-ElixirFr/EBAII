How to include IPython notbook in a Jupyter-book: Environment conda:
`mamba env create -f environment.yaml -p ./rmarkdown`

Load conda environment: `conda activate --no-stack ./rmarkdown`

Convert ipython-notebook to rst:
`R --vanilla -e 'rmarkdown::render("Practice.rmd");'` for Practice.rst,
and `R --vanilla -e 'rmarkdown::render("Theory.rmd");'`

Now go to general section to compile the whole book.
