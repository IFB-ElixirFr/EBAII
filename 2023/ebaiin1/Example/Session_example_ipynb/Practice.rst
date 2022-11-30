Pure Ipynb: practice part
=========================

Inclusion of IPython notebook is quite straightforward.

Environment conda:
``mamba env create -f environment.yaml -p ./ipython-notebook``

Load conda environment: ``conda activate --no-stack ./ipython-notebook``

Convert ipython-notebook to rst:
``jupyter nbconvert --to rst Practice.ipynb`` for Practice.rst, and
``jupyter nbconvert --to rst Theory.ipynb``

Now go to general section to compile the whole book.

PS:
^^^

Note you can activate a conda environment within a notebook. See on
`StackOverflow <https://stackoverflow.com/questions/74597051/activate-conda-environment-inside-r-script>`__.
But I would highly discourage it.
