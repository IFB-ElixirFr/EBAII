if [ ! -d ./jupyter-book ]; then
    # Example of command line used to install jupyter-book
    echo "Creating jupyter-book environment as ${PWD}/jupyter-book"
    mamba env create -f ./environment.yaml -p ./jupyter-book
fi

# Activating jupyterbook environment
echo "Activating conda environment"
conda activate --no-stack ./jupyter-book

if [ -d _build ]; then
    # Cleaning previous version before re-building html
    echo "Removing old version before rendering"
    rm -r _build
fi

echo "Build jupyter-book"
jupyter-book build --builder singlehtml "${PWD}" --toc _toc.yml --config _config.yml