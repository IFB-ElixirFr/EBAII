if [ ! -d ./rmarkdown ]; then
    # Example command line to install env locally:
    echo "Building local conda env"
    mamba env create -f environment.yaml -p ./rmarkdown > rmd.log 2>&1
fi

# Example command line to activate local environment:
echo "Activating local conda env"
conda activate ./rmarkdown

# Example command line to convert Rmarkdown to rst
if [ -f Practice.rst ]; then
    # Remove previous version if exists
    rm --verbose Practice.rst
fi
R --vanilla -e 'rmarkdown::render("Practice.rmd");' >> rmd.log 2>&1

if [ -f Theory.rst ]; then
    # Remove previous version if exists
    rm --verbose Theory.rst
fi
R --vanilla -e 'rmarkdown::render("Theory.rmd");' >> rmd.log 2>&1

echo "RMD over"