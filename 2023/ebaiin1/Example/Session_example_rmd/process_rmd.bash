if [ ! -d ./rmarkdown ]; then
    # Example command line to install env locally:
    echo "Building local conda env"
    mamba env create -f environment.yaml -p ./rmarkdown
fi

# Example command line to activate local environment:
echo "Activating local conda env"
conda activate ./rmarkdown

# Example command line to convert Rmarkdown to html
R --vanilla -e 'rmarkdown::render("Practice.rmd");'
R --vanilla -e 'rmarkdown::render("Theory.rmd");'
