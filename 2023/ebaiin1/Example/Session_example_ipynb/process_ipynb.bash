
if [ ! -d ipython-notebook ]; then
    # Example command line to install env locally:
    echo "Installing dependencies in ${PWD}/ipython-notebook"
    mamba env create -f environment.yaml -p ./ipython-notebook
fi 

# Example command line to activate local environment:
echo "Activating ${PWD}/ipython-notebook environment:"
conda activate --no-stack ./ipython-notebook


# Example command line to convert ipython notebooks to rst
echo "Converting ipynb to rst"
jupyter nbconvert --to rst Practice.ipynb
jupyter nbconvert --to rst Theory.ipynb