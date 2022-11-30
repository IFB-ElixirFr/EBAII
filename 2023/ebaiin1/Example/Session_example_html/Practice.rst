.. code:: ipython3

    from IPython.display import HTML
    HTML(filename='Practice.html')




.. raw:: html

    <html>
        <h2>How to include HTML in Jupyter-book:</h2>
        <p>Environment conda: <code>mamba env create -f environment.yaml -p ./ipython-notebook</code><br />
        Load conda environment: <code>conda activate --no-stack ./ipython-notebook</code><br />
        Convert ipython-notebook to rst: <code>jupyter nbconvert --to rst Practice.ipynb</code> for Practice.rst, and <code>jupyter nbconvert --to rst Theory.ipynb</code></br />
        Now go to general section to compile the whole book.</p>
    </html>


