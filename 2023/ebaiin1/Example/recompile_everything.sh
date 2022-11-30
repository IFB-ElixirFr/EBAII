# Recompile IPython notebooks
cd Session_example_ipynb/ || exit 1
echo "Rebuilding Ipython notebooks"
bash -euiop pipefail process_ipynb.bash > ipynb.log 2>&1

# Recompile R-Markdowns
cd ../Session_example_rmd || exit 1
echo "Rebuilding RMD"
bash -euiop pipefail process_rmd.bash > rmd.log 2>&1

# Recompile HTML
cd ../Session_example_html || exit 1
echo "Rebuilding HTML"
bash -euiop pipefail process_ipynb_to_html.bash > html.log 2>&1

# Re-build book
cd .. || exit 1
echo "Jupyter HTML"
bash -euiop pipefail build_jupyter_book.bash > jupyter.log 2>&1

echo "Process over"