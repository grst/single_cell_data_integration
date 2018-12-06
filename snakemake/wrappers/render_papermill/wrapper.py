"""
Snakemake wrapper to build reports from
jupyter notebooks using papermill.
"""

__author__    = "Gregor Sturm"
__copyright__ = "Copyright 2018, Gregor Sturm"
__email__     = "g.sturm@tum.de"
__license__   = "MIT"


import papermill as pm
import jupytext as jtx
from os import chdir
from os.path import splitext
import nbformat
from tempfile import NamedTemporaryFile
from nbconvert import HTMLExporter

chdir(snakemake.params.root_dir)

input_notebook = snakemake.input.notebook
tmp_notebook = NamedTemporaryFile(suffix=".ipynb")
tmp_notebook_executed = NamedTemporaryFile(suffix=".ipynb")

if splitext(input_notebook)[1] != ".ipynb":
    # step 1: if input != ipynb, convert using jupytext
    nb = jtx.readf(input_notebook, format_name=None)
    jtx.writef(nb, tmp_notebook.name, format_name=None)
    input_notebook = tmp_notebook.name


# step 2: execute using papermill
pm.execute_notebook(input_notebook, tmp_notebook_executed.name, parameters=snakemake.params.nb_params, log_output=True)


# step 3: convert to html using nbconvert
with open(tmp_notebook_executed.name) as f:
    nb = nbformat.read(f, as_version=4)

html_exporter = HTMLExporter()
html_exporter.template_file = 'full'

html, resources = html_exporter.from_notebook_node(nb)

with open(snakemake.output.report, 'w') as f:
    f.write(html)

tmp_notebook.close()
tmp_notebook_executed.close()
