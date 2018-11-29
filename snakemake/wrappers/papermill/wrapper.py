"""
Snakemake wrapper to build rmarkdown reports with
bookdown and custom parameters.
"""

__author__    = "Gregor Sturm"
__copyright__ = "Copyright 2018, Gregor Sturm"
__email__     = "g.sturm@tum.de"
__license__   = "MIT"


import papermill as pm
import jupytext as jtx
from os import chdir
import nbformat
from nbconvert import HTMLExporter

snakemake.input.notebook
snakemake.output.report
snakemake.params.root_dir
snakemake.params.nb_params # dict

chdir(snakemake.params.root_dir)

# step 1: if input != ipynb, convert using jupytext
nb = jtx.readf(FILE, format_name=None)
jtx.writef(nb, FILE, format_name=None)

# step 2: execute using papermill
pm.execute_notebook(IN, OUT, parameters={})

# step 3: convert to html using nbconvert
with open(IN) as f:
    nb = nbformat.read(f, as_version=4)

html_exporter = HTMLExporter()
html_exporter.template_file = 'full'

html, resources = html_exporter.from_notebook_node(nb)

with open(OUT, 'w') as f:
    f.write(html)

render_rmarkdown(snakemake.input.script, snakemake.output.report,
        snakemake.params.root_dir, snakemake.params.rmd_params)
