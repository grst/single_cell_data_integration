"""
Module containing helper functions to build Snakemake rules
"""

from snakemake import shell
import os


def _literal_to_r_str(value):
    """Convert a python value to a corresponding R string"""
    _literal_to_str = {
	True: "TRUE",
	False: "FALSE",
	None: "NULL"
    }
    try:
        return _literal_to_str[value]
    except KeyError:
        # quote a string
        if isinstance(value, str):
            return "'{}'".format(value)
        else:
            return str(value)


def render_rmarkdown(input_file, output_file, root_dir, params=None):
    """
    Snakemake wrapper function to render an Rmarkdown document the way I want it to.

    In particular, this function uses bookdown instead of rmarkdown to
    enable figure/table enumeration and allows to pass
    parameters to a parametrized report.

    Args:
        input_file: path to input (Rmd) file
        output_file: path to output (html) file
        root_dir: knitr working directory (python/R will be executed in this directory)
        params: dictionary that will be passed to `params` arg of `rmarkdown::render`.

    """
    param_str = ""
    if params is not None:
        param_str = ", ".join(["{}={}".format(key, _literal_to_r_str(value)) for key, value in params.items()])

    cmd = (
        "MKL_THREADING_LAYER=GNU "  # was necessary to circumvent incompatibilities of Intel mkl with libgomp.
        "Rscript -e \"rmarkdown::render('{input_file}', "
        "   output_file='{output_file}', "
        "   output_format=bookdown::html_document2(), "
        "   knit_root_dir='{root_dir}', "
        "   params = list({params}))\""
    ).format(
        input_file=os.path.abspath(input_file),
        output_file=os.path.abspath(output_file),
        root_dir=os.path.abspath(root_dir),
        params=param_str
    )

    shell(cmd)
