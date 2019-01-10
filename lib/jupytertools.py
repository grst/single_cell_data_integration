"""
Collection of helper functions
that are convenient in jupyter
notebooks
"""
import pandas as pd
import os
import logging

PD_DEFAULT_ROWS = 6

pd.set_option('display.max_rows', PD_DEFAULT_ROWS)

def display(obj, n=PD_DEFAULT_ROWS, ncol=10, *args, **kwargs):
    """
    extension of the display function, allows to
    specify the number of rows/cols for pandas.
    """
    from IPython.display import display as ipy_display
    with pd.option_context('display.max_rows', n, 'display.max_columns', ncol):
        ipy_display(obj)


def setwd(max_levels=7):
    """
    set the current working directory
    to the first parent directory that
    contains a Snakefile.

    Args:
        max_levels: If no Snakefile is found after going `max_levels` up,
            an error is raised.

    """

    cnt = 0
    while "Snakefile" not in os.listdir(os.getcwd()):
        if cnt > max_levels:
            raise FileNotFoundError("Snakefile not found in the top {}"
                                    "directories".format(max_levels))
        os.chdir("..")
        cnt += 1

    if cnt:
        print("Changed to directory {}".format("/".join([".."] * cnt)))

    print("Working directory is {}".format(os.path.abspath(os.getcwd())))


def fix_logging():
    """
    Workaround for rstudio/reticulate#386.

    Get all logger instances and add a custom PrintHandler,
    as writing to Stdout doesn't seem to work.
    """
    class PrintHandler(logging.Handler):
        def emit(self, record):
            print(self.format(record))

    for name, logger in logging.Logger.manager.loggerDict.items():
        try:
            logger.handlers = []
            logger.addHandler(PrintHandler())
        except AttributeError:
            pass
