"""
Collection of helper functions
that are convenient in jupyter
notebooks
"""
from IPython.display import display as ipy_display
import pandas as pd

PD_DEFAULT_ROWS = 6

pd.set_option('display.max_rows', PD_DEFAULT_ROWS)

def display(obj, n=PD_DEFAULT_ROWS, ncol=10, *args, **kwargs):
    with pd.option_context('display.max_rows', n, 'display.max_columns', ncol):
        ipy_display(obj)


