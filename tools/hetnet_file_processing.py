import os
import pandas as pd
from .processing import head, regularize_colnames

def read_ctd(filename, nrows=None):
    """Read in a file from CTD"""

    # Look at header for the file
    output = head(filename, n_lines=30, print_out=False)
    # Determine which row supplies the column names
    skip_idx = output.split('\n').index('# Fields:')+1

    df = pd.read_csv(filename, skiprows=lambda x: x in list(range(skip_idx))+[skip_idx+1], nrows=nrows)
    df.columns = regularize_colnames(df.columns)
    return df


def read_reactome(filename, nrows=None):
    """Read a reactome mappings file"""
    df = pd.read_csv(filename, sep='\t', header=None, nrows=nrows, dtype=str)
    if df.shape[1] == 8:
        col_names = ['external_id', 'pe_reactome_id', 'pe_name', 'reactome_id', 'reactome_url', 'reactome_name', 'evidence_code', 'taxon']
        df.columns = col_names
    elif df.shape[1] == 6:
        col_names = ['external_id', 'reactome_id', 'reactome_url', 'reactome_name', 'evidence_code', 'taxon']
        df.columns = col_names

    return df

