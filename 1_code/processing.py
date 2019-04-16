import os
import subprocess
import numpy as np
import pandas as pd

def head(file_name, n_lines=20, print_out=True, line_nums=False):
    """Get the first n_lines lines of a file. Print if print_out=True, else return as list"""

    assert type(n_lines) == int
    n_lines = str(n_lines)

    if os.path.splitext(file_name)[-1] in ['.gz', '.zip']:
        zcat_p = subprocess.Popen(['zcat', file_name], stdout=subprocess.PIPE)
        head_p = subprocess.Popen(['head', '-n', n_lines], stdin=zcat_p.stdout, stdout=subprocess.PIPE)
        zcat_p.stdout.close()
        output = head_p.communicate()[0].decode('utf-8')
    else:
        output = subprocess.Popen(['head', '-n', n_lines, file_name], stdout=subprocess.PIPE).stdout.read().decode('utf-8')

    if not print_out:
        return output

    if line_nums:
        pad = int(np.log10(int(n_lines)))+1
        fmt_str = '{:' + '{}'.format(pad) + '.0f}:'

        split_out = output.rstrip('\n').split('\n')
        for i, line in enumerate(split_out):
            print(fmt_str.format(i), line)
    else:
        print(output)


def cammel_to_underscore(name):
    """Converts CammelCaseNames to underscore_separated_names"""
    # Get the indicies of capital letters
    idx_to_change = []
    for i, letter in enumerate(name):
        if letter != letter.lower():
            idx_to_change.append(i)

    # Ensure first section is grabbed when starting with lower case e.g. cammelCaseExamples
    idx_to_change = [0] + idx_to_change if 0 not in idx_to_change else idx_to_change

    # If strings of caplital letters in a reow (e.g. OrchidIDs)
    # Start with only the first of the string and remove the rest...
    prev_idx = 0
    to_remove = []
    for idx in idx_to_change:
        if idx - prev_idx == 1:
            to_remove.append(idx)
        prev_idx = idx
    idx_to_change = [x for x in idx_to_change if x not in to_remove]

    # Build the new name
    out_name = ''
    for i, start_idx in enumerate(idx_to_change):
        if i+1 < len(idx_to_change):
            end_idx = idx_to_change[i+1]
            out_name += name[start_idx:end_idx].lower() + '_'
        else:
            out_name += name[start_idx:].lower()
    if not out_name:
        return name
    return out_name


def remove_lead_chars(name):
    """Strips whitespace and special characters from the start of a string"""
    start = 0
    while not name[start].isalpha():
        start+=1
    return name[start:]


def remove_end_chars(name):
    """Strips whitespace and special characters from the end of a string"""
    end = len(name)
    while not name[end-1].isalnum():
        end-=1
    return name[:end]


def strip_special_chars(name):
    """Removes most special characters from a string, keeps [' ', '-', or '_']"""
    return ''.join([c for c in name if c.isalnum() or c in [' ', '-', '_']])


def regularize_colname(name):
    """Regularize a string to a Pandas queryable name, with underscore separation"""
    out_name = remove_lead_chars(name)
    out_name = remove_end_chars(out_name)
    out_name = strip_special_chars(out_name)
    if out_name.isidentifier():
        out_name = cammel_to_underscore(out_name)
    out_name = out_name.replace(' ', '_').replace('-', '_')
    return out_name.lower()


def regularize_colnames(col_names):
    """Regularize a list of column names. See regularize_colname"""
    return [regularize_colname(c) for c in col_names]


def split_col(col, char):
    """Splits a Series by a characters, into a new series of equal length with each element as a list"""
    return col.apply(lambda s: s.split(char))


def expand_split_col(col_split):
    """
    Expands a Series that contains a list into a new DataFrame with single item per row, and the original index
    contained in a new column 'old_idx'

    Example

        1 [5, 4, 3]
        3 [2, 5]
        5 [1]
        name='numbers'

    would become
        old_idx, numbers
        1, 5
        1, 4
        1, 3
        3, 2
        3, 5
        5, 1
    """
    col_name = col_split.name

    col_split = col_split.to_frame().reset_index()

    old_idx = []
    new_col = []
    for row in col_split.itertuples():
        for item in row[2]:
            old_idx.append(row[1])
            new_col.append(item)

    return pd.DataFrame({'old_idx': old_idx, col_name: new_col})


def expand_col_on_char(df, col_name, char):
    """
    Expands rows in a dataframe due to a column where elements contain multiple values separated by a character,
    resulting in only one element per row in the new column


    For example, a column of pipe separated Pubmed IDs (in .csv format):
        letter, pmid
        a, 12805067|17698565|18279049|21329777
        b, 10072544|12721113

    Would become:
        letter, pmid
        a, 12805067
        a, 17698565
        a, 18279049
        a, 21329777
        b, 10072544
        b, 12721113
    """

    # Copy df and get column order
    df_out = df.copy()
    col_order = df_out.columns

    # Split the desired column on the desired character and expand the rows
    col_split = split_col(df_out[col_name], char)
    col_split = expand_split_col(col_split)

    # Make the old index avaliable for merging
    df_out = df_out.reset_index()
    df_out = df_out.rename(columns={'index': 'old_idx'})
    df_out = pd.merge(df_out.drop(col_name, axis=1), col_split, on='old_idx', how='outer')
    return df_out.drop('old_idx', axis=1)[col_order]


