import os
import subprocess
import numpy as np
import pandas as pd
from copy import deepcopy
from itertools import chain
from collections import defaultdict

def head(file_name, n_lines=20, print_out=True, line_nums=False):
    """Get the first n_lines lines of a file. Print if print_out=True, else return as list. Works on UNIX systems"""

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


def map_to_parent(df, parent_map, col_name):
    """Map a column to its parent value"""
    df_out = df.copy()
    out_name = 'parent_'+col_name
    df_out[out_name] = df_out[col_name].map(parent_map)
    return df_out


def make_parent_map(df, child_id_col, parent_id_col, name_col):
    """
    Takes dataframe with identifiers for parent and child, and name for child and produces a map dict from
    child name to parent name
    """
    map_df = pd.merge(df, df[[name_col, child_id_col]], how='left', left_on=parent_id_col, right_on=child_id_col,
                      suffixes=('', '_parent'))

    par_name_col = name_col+'_parent'
    # Fills in root nodes with its own name
    map_df[par_name_col] = map_df[par_name_col].fillna(map_df[name_col])
    return map_df.set_index(name_col)[par_name_col].to_dict()


def prepend_direction_to_map(mapper, directions=iter(()), char=''):
    """
    In a dict, will prepend all elements of `directions` to all keys and values, sparated by `char`

    e.g. a dict of  {'metabolic processing': 'metabolic processing',
                     'acetylation': 'metabolic processing'}
    And directions ['increasing', 'decreasing'] with char '^'

    becomes:
         {'increases^metabolic processing': 'increases^metabolic processing',
          'increases^acetylation': 'increases^metabolic processing',
          'decreases^metabolic processing': 'decreases^metabolic processing',
          'decreases^acetylation': 'decreases^metabolic processing'}
    """
    out_map = dict()
    for k, v in mapper.items():
        for d in directions:
            new_k = d+char+k
            new_v = d+char+v
            out_map[new_k] = new_v
    return out_map


def get_parent_to_child_map(df, parent_col, child_col, sep='|', verbose=False):
    """
    Gets a map from the parent ID to IDs of all children (and children of children etc.)
    Source dataframe (EG, CTD)
    """
    # Initialize the map
    par_to_children = defaultdict(set)
    expanded = expand_col_on_char(df.dropna(subset=[parent_col]), parent_col, sep)

    # Get mappings from child to local parents
    for row in expanded.itertuples():
        par_to_children[getattr(row, parent_col)].add(getattr(row, child_col))

    # Extend the map... parents map to local children, but need mapping to the children of children, too
    par_to_children = extend_map(par_to_children, verbose)

    # Discard empty values (default_dict will add an empty one when simply accessed)
    par_to_children = {k: v for k, v in par_to_children.items() if len(v) > 0}

    return par_to_children


def extend_map(mapper, verbose=False):
    """
    Takes a dict mapper of a tree (e.g. parents to children) a extends down so parents map to all
    intermediate children, terminating at a leaf.
    """

    prev_total_mappings = 0
    total_mappings = len(list(chain(*mapper.values())))
    num_iter = 0
    i = 0

    while (total_mappings != prev_total_mappings):
        prev_total_mappings = len(list(chain(*mapper.values())))
        current_map = deepcopy(mapper)
        for k, v in current_map.items():
            for val in v:
                mapper[k].update(mapper[val])
        total_mappings = len(list(chain(*mapper.values())))


        if verbose:
            num_iter += total_mappings
            # Print an update every 10 Million iterations (approx)
            if i == 0 or np.floor(num_iter / 1e7) - np.floor((num_iter - total_mappings) / 1e7) > 0:
                print('Iter: {}, Total Mappings: {:,}'.format(i, total_mappings))
        i += 1

    return mapper


def char_combine_col(col, char='|'):
    """
    Converts a column to a cell, splitting elements within that column along a character, dedupcating all elements,
    then joining across that character on that charcter.

    e.g. 1   "123456|102934"
         2   "123456"
         3   "102934|432452"
         4   "201945"

         becomes:
            "123456|102934|432452|201945"
    """
    elems = []
    for elem in col:
        elems += str(elem).split(char)
    return char.join(list(set(elems)))


def char_combine_dataframe_rows(df, char='|'):
    return df.apply(lambda col: char_combine_col(col, char))


def combine_group_rows_on_char(df, group_on, combine_cols=None, char='|'):
    """
    Performs a Groupby on a dataframe and then converts each group into a single row, joinned by a character `char`

    Primarly suppports grouping on columns, other methods have not been tested.

    :param df:  The dataframe to group
    :param group_on: the column name or list of column names to group by
    :param combine_cols: a list of column names to combine with a character, if None, will combine all columns.
        can save computation time to provide only the columns of interest for combination
    :param char: the character to combine the columns with. Defaults to a `|` character.

    :return: Dataframe with 1 row per group, and information of different rows joined by given character.

    """
    col_order = df.columns

    if type(group_on) in (str, int, float):
        group_on = [group_on]

    grouped = df.groupby(group_on)

    if combine_cols is None:
        combine_cols = find_cols_with_multi_values(grouped)

    out_df = grouped.first()
    for col in combine_cols:
        out_df[col] = grouped[col].apply(char_combine_col, char=char)

    return out_df.reset_index()[col_order]


def find_cols_with_multi_values(grouped):
    """
    In a Pandas Groupby object, determines which columns have more than one value per group
    """
    multi_val_cols = (grouped.nunique() > 1).sum()
    multi_val_cols = multi_val_cols[multi_val_cols > 0].index.tolist()
    return multi_val_cols


def make_child_to_root_map(df, parent_col, child_col, sep='|', verbose=False):
    """
    In a dataframe where one column is child identfiers, and another parent,
    creates a map from child to root.
    """

    root_concepts = df[df[parent_col].isnull()][child_col].tolist()
    parent_to_child = get_parent_to_child_map(df, parent_col, child_col, sep, verbose)

    child_to_parent = dict()

    for k, v in parent_to_child.items():
        for val in v:
            if k in root_concepts:
                child_to_parent[val] = k

    for r in root_concepts:
        child_to_parent[r] = r

    return child_to_parent


def convert_abbrev_mapper_to_full(mapper, map_df, abbrev, full):
    """
    With a mapping file from abbrev to abbrev, and a datafame that has abbrev to full mappings, produces
    a mapper from full to full.
    """
    abbrev_to_name = map_df.set_index(abbrev)[full].to_dict()
    return {abbrev_to_name[k]: abbrev_to_name[v] for k, v in mapper.items()}

