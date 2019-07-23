import os
import pandas as pd
from copy import deepcopy
from inflection import singularize
from .processing import head, regularize_colnames


CHEM_XREF_KEY_ORDER = ['UNII', 'DRUGBANK', 'MESH', 'CHEMBL.COMPOUND',  'PUBCHEM', 'CompoundCID',
                       'InChIKey', 'NDF-RT', 'RXCUI', 'NCI_THESAURUS', 'CAS', 'IUPHAR',
                       'Iupac', 'KeggId', 'SMILES']


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


# The following functions are for converting the parsed inxight drugs dump into a hetnet
def extract_nodes_from_records(records, node_key, extract_keys):
    nodes = set()
    out_keys = [singularize(k) for k in extract_keys]

    for record in records:
        rec_nodes = record[node_key]

        num_nodes = len(rec_nodes[extract_keys[0]])
        for n in range(num_nodes):
            curr_node = []

            for k in extract_keys:
                curr_node.append(rec_nodes[k][n])
            nodes.add(tuple(curr_node))

    nodes = pd.DataFrame(list(nodes), columns=out_keys)

    nodes['label'] = singularize(node_key).capitalize()
    col_order = ['id', 'name', 'label']
    col_order = [c for c in col_order if c in nodes.columns] + [c for c in nodes.columns if c not in col_order]

    return nodes[col_order]

def records_to_target_nodes(records):

    out = extract_nodes_from_records(records, 'targets', ['ids', 'sources', 'names'])
    return out.dropna(subset=['id', 'source', 'name'], how='all').rename(columns={'source': 'id_source'})

def records_to_condition_nodes(records):

    out = extract_nodes_from_records(records, 'conditions', ['doids', 'meshids', 'names'])
    out = out.dropna(subset=['meshid', 'name', 'doid'], how='all')

    return out.rename(columns={'meshid': 'id'})[['id', 'name', 'label', 'doid']]

def records_to_chem_nodes(records):
    all_c_ids = []
    all_c_names = []
    all_c_sources = []

    for rec in records:
        for o in CHEM_XREF_KEY_ORDER:
            xrefs = rec['xrefs'].get(o, [])
            if len(xrefs) > 0:
                for x in xrefs:
                    all_c_ids.append(x)
                    all_c_names.append('GET')
                    all_c_sources.append(o)
                break

    out = pd.DataFrame({'id': all_c_ids, 'name': all_c_names, 'source': all_c_sources})
    out['label'] = 'Compound'
    return out

def extract_edges_from_records(records, edge_key, extract_keys, is_empty=lambda f: False):
    out_data = []
    out_keys = [singularize(k) for k in extract_keys]

    for rec in records:
        rec_edges = rec[edge_key]

        # Make sure there's an edge
        if is_empty(rec_edges):
            continue

        # Get the proper X-ref
        for o in CHEM_XREF_KEY_ORDER:
            xrefs = rec['xrefs'].get(o, [])
            if len(xrefs) > 0:
                break

        # Grab the important info
        num_edges = len(rec_edges[extract_keys[0]])
        for n in range(num_edges):
            curr_edge = dict()
            for k, ok in zip(extract_keys, out_keys):
                curr_edge[ok] = rec_edges[k][n]
            for xref in xrefs:
                curr_edge['comp_id'] = xref
                out_data.append(deepcopy(curr_edge))

    return pd.DataFrame(out_data)


def is_empty_condition(condition):
    return all([x is None for x in condition['meshids']]) and \
        all([x is None for x in condition['doids']]) and \
        all([x is None for x in condition['names']])


def is_empty_target(target):
    return all([x is None for x in target['ids']]) or all([x == "Unknown" for x in target['ids']])


def records_to_treats_edges(records):

    extract_keys = ['doids', 'meshids', 'names', 'phases', 'modalities', 'dates']
    out = extract_edges_from_records(records, 'conditions', extract_keys, is_empty_condition)

    return out


def records_to_target_edges(records):

    out = extract_edges_from_records(records, 'targets', ['ids', 'pharmas'], is_empty_target)

    return out.rename(columns={'id': 'target_id', 'pharma': 'interaction'})

# The following functions are for interfacing with the chembl targets api.
def query_chembl(ids, offset, limit):

        query_url = 'https://www.ebi.ac.uk/chembl/api/data/target/set/{}?format=json'
        this_q = query_url.format(';'.join(chembl_targets[offset:offset+limit]))
        r = requests.get(this_q)
        return json.loads(r.text)['targets']


def download_chembl_targets(target_ids):
    limit = 100
    all_res = []

    for i in tqdm(range(len(target_ids) // 100)):
        offset = i*limit
        all_res += query_chembl(target_ids, offset, limit)
        time.sleep(1)
    offset = (i+1)* limit
    limit = len(target_ids) - offset
    all_res += query_chembl(target_ids, offset, limit)

    return all_res

def process_chembl_targets_api(chembl_results):

    target_keys = ['organism', 'target_chembl_id', 'tax_id', 'target_type', 'pref_name']
    component_keys = ['component_type', 'component_id']
    xref_keys = ['xref_name', 'xref_src_db', 'xref_id']

    target_xrefs = []
    for elem in chembl_results:
        this_xref = dict()
        for tk in target_keys:
            this_xref[tk] = str(elem[tk])

        tcs = elem['target_components']

        if len(tcs) == 0:
            for xk in xref_keys:
                this_xref[xk] = float('nan')
            for ck in component_keys:
                this_xref[ck] = float('nan')
            target_xrefs.append(this_xref)

        for tc in tcs:
            xrefs = tc['target_component_xrefs']
            for ck in component_keys:
                this_xref[ck] = str(tc[ck])

            if len(xrefs) == 0:
                for xk in xref_keys:
                    this_xref[xk] = float('nan')
                target_xrefs.append(deepcopy(this_xref))

            for xr in xrefs:
                for xk in xref_keys:
                    this_xref[xk] = str(xr[xk])
                target_xrefs.append(deepcopy(this_xref))
    return pd.DataFrame(target_xrefs)


def select_chembl_target_xref(chembl_target_df):
    # Only get the following xrefs for chembl targets
    chembl_protein_xrefs = ['EnsemblGene', 'UniProt']
    non_protein_target_types = ['ORGANISM','CELL-LINE','NUCLEIC-ACID','METAL','SMALL MOLECULE',
                                'SUBCELLULAR','UNKNOWN','MACROMOLECULE']

    # Filter by the xrefs
    chembl_info = (chembl_target_df.query('xref_src_db in @chembl_protein_xrefs')
                                   .sort_values('xref_src_db')
                                   .drop_duplicates(subset=['component_id', 'target_chembl_id'])
                                   .sort_index())
    # Add back in the non-protein targets
    chembl_info = pd.concat([chembl_info, chembl_target_df.query('target_type in @non_protein_target_types')],
                            ignore_index=True)

    return chembl_info


def determine_evidence(code):
    """
    Computationaly derived codes are taken from GO's webiste:
    http://geneontology.org/docs/guide-go-evidence-codes/

    :param code: str, 3 (or 2) letter code for annotation evidence.
    :retrun: str, 'curated' or 'computed' status of the code.
    """
    comp_codes = ['ISS', 'ISO', 'ISA', 'ISM', 'IGC', 'RCA', 'IEA']
    return 'computed' if code.upper() in comp_codes else 'curated'

