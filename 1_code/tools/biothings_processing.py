import biothings_client
import pandas as pd

def download_mychem_ids(chem_ids):
    mc = biothings_client.get_client('chem')

    fields = ['aeolus.drug_name', 'chebi.name', 'chebi.id', 'chembl.pref_name', 'chembl.molecule_chembl_id',
              'drugbank.name', 'drugbank.id ', 'drugcentral.pharmacology_class.mesh_pa.code ',
              'drugcentral.xrefs.mesh_descriptor_ui', 'drugcentral.xrefs.mesh_supplemental_record_ui',
              'ginas.xrefs.MESH', 'ginas.preferred_name', 'unii.unii', 'unii.preferred_term', 'pharmgkb.name',
              'pharmgkb.xrefs.mesh', 'drugbank.targets', 'drugcentral.bioactivity']


    res = mc.getchems(chem_ids, fields=fields)
    return res


def get_dotval(item, dotfield, list_ok=True):
    fields = dotfield.split('.')
    for field in fields:
        if type(item) == dict:
            item = item.get(field, dict())
        elif type(item) == list:
            item = [i.get(field, dict()) for i in item]

    if item != dict() and type(item) != list:
        return item
    elif type(item) == list and not all([x == dict() for x in item]):
        if len(item) == 1:
            return item[0]
        elif not list_ok:
            # Uniquify and resturn as | joined string
            return '|'.join(list(set([str(x) for x in item])))
        return [float('nan') if i == dict() else i for i in item]
    else:
        return float('nan')

def get_best_name(res_entry):
    names = ['unii.preferred_term', 'chembl.pref_name', 'chebi.name', 'drugbank.name',
         'aeolus.drug_name', 'pharmgkb.name', 'ginas.preferred_name']

    for n in names:
        out_name = get_dotval(res_entry, n)
        if type(out_name) == str:
            return out_name.capitalize()

    return float('nan')

def get_xrefs(res_entry):
    other_xrefs = ['chebi.id',  'chembl.molecule_chembl_id',  'drugbank.id', 'unii.unii']
    mesh_xrefs = ['drugcentral.xrefs.mesh_descriptor_ui',   'drugcentral.xrefs.mesh_supplemental_record_ui',
                     'ginas.xrefs.MESH',  'pharmgkb.xrefs.mesh']
    out = dict()

    out['q'] = res_entry['query']

    if not res_entry.get('notfound', False):
        res_id = res_entry.get('_id')
        # Definition of inchi-ikey is 14-10-1
        if tuple(len(x) for x in res_id.split('-')) == (14, 10, 1):
            out['ikey'] = res_id
        else:
            out['ikey'] = float('nan')

    for o in other_xrefs:
        source = o.split('.')[0] + '_id'
        out[source] = get_dotval(res_entry, o)

    mesh = []
    for m in mesh_xrefs:
        mesh_ref = get_dotval(res_entry, m)
        if type(mesh_ref) == str:
            mesh.append(mesh_ref)
        elif type(mesh_ref) == list:
            mesh += mesh_ref

    if len(mesh) == 0:
        out['mesh_id'] = float('nan')
    elif len(set(mesh)) == 1:
        out['mesh_id'] = 'MESH:'+mesh[0]
    else:
        out['mesh_id'] = list(set(['MESH:'+ m for m in mesh]))

    return out

def join_list_char(l, c='|'):
    if type(l) == list:
        return c.join([str(e) for e in l])
    else:
        return l

def col_lists_to_str(col, c='|'):
    return col.apply(lambda e: join_list_char(e, c))

def make_mg_uniprot_map(mg_res):
    # Some one to many found, so we will de-duplicate
    uniprot_map = {x['query']: x.get('entrezgene', float('nan')) for x in mg_res['out']}

    # Simple strategy of just taking the lowest numbered cross reference
    dup_xref = [x[0] for x in mg_res['dup']]
    for x in mg_res['out']:
        # Find the duplicated Xrefs... Some may not have an entrezgene, but multiples of other ids, so check
        if x['query'] in dup_xref and not pd.isnull(uniprot_map[x['query']]):
            # take the smallest numbered xref....
            if int(x['entrezgene']) < int(uniprot_map[x['query']]):
                uniprot_map[x['query']] = x['entrezgene']
    return uniprot_map

def parse_drugbank_target(res_entry):

    cid = get_dotval(res_entry, 'query')
    db_id = get_dotval(res_entry, 'drugbank.id')
    db_targets = get_dotval(res_entry, 'drugbank.targets')

    if pd.isnull(db_id) or type(db_targets) == float:
        return [], []

    if type(db_targets) != list:
        db_targets = [db_targets]

    parsed_target_nodes = []
    parsed_target_edges = []

    node_fields = ['gene_name', 'name', 'organism', 'source', 'uniprot']
    edge_fields = ['gene_name', 'name', 'actions', 'pmids', 'uniprot']

    for t in db_targets:
        n_vals = dict()
        e_vals = dict()
        for f in node_fields:
            n_vals[f] = get_dotval(t, f, list_ok=False)

        e_vals['drugbank_id'] = db_id
        e_vals['compound_q_id'] = cid
        for f in edge_fields:
            e_vals[f] = get_dotval(t, f, list_ok=False)

        parsed_target_nodes.append(n_vals)
        parsed_target_edges.append(e_vals)

    return parsed_target_nodes, parsed_target_edges


def get_drugbank_targets(mychem_dump):
    nodes = list()
    edges = list()

    for i, res_entry in enumerate(mychem_dump):
        try:
            n, e = parse_drugbank_target(res_entry)
        except:
            print(i)
            return None
        nodes += n
        edges += e

    return pd.DataFrame(nodes), pd.DataFrame(edges)

def parse_drugcentral_target(res_entry):
    cid = get_dotval(res_entry, 'query')
    dc_targets = get_dotval(res_entry, 'drugcentral.bioactivity')

    if pd.isnull(cid) or type(dc_targets) == float:
        return [], []

    if type(dc_targets) != list:
        dc_targets = [dc_targets]

    parsed_target_nodes = []
    parsed_target_edges = []

    node_fields = ['uniprot.gene_symbol', 'uniprot.swissprot_entry', 'target_name', 'target_class',
                   'organism', 'uniprot.uniprot_id']
    edge_fields = ['uniprot.gene_symbol', 'action_type', 'moa_source', 'moa', 'act_source', 'uniprot.uniprot_id',
                  'target_name']

    for t in dc_targets:
        n_vals = dict()
        e_vals = dict()
        for f in node_fields:
            n_vals[f.split('.')[-1]] = get_dotval(t, f, list_ok=False)

        e_vals['compound_q_id'] = cid
        for f in edge_fields:
            e_vals[f.split('.')[-1]] = get_dotval(t, f, list_ok=False)

        parsed_target_nodes.append(n_vals)
        parsed_target_edges.append(e_vals)

    return parsed_target_nodes, parsed_target_edges


def get_drugcentral_targets(mychem_dump):
    nodes = list()
    edges = list()

    for i, res_entry in enumerate(mychem_dump):
        try:
            n, e = parse_drugcentral_target(res_entry)
        except:
            print(i)
            return None
        nodes += n
        edges += e

    return pd.DataFrame(nodes), pd.DataFrame(edges)


def process_chemicals_in_mychem_results(results):

    parsed_mychem = []
    for r in results:
        out = get_xrefs(r)
        out['name'] = get_best_name(r)
        parsed_mychem.append(out)

    chem_node_df = pd.DataFrame(parsed_mychem)
    return chem_node_df
