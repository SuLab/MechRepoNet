import json
import base64
import datetime
import pickle
import time
from collections import defaultdict
from itertools import chain
import requests
from tqdm import tqdm
from itertools import cycle
#from processing import load_api_results

# url = "https://stitcher.ncats.io/api/stitches/v1/9ZOQ3TZI87"
# d = requests.get(url).json()
xref_keys = ['CAS', 'Cas', 'ChEMBL', 'ChemblId', 'CompoundCID', 'CompoundName', 'CompoundSmiles',
             'CompoundUNII', 'DRUG BANK', 'DrugbankId', 'IUPHAR', 'InChIKey', 'Iupac',
             'KeggId', 'MESH', 'NCI_THESAURUS', 'NDF-RT', 'PUBCHEM', 'RXCUI', 'SMILES', 'UNII', 'Unii',
             'drugbank-id', 'smiles', 'unii']
combine_keys = [
    ('UNII', 'unii', 'Unii', 'CompoundUNII'),
    ('CAS', 'Cas'),
    ('DRUGBANK', 'DrugbankId', 'DRUG BANK', 'drugbank-id'),
    ('CHEMBL.COMPOUND', 'ChEMBL', 'ChemblId'),
    ('SMILES', 'smiles', 'CompoundSmiles')
]


def download_stitcher():
    url = "https://stitcher.ncats.io/api/stitches/v1?top=100&skip={}"
    skip = 0
    contents = []
    t = tqdm(total=102012 / 100)
    while True:
        t.update()
        d = requests.get(url.format(skip)).json()
        if not d['contents']:
            break
        contents.extend(d['contents'])
        skip += 100
        time.sleep(1)

    return contents


def alwayslist(value):
    """If input value if not a list/tuple type, return it as a single value list."""
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return value
    else:
        return [value]


def organize_one_record(d, rec_num=0, compressed=False):
    node_xref = defaultdict(lambda: defaultdict(set))
    for xref_key in xref_keys:
        for node in alwayslist(d['sgroup']['properties'].get(xref_key, [])):
            node_xref[node['node']][xref_key].add(node['value'])
    node_xref = {k: dict(v) for k, v in node_xref.items()}

    # combine these sets of keys together
    # the first will be the prefix for curieutil and is what they get combined into

    for xrefs_kv in node_xref.values():
        for combine_key in combine_keys:
            xrefs_kv[combine_key[0]] = set(chain(*[xrefs_kv.get(k, set()) for k in combine_key]))
            for k in combine_key[1:]:
                if k in xrefs_kv:
                    del xrefs_kv[k]

    # Drug Products
    targets_nv = alwayslist(d['sgroup']['properties'].get('Targets', []))
    conditions_nv = alwayslist(d['sgroup']['properties'].get('Conditions', []))
    for t in targets_nv:
        t['actualvalue'] = json.loads(base64.b64decode(t['value']).decode())
    for t in conditions_nv:
        t['actualvalue'] = json.loads(base64.b64decode(t['value']).decode())

    conditions = defaultdict(list)
    for condition in conditions_nv:
        conditions[condition['node']].append(condition['actualvalue'])
    conditions = dict(conditions)

    targets = defaultdict(list)
    for target in targets_nv:
        targets[target['node']].append(target['actualvalue'])
    targets = dict(targets)

    records = []
    for k in set(conditions.keys()) | set(targets.keys()) | set(node_xref.keys()):
        records.append({'rec_num': rec_num, 'node_num': k, 'xrefs': node_xref.get(k, dict()),
                        'conditions': conditions.get(k, []), 'targets': targets.get(k, [])})

    if compressed:
        return compress_record(records)

    return records


def compress_record(d):
    nodes = [n['node_num'] for n in d]
    x_refs = dict()
    for n in d:
        for k, v in n['xrefs'].items():
            x_refs[k] = x_refs.get(k, set()).union(v)
    conditions = []
    targets = []
    for n in d:
        conditions += n['conditions']
        targets += n['targets']

    return {'nodes': nodes, 'xrefs': x_refs, 'conditions': conditions, 'targets': targets}


def unknown_to_none(items):
    return [None if i == 'Unknown' else i for i in items]


def extract_conditions(d):
    meshids = []
    phases = []
    modalities = []
    doids = []
    do_precision = []
    mesh_precision = []
    names = []
    dates = []

    for condition in d['conditions']:
        meshids.append(condition.get('ConditionMeshId'))
        phases.append(condition.get('HighestPhase'))
        modalities.append(condition.get('TreatmentModality'))
        doids.append(condition.get('ConditionDoId'))
        do_precision.append(condition.get('isConditionDoImprecise'))
        mesh_precision.append(condition.get('isConditionMeshImprecise'))
        names.append(condition.get('ConditionName'))
        dates.append(condition.get('ConditionProductDate'))

    # Unknown and None should be Equivalent
    doids = unknown_to_none(doids)
    meshids = unknown_to_none(meshids)

    # Add CURI to disease identifiers
    doids = ['DOID:'+str(i) if i is not None else i for i in doids]
    meshids = ['MESH:'+i if i is not None else i for i in meshids]

    out = {'meshids': meshids, 'phases': phases, 'modalities': modalities, 'doids': doids,
                'do_precision': do_precision, 'mesh_precision': mesh_precision, 'names': names, 'dates': dates}

    # Unknown and None should be Equivalent
    return {k: unknown_to_none(v) for k, v in out.items()}


def expand_data(data, length):
    if type(data) == str:
        return [data]*length
    if type(data) == list:
        data_len = len(data)
        if data_len == length:
            return data
        elif data_len < length:
            return [x for i, x in zip(range(length), cycle(data))]
        else:
            return data[:length]
    else:
        return [data]*length


def extract_targets(d):
    t_ids = []
    id_sources = []
    pharmas = []
    names = []

    # IF WE GET A LIST, ENSURE EACH OF SAME LENGTH
    for t in d['targets']:
        if type(t['PrimaryTargetId']) == str:
            t_ids.append(t.get('PrimaryTargetId'))
            id_sources.append(t.get('PrimaryTargetType'))
            pharmas.append(t.get('TargetPharmacology'))
            names.append(t.get('PrimaryTargetLabel'))

        elif type(t.get('PrimaryTargetId')) == list:
            l = len(t.get('PrimaryTargetId'))
            t_ids += t.get('PrimaryTargetId')
            id_sources += expand_data(t.get('PrimaryTargetType'), l)
            pharmas += expand_data(t.get('TargetPharmacology'), l)
            names += expand_data(t.get('PrimaryTargetLabel'), l)

    out = {'ids': t_ids, 'sources': id_sources, 'pharmas': pharmas, 'names': names}
    # Unknown and None should be Equivalent
    return {k: unknown_to_none(v) for k, v in out.items()}


def organize_data(contents):
    records = []
    for i, d in tqdm(enumerate(contents), total=len(contents)):
        cur_rec = organize_one_record(d, i, compressed=True)
        cur_rec['conditions'] = extract_conditions(cur_rec)
        cur_rec['targets'] = extract_targets(cur_rec)

        records.append(cur_rec)

    return records


def load_parsed_data():
    d = load_api_results("stitcher_parsed_{}.pkl")
    if d is not None:
        return d
    else:
        raise FileNotFoundError('Parsed Sticher file not found: Please run main() to get this file.')


def main():
    contents = load_api_results('stitcher_dump_{}.pkl', False, download_stitcher)
    d = organize_data(contents)
    with open("stitcher_parsed_{}.pkl".format(datetime.datetime.now().strftime("%Y-%m-%d")), "wb") as f:
        pickle.dump(d, f)


if __name__ == '__main__':
    main()
