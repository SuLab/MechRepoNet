import shlex
import pandas as pd
from collections import defaultdict


def read_lines(filename):
    """Generator to read lines of a file"""
    with open(filename, 'r') as fin:
        for line in fin.readlines():
            yield line


def get_term_keys(obo_file):
    """
    Returns all the keys in terms in the .obo file.

    :param obo_file: The location of the obo file to be examined.
    :return: set, the keys for all terms within the obo file.
    """
    term_keys = set()
    is_term = False

    for line in read_lines(obo_file):
        # Determine if term or a typedef
        if '[Term]' in line:
            is_term = True
            continue
        if line == '\n':
            is_term = False

        if is_term:
            term_keys.add(line.split(':')[0])
    return term_keys


def get_prop_types(obo_file):
    """Determines all the `property_value` types within a .obo file.

    E.G. for the following line:
        property_value: http://purl.obolibrary.org/obo/chebi/inchikey "SBLSYFIUPXRQRY-UHFFFAOYSA-N" xsd:string
    The property_value wouls be http://purl.obolibrary.org/obo/chebi/inchikey

    :param obo_file: the locaiton of the obo file to examine.
    :return: set, with all the property_value predicates.
    """
    props = set()
    for line in read_lines(obo_file):
        if line.startswith('property_value'):
            props.add(parse_prop_val(line)[0])
    return props


def parse_single_value(line, key):
    """When the key is already known, returns the value as parsed from the line"""
    # Remove key from the line
    val = line.lstrip(key+':').strip()
    # Convert boolean text to values
    if val == 'true':
        val = True
    elif val == 'false':
        val = False
    return val


def parse_prop_val(line):
    """
    Parses a property_value line in an ontology. Returns property type and the value

    Examples of lines parsed:
        property_value: http://purl.obolibrary.org/obo/chebi/inchikey "SBLSYFIUPXRQRY-UHFFFAOYSA-N" xsd:string
        property_value: IAO:0000589 "cell part (CARO)" xsd:string
        property_value: IAO:0000412 http://purl.obolibrary.org/obo/caro.owl
    """

    # Shlex is too slow to run on every property_value, so only do as last resort...
    spl = line.split(' ')
    if len(spl) != 4:
        # splits on spaces if not in "quotes like this"
        spl = shlex.split(line)

    prop_type = spl[1]

    # Most property values encapsulated in quotes, but not all:
    try:
        # Removes quotes
        prop_val = eval(spl[2])
    except:
        # Keeps as-is
        prop_val = spl[2]
    return prop_type, prop_val


def get_ontology_nodes(obo_file, prefix=None, props=None):
    """
    Converts an ontology's .obo file into a DataFrame of nodes.

    :param obo_file: Location of the obofile to be read.
    :param prefix: if not none, the prefix for terms to extract as nodes. Terms with a prefix
        not matching this value will be skipped.
    :param props: List, the `property_value` properties to extract for the nodes. See `get_prop_types` to
        get a list of avaliable properties to extract.

    :return: DataFrame with ids, names xrefs, and other relevatent data for terms in an ontology.
    """
    nodes = list()
    is_term = False
    single_value = ['id','def','name', 'is_obsolete', 'creation_date']
    multi_value = ['replaced_by', 'alt_id', 'xref', 'subset']


    for line in read_lines(obo_file):
        # Terms start with specific Flag...
        if '[Term]' in line:
            is_term = True
            term = dict()
            mv_items = defaultdict(list)
            continue
        # Empty lines indicate new terms
        if line == '\n':
            if is_term:
                nodes.append({**{k: '|'.join(v) for k, v in mv_items.items()}, **term})
            is_term = False
            continue

        # Get single values
        key = line.split(':')[0]

        if key in single_value:
            # If we only want one Ontology's worth of terms,
            # ensure that the correct Prefex is on the identifier
            if prefix and key == 'id':
                _id = parse_single_value(line, key)
                if not _id.startswith(prefix):
                    is_term = False
                    continue
                else:
                    term[key] = _id
            else:
                term[key] = parse_single_value(line, key)

        if key in multi_value:
            mv_items[key].append(parse_single_value(line, key))
        if props is not None and key == 'property_value':
            prop_type, prop_val = parse_prop_val(line)
            if prop_type in props:
                term[prop_type] = prop_val

    out = pd.DataFrame(nodes)
    if 'is_obsolete' in out:
        out['is_obsolete'] = out['is_obsolete'].fillna(False)

    if prefix is None:
        out['id_src'] = out['id'].str.split(':', expand=True)[0]

    return out


def parse_edge_line(line):
    """
    Parses lines thats represent an edge in an .obo file. Returns object ID, name(if available) and predicate.
    Currently only supports edges starting with 'is_a:' and 'relationship:'

    Examples of edge test-cases parsed:
      is_a: HP:0001392 ! Abnormality of the liver
      relationship: has_component CHEBI:16412 {cardinality="2"} ! lipopolysaccharide
      relationship: has_part MOD:00160 ! N4-glycosyl-L-asparagine
      relationship: is_conjugate_base_of CHEBI:17883
    """
    line_spl = line.split(' ')

    # Get the targets name
    if '!' in line:
        name = line[line.index('!')+1:].strip()
    else:
        name = float('nan')

    # is_a relationship: 'is_a: tgt_id ! tgt_name'
    if line_spl[0] == 'is_a:':
        target = line_spl[1].strip()
        rel_type = 'is_a'
    else:
        # relationship: 'relationship: rel_type tgt_id {possible_other_junk} ! tgt_name'
        target = line_spl[2].strip()
        rel_type = line_spl[1].strip()

    return {'tgt_id': target, 'rel_type': rel_type, 'tgt_name': name}


def get_ontology_edges(obo_file, prefix=None):
    """
    Produces a DataFrame of s, p, o triples for edges in a .obo file. (Names are also returned)

    :param obo_file", file name of the obofile to be read
    :param prefix: str, ontology prefix. Will only return edges where the subject's prefix matches.
        E.G. prefix='CHEBI:' will not return any edges for Terms starting with 'BFO:'
    :return: DataFrame, cols=['src_id', 'src_name', 'rel_type', 'tgt_id', 'tgt_name']
    """
    # Initialize edges and Term keys to extract
    edges = list()
    source_info = ['id', 'name']
    edge_lines = ['relationship', 'is_a']
    is_term = False

    for line in read_lines(obo_file):
        # Terms start with specific Flag...
        if '[Term]' in line:
            is_term = True
            term_info = dict()
            continue
        # Empty lines indicate new terms
        if line == '\n':
            is_term = False
            continue

        # Get single values
        key = line.split(':')[0]

        if key in source_info:
            # If we only want one Ontology's worth of terms,
            # ensure that the correct Prefex is on the identifier
            if prefix and key == 'id':
                _id = parse_single_value(line, key)
                if not _id.startswith(prefix):
                    is_term = False
                    continue
                else:
                    term_info['src_'+key] = _id
            else:
                term_info['src_'+key] = parse_single_value(line, key)

        # Extract teh edge specific info
        if key in edge_lines and is_term:
            edge_info = parse_edge_line(line)
            edges.append({**term_info, **edge_info})


    return pd.DataFrame(edges)[['src_id', 'src_name', 'rel_type', 'tgt_id', 'tgt_name']]
