import os
import json
from collections import Counter
from wikidataintegrator.wdi_core import WDItemEngine
import sys
sys.path.append('../hetnet-ml/src')
import graph_tools as gt


def parse_result_uris(result):
    for c in result:
        idx = result[c].str.startswith('http://www.wikidata.org/entity')
        if sum(idx) > 0:
            result[c][idx] = result[c][idx].str.split('/', expand=True).iloc[:,-1]
    return result.drop_duplicates()


def get_col_order(query_text):
    lines = iter(query_text.split('\n'))
    line = next(lines)
    while not line.lower().lstrip().startswith('select'):
        line = next(lines)
    return [l.strip(' ') for l in line.split(' ?')[1:]]


def execute_sparql_query(query_text, endpoint='http://avalanche.scripps.edu:9999/bigdata/sparql'):
    # Execute the qurey
    result = WDItemEngine.execute_sparql_query(query_text, endpoint=endpoint, as_dataframe=True)
    # Don't do any processing if empy result
    if len(result) == 0:
        return result
    # Enforce the proper column order and return
    col_order = get_col_order(query_text)
    return parse_result_uris(result)[col_order]


def generate_instance_tag(name, p, n_id):
    return "\n    {n} {p} wd:{n_id} .".format(n=name, p=p, n_id=n_id)


def to_query_name(name):
    return '?' + name.lower().replace(' ', '_').replace('-', '')


def append_node_counts(query_names):
    counter = Counter(query_names)
    for name, count in counter.items():
        if count > 1:
            for i in range(1, count+1):
                idx = query_names.index(name)
                query_names[idx] = query_names[idx] + str(i)
    return query_names


def get_query_numb(qname):
    try:
        numb = str(int(qname[-1]))
    except ValueError:
        numb = ''
    return numb


def build_header(query_names):
    header_text = "SELECT DISTINCT "
    for name in query_names:
        header_text += '{n} {n}Label '.format(n=name)
    return header_text


def where_tag():
    return "\nWHERE {"


def footer():
    return '\n\n    SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGAGE],en" }\n}'


def build_node(name, qname, n_id, p, qualifiers=iter(())):
    """

    :param name:
    :param qname:
    :param n_id:
    :param p:
    :param qualifiers: list of dicts, keys are 'node_q_names', a list of query_names for the qualifier (e.g. ['?gene1'])
        'node_q_types' list of typing information for node in the qualifying statement,(e.g. ['wd:Q7187']),
        'p', the predicate in the qualifier and 'not', if True, will filter relationships that meet this condition,
        if False, will filter relationships that do not meet this conidtion.
        Qualifying statements that further define the node type. See WDQueryBuilder._parse_qualifier() for more
        example of one qualifier:   {'names': {'s': '?gene1', 'o': '?protein1'},
                                    'types': {{'?gene1': {'id': 'Q7187', 'p': 'wdt:P31'}}},
                                    'p': 'P688',
                                    'not': False}
    :return:
    """
    numb = get_query_numb(qname)

    # Build the appropriate info to get the correct subject and object node-types
    # Use set so not repeated in cases of self-referential edges
    node_text = '\n\n    # Initial typing for {name}{numb}'.format(name=name, numb=numb)
    node_text += generate_instance_tag(qname, p, n_id)

    visited = {qname}

    # Get any qualifiers for subjects and objects
    for qual in qualifiers:
            node_text += '\n    # Qualifier for {name}{numb}'.format(name=name, numb=numb)
            node_text += build_qualifier(visited, **qual)
            visited.update(set(qual['types'].keys()))

    return node_text


def build_qualifier(visited, names, types, p, q_not):
    qual_out = ''

    for qname, vals in types.items():
        if qname not in visited:
            qual_out += generate_instance_tag(qname, vals['p'], vals['id'])

    # Negative qualifiers need to filter out edges
    if q_not:
        qual_out += '\n    FILTER NOT EXISTS {{ {s} wdt:{p} {o} }} .'.format(s=names['s'],
                                                                             p=p, o=names['o'])
    # Otherwise add the edges
    else:
        qual_out += '\n    {s} wdt:{p} {o} .'.format(s=names['s'],
                                                     p=p, o=names['o'])
    return qual_out


def get_edge_statement_text(s_qname, p_id, p_x, o_qname, qual=None):
    if not qual:
        out_text = "{s} wdt:{p}{p_x} {o}".format(s=s_qname, p=p_id, p_x=p_x, o=o_qname)
    else:
        out_text = '{s} p:{p}{p_x} [ps:{p}{p_x} {o};'.format(s=s_qname, p=p_id, p_x=p_x, o=o_qname)
        idx = out_text.index('[')
        shift = len(out_text[:idx+1])
        out_text += '\n' + ' '*shift + 'pq:{qual} ?qualifier]'.format(qual=qual)
    return out_text


def build_single_edge(s_qname, o_qname, sp_x, op_x, id, qual=None, rev_id=None, rev_qual=None, abbrev=None):

    # Some edges are bi-directonal e.g. drug-used-for-treatment and medical-condition-treated
    # Get edges in eitehr dircetion
    if rev_id is not None:
        edge_text = "\n\n    {{ {fwd} }}".format(fwd=get_edge_statement_text(s_qname, id, op_x, o_qname, qual))
        edge_text += "\n    UNION {{ {rev} }}".format(rev=get_edge_statement_text(o_qname, rev_id, sp_x, s_qname, rev_qual))
    else:
        edge_text = "\n\n    " + get_edge_statement_text(s_qname, id, op_x, o_qname, qual)

    return edge_text


def build_multi_edge(node_q_names, ps, p_extends):
    """No qualifier support for now"""

    edge_text = '\n'
    for i, (p, p_x) in enumerate(zip(ps, p_extends)):
        edge_text += "\n    " + get_edge_statement_text(node_q_names[i], p, p_x, node_q_names[i+1]) + '.'
    return edge_text


class WDHetnetQueryBuilder:

    def __init__(self, node_info, edge_info):

        self.node_info = self._get_info(node_info)
        self.edge_info = self._get_info(edge_info)

        self._validate_nodes(self.node_info)
        self._validate_edges(self.edge_info)

        self.node_abbrevs = self.generate_abbrevs_from_info(self.node_info)
        self.edge_abbrevs = self.generate_abbrevs_from_info(self.edge_info)

        self.node_abv_to_full = self.abbrev_to_full(self.node_abbrevs)
        self.edge_abv_to_full = self.abbrev_to_full(self.edge_abbrevs)

        self.subclass = self._extract_node_key('subclass')
        self.extend = self._extract_node_key('extend')

    def _get_info(self, info):
        # if fileame for json file is given, read the file.
        if os.path.isfile(info):
            info = json.load(open(info, 'r'))

        return info

    def _validate_nodes(self, info):
        self._validate_info(info)
        for v in info.values():
            assert 'subclass' in v
            assert 'extend' in v

        # ids for nodes can either be Q-values or None
        id_first_char = set(v['id'][0] for v in info.values() if v['id'] is not None)
        assert len(id_first_char) == 1
        id_first_char = id_first_char.pop()
        assert id_first_char == 'Q'

    def _validate_edges(self, info):
        self._validate_info(info)

        # ids for edges must begin with a P
        id_first_char = [v['id'][0] for v in info.values() if type(v['id']) == str]

        # Number of inner types == number of edges - 1  e.g.:
        # (start_type) edge_1 (inner_type1) edge_2 (inner_type2) edge_3 (end_type)
        for v in info.values():
            if type(v['id']) == list:
                assert len(v['id']) == len(v['inner_types'])+1
                # Grab the first characters to make sure they conform
                for i in v['id']:
                    id_first_char.append(i[0])

        # all list ids should start with P
        id_first_char = set(id_first_char)
        assert len(id_first_char) == 1
        id_first_char = id_first_char.pop()
        assert id_first_char == 'P'

    def _validate_info(self, info):
        # Make sure its a dictionary
        assert type(info) == dict

        # Ensure all values have an id and an abbreviation
        for v in info.values():
            assert 'id' in v
            assert 'abbrev' in v

        # Ensure unique abbrevations
        all_abbrevs = set(v['abbrev'] for v in info.values())
        assert len(all_abbrevs) == len(info)

    @staticmethod
    def generate_abbrevs_from_info(info):
        return {k: v['abbrev'] for k, v in info.items()}

    @staticmethod
    def abbrev_to_full(abbrevs):
        return {v: k for k, v in abbrevs.items()}

    def _extract_node_key(self, key):
        return {v['id']: v[key] for v in self.node_info.values() if v['id'] is not None}

    @staticmethod
    def determine_p(use_subclass, extend=True):
        p = "wdt:P279*" if use_subclass else "wdt:P31/wdt:P279*"
        # Option to not extend down 'subclass_of' edges (useful for highly populated node types)
        if not extend:
            p = p.replace('/wdt:P279*', '').replace('*', '')
        return p

    def parse_qualifiers(self, qualifiers, numb=''):
        parsed = list()
        for qual in qualifiers:
            parsed.append(self._parse_qualifier(qual, numb))
        return parsed

    def _parse_qualifier(self, qual, qual_num=''):
        node_q_names = {}

        for node in ['s', 'o']:
            node_name = qual[node]

            if node_name in self.node_info.keys():
                node_q_names[node] = to_query_name(node_name) + qual_num
            elif node_name is None:
                node_q_names[node] = '?na' + qual_num
            else:
                node_q_names[node] = 'wd:' + node_name

        node_q_types = dict()
        for semtype, qname in node_q_names.items():
            if not qname.startswith('wd:') and qname != '?na' + qual_num:
                current_type = {}
                qid = self.node_info[qual[semtype]]['id']

                current_type['id'] = qid
                current_type['p'] = self.determine_p(self.subclass.get(qid, False), self.extend.get(qid, False))
                node_q_types[qname] = current_type

        return {'names':node_q_names, 'types':node_q_types, 'p': qual['p'], 'q_not': qual.get('not', False)}

    def build_query_from_abbrev(self, abbrev, target=None):
        # Split the edge abbreviation in to subject, predicate, object
        s_abv, p_abv, o_abv = gt.parse_edge_abbrev(abbrev)

        # Get the full versions of names
        s_name = self.node_abv_to_full[s_abv]
        p_name = self.edge_abv_to_full[p_abv]
        o_name = self.node_abv_to_full[o_abv]

        # Get all the node types from node and edge info
        nodes_in_query = [s_name]
        nodes_in_query += [self.node_abv_to_full[n] for n in self.edge_info[p_name].get('inner_types', [])]
        nodes_in_query += ['qualifier'] if 'qual' in self.edge_info[p_name] else []
        nodes_in_query += [o_name]

        # Convert to a queryname (e.g. ?compound1)
        query_names = [to_query_name(n) for n in nodes_in_query]
        if len(set(query_names)) < len(query_names):
            query_names = append_node_counts(query_names)

        query_text = build_header(query_names)
        query_text += where_tag()

        # Build definitions for nodes to capture
        for qname, nname in zip(query_names, nodes_in_query):
            # Typing for qualifiers comes later
            if nname == 'qualifier':
                continue
            numb = get_query_numb(qname)
            nid = self.node_info[nname]['id']
            if nid is None:
                continue
            p = self.determine_p(self.subclass.get(nid, False), self.extend.get(nid, False))
            quals = self.parse_qualifiers(self.node_info[nname].get('qualifiers', iter(())), numb)
            query_text += build_node(nname, qname, nid, p, quals)

        # Add in the specific edge we're after
        p_info = self.edge_info[p_name]
        s_qname = query_names[0]
        o_qname = to_query_name(target) if target is not None else query_names[-1]
        sp_x = self.node_info[s_name].get('p_extend', '')
        op_x = self.node_info[o_name].get('p_extend', '')

        if type(p_info['id']) == list:
            p_extends = [self.node_info[nname].get('p_extend', '') for nname in nodes_in_query]
            query_text += build_multi_edge(query_names, p_info['id'], p_extends)
        else:
            query_text += build_single_edge(s_qname, o_qname, sp_x, op_x, **p_info)
        query_text += footer()

        return query_text
