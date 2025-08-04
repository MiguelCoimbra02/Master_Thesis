import itertools
import networkx as nx
import math
import ast
import matplotlib.pyplot as plt
import copy

'''def remove_edges_for_display(g):
    remove_edge = []
    checked_edges = set()
    for s, t in g.edges():
        if (t, s) not in checked_edges:
            checked_edges.add((s, t))
        else:
            remove_edge.append((s, t))
    g.remove_edges_from(remove_edge)'''

#def used to transform the tf_rank into tuples

def create_graph(nodes_data, edges_data):
    '''
    Creates a NetworkX graph from the given nodes and edges data from neo4j.
    '''
    g = nx.MultiGraph()
    for node_name, node_data in nodes_data.items():
        if not g.has_node(node_name):  # Check if the node already exists
            g.add_node(node_name, **node_data)

    for edge in edges_data:
        source, target, data = edge
        data_dict = {key: value for key, value in data}

        # Check if an edge already exists between source and target
        if not g.has_edge(source, target):
            g.add_edge(source, target, **data_dict)
    
    
    subgraph = autocomplete_graph_data(g)
    # ---- complete edge data with direct edges from the co-exp edges ----
    return subgraph

def autocomplete_graph_data(g, edges_to_check=None, min_co_exp_width = 1, max_co_exp_width = 10): 
    directed_edges_to_add = []
    edges_to_remove = []

    #print('EDGES_TO_CHECK: ',edges_to_check)

    # edges to check are important to avoid checking all edges when expanding
    edges_to_check = {(edge[0], edge[1]) for edge in edges_to_check} if edges_to_check else None
    
    for source, target, data in g.edges(data=True):
        # If edges_to_check is provided, only process those specific edges
        if edges_to_check and (source, target) not in edges_to_check and (target, source) not in edges_to_check:
            continue  # Skip edges that are not in the provided list
        
        #i Need this bc when expanding some edges will already be completed while others wont. if i run ast.literal_eval and key is already int gives error
        if isinstance(data['directed'], str):
            key = ast.literal_eval(data['directed'])
        else:
            key = data['directed']

        if isinstance(data['tf_rank'], str):
            tf_rank = ast.literal_eval(data['tf_rank'])
        else:
            tf_rank = data['tf_rank']
            

        if key == 1:
            data['tf_rank'] = tf_rank[1:] 
            directed_edges_to_add.append((source, target, key, dict(data))) #check later if here should be dict(data)
            data['tf_rank'] = [0] 
            data['directed'] = 0
    
    for source, target, key, data in directed_edges_to_add:
        g.add_edge(source, target, key, **data)
        g[source][target][key]['interaction'] = 'regulates expression'
        print('source', source, '-> target', target)


        source_is_tf = 'isTF' in g.nodes[source] and g.nodes[source]['isTF'] == '1'
        target_is_tf = 'isTF' in g.nodes[target] and g.nodes[target]['isTF'] == '1'

        # had to make undirect edges in DB, i saved source and target as attributes
        # Case 1: source regulates target
        if str(g[source][target][key]['source']) == str(source) and source_is_tf:

            g[source][target][key]['id'] = f'{source} modulates the expression of {target}'
            g[source][target][key]['arrows'] = {'to': {'enabled': True, 'scaleFactor': 0.5}}
            g[source][target][key]['source'] = source  # Add source metadata
            g[source][target][key]['target'] = target  # Add target metadata
            print('aliiiii')
        
        # Case 2: target regulates source
        elif str(g[source][target][key]['source']) == str(target) and target_is_tf:
            g[source][target][key]['id'] = f'{target} modulates the expression of {source}'
            g[source][target][key]['arrows'] = {'to': {'enabled': True, 'scaleFactor': 0.5}}  
            g[source][target][key]['source'] = target  # Add source metadata
            g[source][target][key]['target'] = source  # Add target metadata
            print('aquiiiiii')
        
        else:
            edges_to_remove.append((source, target, key))
    #removing edges
    print('edges_to_remove', edges_to_remove)
    
    for source, target, key in edges_to_remove:
        g.remove_edge(source, target, key=key)

    #dict tf_rank : edge width
    #when converting a list of len 1 into tuple ill get (x,)
    tf_width_dict = {
        (1,): 7,
        (2,): 7 ,
        (1,2): 9,
        (3,): 11,
        (1,3): 13,
        (2,3): 13,
        (1,2,3): 15, 
    }

    #defining width
    for fr, to, key, attrs in g.edges(keys=True, data=True):
        #co_exp edges
        if key == 0:
            #scaling the width by the irp_score in a defined interval
            width = min_co_exp_width + (float(attrs['irp_score']) * (max_co_exp_width - min_co_exp_width))
            attrs['width'] = width
        
        #directed edges
        else:
            #width of directed edges by the dict
            tf_rank = tuple(attrs['tf_rank'])
            if tf_rank in tf_width_dict:
                attrs['width'] = tf_width_dict[tf_rank]
    return g


def filter_edges_nodes(graph, query_nodes, tf_ranks, rangeSlider_co_exp):
    g = graph.copy()
    to_remove_e1 = []
    to_remove_n1 = set()  # Set to hold nodes to be removed
    if tf_ranks == [] and rangeSlider_co_exp == 0:
        return g
        
    if len(tf_ranks) == 0:
        pass
    else:  
        for s, t, k, d in g.edges(keys=True, data=True):  
            tf_ranks_copy = tf_ranks.copy()
            print('rf rank of edge: ', d['tf_rank'])
            if k == 0:  
                if 0 not in tf_ranks_copy:
                    to_remove_e1.append((s, t, k))
            else:
                 #garantee 1,2,3 not in tf_ranks
                 
                 if 4 in tf_ranks_copy and all(rank not in tf_ranks_copy for rank in [1, 2, 3]):
                     continue
                 else:
                    if 4 in tf_ranks_copy:
                         tf_ranks_copy.remove(4)
                    if 0 in tf_ranks_copy:
                         tf_ranks_copy.remove(0)

                    if not tf_ranks_copy:
                        to_remove_e1.append((s, t, k))
                        continue

                    #RETOMAR AQUI
                    
                    tf_rank_list = ast.literal_eval(d['tf_rank']) if isinstance(d['tf_rank'], str) else d['tf_rank']
                    print('rf rank of edge: ',tf_rank_list)
                    if not all(rank in tf_rank_list for rank in tf_ranks_copy):  
                        to_remove_e1.append((s, t, k))

        print('tf_rank of the search',tf_ranks_copy)
        g.remove_edges_from(set(to_remove_e1))
        
        for node in g.nodes():
            # Skip query nodes
            if node in query_nodes:
                continue
            
            # Check if there is a path to any query node
            is_connected = any(nx.has_path(g, query_node, node) for query_node in query_nodes)
            if not is_connected:
                to_remove_n1.add(node)
        g.remove_nodes_from(to_remove_n1)
   
    #range Slider only for undirected edges
    to_remove_e = []
    for s, t, k, d in g.edges(keys=True, data=True):
        if k == 0:
            if ast.literal_eval(d["irp_score"]) < rangeSlider_co_exp:
                to_remove_e.append((s, t, k))

    g.remove_edges_from(set(to_remove_e))


    to_remove_n = set()  # Set to hold nodes to be removed

    for node in g.nodes():
        # Skip query nodes
        if node in query_nodes:
            continue
        
        # Check if there is a path to any query node
        is_connected = any(nx.has_path(g, query_node, node) for query_node in query_nodes)
        if not is_connected:
            to_remove_n.add(node)
    g.remove_nodes_from(to_remove_n)   

    # print(f"Final nodes: {list(g.nodes)}")
    # print(f"Final edges: {list(g.edges(keys=True))}") 
    return g


def get_top_edges_by_irpscore(graph, top_n = 100):
    g = nx.MultiGraph(graph.copy())
    print('number_viz_nodes', top_n)
    print(g.number_of_nodes())
    if g.number_of_nodes() <= top_n:
        return g
    
    all_edges = list(g.edges(keys=True, data=True))

    #add direct edges
    nodes_2_keep = set()
    for u, v, k, data in all_edges:
        if k == 1 and len(nodes_2_keep) < top_n:
            nodes_2_keep.add(u)
            nodes_2_keep.add(v)
    
    irp_score_sorted = sorted(all_edges, key=lambda x: x[3].get('irp_score', 0), reverse=True)
    
    for u, v, k, data in irp_score_sorted:
        nodes_2_keep.add(u)
        nodes_2_keep.add(v)

        if len(nodes_2_keep) > top_n:
            break
    for node in list(g.nodes):
        if node not in nodes_2_keep:  
            g.remove_node(node)

    return g
  


def expanded_graph4display(g_b4_expand, g_after_expand, top_n=50):
    g_b4_expand = g_b4_expand.copy()
    g_after_expand = g_after_expand.copy()

    nodes_g_b4_expand = set(g_b4_expand.nodes())
    nodes_g_after_expand = set(g_after_expand.nodes())

    # Extract edges as (u, v, k, irp_score)
    edges_g_b4_expand = set((u, v, k, d.get('irp_score', 0)) for u, v, k, d in g_b4_expand.edges(keys=True, data=True))
    edges_g_after_expand = set((u, v, k, d.get('irp_score', 0)) for u, v, k, d in g_after_expand.edges(keys=True, data=True))

    # Correctly initialize sets!
    new_edges = set()
    for edge_after in edges_g_after_expand:
        if edge_after not in edges_g_b4_expand:
            new_edges.add(edge_after)

    new_nodes = set()
    for node_after in nodes_g_after_expand:
        if node_after not in nodes_g_b4_expand:
            new_nodes.add(node_after)

    print('new_edges', new_edges)
    print('len_new edges', len(new_edges))
    print('len_new nodes', len(new_nodes))

    if len(new_nodes) <= top_n:
        return g_after_expand
    else:
        nodes_2_keep = set()

        # Include nodes connected by edges with k == 1
        for u, v, k, irp_score in new_edges:
            if k == 1:
                nodes_2_keep.add(u)
                nodes_2_keep.add(v)

        # Sort by irp_score
        irp_score_sorted = sorted(new_edges, key=lambda edge: edge[3], reverse=True)
        print('irp_score_sorted', irp_score_sorted)
        print('len irp_score_sorted', len(irp_score_sorted))

        # Add top scoring nodes
        for u, v, k, irp_score in irp_score_sorted:
            if len(nodes_2_keep) >= top_n:
                break
            nodes_2_keep.add(u)
            nodes_2_keep.add(v)

    # Remove all new nodes not in the final kept set
    for node in new_nodes:
        if node not in nodes_2_keep:
            g_after_expand.remove_node(node)

    print('len nodes2keep', len(nodes_2_keep))
    return g_after_expand


def graph2json(g, query_nodes=[]):
    nlist = []
    for nodeid, attrs in g.nodes(data=True):
        nodeData = copy.deepcopy(attrs)
        nodeData['id'] = nodeid
        nodeData['label'] = nodeid
        
        if nodeData['isTF'] == '1': #in case isTF
            nodeData['color'] = {'border': '#436e67',
                                 'background': '#6c9585'}
            nodeData['shape'] = 'dot'
        
        elif nodeData['isTR'] == '1': #in case isTR
            nodeData['color'] = {'border': '#abbb95 ',
                                 'background': '#bed0a6'}
            nodeData['shape'] = 'dot'
                
        for key, value in nodeData.items():
            if isinstance(value, float) and math.isnan(value):
                nodeData[key] = None  # Use None instead of null for JSON serialization

        if nodeid in query_nodes:
            nodeData['color'] = {'border': "#233a48",
                                 'background': '#ffa500'
                                } 
            nodeData['borderWidth'] = 2
        nlist.append(nodeData)

    elist = []            
    for fr, to, attrs in g.edges(data=True): 
        if ast.literal_eval(attrs['cis_elements']) == 0.0:
            attrs['cis_elements'] = '--'
        elist.append({
                    'from': attrs.get('source', fr),
                    'to': attrs.get('target', to),
                    'id': attrs.get('id', ''),
                    'label': attrs.get('interaction', ''),
                    'irp_score': attrs.get('irp_score', 0),
                    'ConnecTF': attrs.get('connecTF', '--'),
                    'edge_type': attrs.get('edge_type', '--'),
                    'gene_name': attrs.get('gene_name', '--'),
                    'cis_elements': attrs.get('cis_elements', 0),
                    'cis_value': attrs.get('cis_value', 0.0),
                    'dap_seq': attrs.get('dap_seq', 0),
                    'tf_rank': attrs.get('tf_rank', [0]),  
                    'width': attrs.get('width', 1),  
                    'directed': attrs.get('directed', 'no'),
                    'arrows': attrs.get('arrows', 'undefined'),
                    'hidden': attrs.get('hidden', False)
                })
    
        
    result =  {'network': {'nodes': nlist, 'edges': elist}}
    return result
