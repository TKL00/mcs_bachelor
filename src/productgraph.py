import networkx as nx
import itertools

def _has_node_in_common(u, v, n_coordinates):
    for i in range(n_coordinates):
        if u[i] == v[i]:
            return True
    return False

### Authors: Tobias Klink Lehn (toleh20@student.sdu.dk) and Kasper Halkjær Beider (kbeid20@student.sdu.dk)
def product_graph(L):
    """
        Computes the modular product of a list of graphs.
    """

    product_graph = nx.Graph()
    ## Computes the cartesian products of the two node sets (sorted)
    node_list = [sorted(list(g.nodes)) for g in L]
    product_nodes = list(itertools.product(*node_list))
    n_graphs = len(L)

    product_graph.add_nodes_from(product_nodes)
    node_count = len(product_graph.nodes)

    product_graph_nodes = list(product_graph.nodes)

    for i in range(node_count):
        ## v = (v_1, v_2, ..., v_n)
        node_i = product_graph_nodes[i]

        neighbourhoods = []
        for L_i_node in range(n_graphs):
            ## Add the neighbourhood of v_i in product node from its corresponding graph L[i]
            neighbourhoods.append(L[L_i_node].adj[node_i[L_i_node]])

        for j in range(i + 1, node_count):
            ## u = (u_1, u_2, ..., u_n)
            node_j = product_graph_nodes[j]

            if not _has_node_in_common(node_i, node_j, n_graphs):

                all_agree_adj = True
                all_agree_not_adj = True 
                for index in range(n_graphs):
                    ## Determine node v_i's neighbourhood
                    v_i_neighbourhood = neighbourhoods[index]
                    
                    ## If node v_i is connected to u_i, then all nodes can't agree on being non-adjacent
                    if node_j[index] in v_i_neighbourhood:
                        all_agree_not_adj = False
                    ## Opposite argument
                    elif node_j[index] not in v_i_neighbourhood:
                        all_agree_adj = False
                    
                    ## Stop looking as no consensus was found
                    if not all_agree_adj and not all_agree_not_adj: break
                
                ## add edge between this node and the anchor point
                if all_agree_adj:
                    product_graph.add_edge( node_i, node_j, color="blue")

                elif all_agree_not_adj:
                    product_graph.add_edge( node_i, node_j, color="red")

    return product_graph 

def product_graph_list(L, lg_node_anchor={}):
    """
        Computes the modular product of all NetworkX graphs contained in L.
        lg_node_anchor maps a node, u, of L[0] to all nodes, v_i, of lg_node_anchor[u] .
        With that, node lg_node_anchor[u][0] is also mapped to lg_node_anchpr[u][1] ... and so forth.

        Additionally, the product graph nodes are filtered with respect to the lg_node_anchor such that no
        nodes outside of the anchor nodes in the product graph are included.
    """

    product_graph = nx.Graph()

    ## list of node lists
    node_list = [sorted(list(g.nodes)) for g in L]
    n_graphs = len(L)
    ## Computes the cartesian products of all the node sets of the graphs in L.
    ## NOTE: Instead of cartesian product of all node sets, do cartesian product on sets 
    ## of single bonds, double bonds etc. and unify these.
    product_nodes = list(itertools.product(*node_list))
    
    node_count = 1
    for i in range(n_graphs):
        node_count = node_count * len(node_list[i])

    ## unpack all anchor_nodes from the node anchor
    anchor_nodes = [ (v, *lg_node_anchor[v]) for v in lg_node_anchor]
    ## List of vectors. Each vector 'i' specifies which edges in graph L[i] are already included in the anchor.
    anchor_vector = [ [anchor_point[i] for anchor_point in anchor_nodes] for i in range(n_graphs) ]

    ## Firstly, nodes are only added if they're in the neighbourhood of any anchor point
    ## and are added based on their edge to an anchor
    for i in range(node_count):
        node_i = product_nodes[i]
        
        ## check for mixed
        contain_anchor_nodes = False 
        for i in range(n_graphs):
            if node_i[i] in anchor_vector[i]:
                contain_anchor_nodes = True
        
        ## Only look at nodes that either are anchor nodes, or do not try and match 
        ## an anchor node from one graph to a non-anchor node in other graph
        if (contain_anchor_nodes and node_i in anchor_nodes) or not contain_anchor_nodes:

            ## compute a list where neighbourhood[i] is the neighbourhood of node i
            neighbourhoods = []
            for L_i_node in range(n_graphs):
                ## Add the neighbourhood of v_i in product node from its corresponding graph L[i]
                neighbourhoods.append(L[L_i_node].adj[node_i[L_i_node]])
            
            ## for each anchor node, consider whether ALL (v_1, v_2, ..., v_n) are connected to the anchor 
            ## or ALL NOT connected to the anchor.
            for anchor_node in anchor_nodes:

                ## do not add "loops" on anchor nodes
                if node_i != anchor_node:
                    all_agree_connected = True
                    all_agree_disconnected = True
                    for index in range(n_graphs):
                        ## Determine node v_i's neighbourhood
                        v_i_neighbourhood = neighbourhoods[index]
                        
                        ## Look for the anchor in neighbourhood and 
                        ## exclude the other option when connection found/not found.
                        if anchor_node[index] in v_i_neighbourhood:
                            all_agree_disconnected = False
                        elif anchor_node[index] not in v_i_neighbourhood:
                            all_agree_connected = False
                        
                        ## stop looking as no consensus was found
                        if not all_agree_connected and not all_agree_disconnected: break
                        
                    ## add edge between this node and the anchor point
                    if all_agree_connected:
                        product_graph.add_edge( node_i, anchor_node, color="blue")

                    elif all_agree_disconnected:
                        product_graph.add_edge( node_i, anchor_node, color="red")

    ## Secondly, among all nodes in the neighbourhood, add internal edges
    added_nodes = list(product_graph.nodes)
    added_nodes_count = len(added_nodes)

    ## Add all edges between the nodes that are connected to the anchor
    for i in range(added_nodes_count):
        ## v = (v_1, v_2, ..., v_n)
        node_i = added_nodes[i]

        neighbourhoods = []
        for L_i_node in range(n_graphs):
            ## Add the neighbourhood of v_i in product node from its corresponding graph L[i]
            neighbourhoods.append(L[L_i_node].adj[node_i[L_i_node]])

        for j in range(i + 1, added_nodes_count):
            ## u = (u_1, u_2, ..., u_n)
            node_j = added_nodes[j]

            if not _has_node_in_common(node_i, node_j, n_graphs):

                all_agree_adj = True
                all_agree_not_adj = True 
                for index in range(n_graphs):
                    ## Determine node v_i's neighbourhood
                    v_i_neighbourhood = neighbourhoods[index]
                    
                    ## If node v_i is connected to u_i, then all nodes can't agree on being non-adjacent
                    if node_j[index] in v_i_neighbourhood:
                        all_agree_not_adj = False
                    ## Opposite argument
                    elif node_j[index] not in v_i_neighbourhood:
                        all_agree_adj = False
                    
                    ## Stop looking as no consensus was found
                    if not all_agree_adj and not all_agree_not_adj: break
                
                ## add edge between this node and the anchor point
                if all_agree_adj:
                    product_graph.add_edge( node_i, node_j, color="blue")

                elif all_agree_not_adj:
                    product_graph.add_edge( node_i, node_j, color="red")
    
    return product_graph