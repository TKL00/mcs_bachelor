import networkx as nx
import itertools

def _has_node_in_common(u, v, n_coordinates):
    for i in range(n_coordinates):
        if u[i] == v[i]:
            return True
    return False

### Authors: Tobias Klink Lehn (toleh20@student.sdu.dk) and Kasper Halkjær Beider (kbeid20@student.sdu.dk)
def product_graph_no_limit(L):
    """
        Computes the modular product of a list of graphs. Edges in the returned modular product are decorated with red/blue colors
        following the definition of the modular product.
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

def product_graph_limit(L, anchor_nodes, molecule=False):
    """
        Computes the modular product of all NetworkX graphs contained in L.
        With that, node lg_node_anchor[u][0] is also mapped to lg_node_anchpr[u][1] ... and so forth.

        Additionally, the product graph nodes are filtered with respect to the lg_node_anchor. Only nodes
        connected to an anchor point by a red/blue edge will be included. Blue/red edges among added
        nodes are also included.

        `Parameters`:
            L: List of graphs
            anchor_nodes: list of anchor nodes in the product graph of the form (v_1, v_2, ..., v_n)
            molecule (Boolean): Indicates whether the graphs are decorated with molecule attributes or not

        `Returns`:
            product_graph (Graph): A NetworkX graph that contains anchor nodes and all nodes connected to anchor. If the 
                                   anchor is only one node and no nodes are connected to the anchor, the product graph is empty.
                                   If the anchor consists of multiple nodes, the product graph will contain at least the anchor nodes.

    """

    def molecule_atom_bond_check(node, dimensions, atom_pair_attributes, bond_type_attributes):
        """
            Returns true if the given node agrees on bond types and atom pairs 
        """

        agree_on_atom_pair = True
        agree_on_bond_type = True
        ## set the target value to be the value of the first node
        target_atom_pair = atom_pair_attributes[0][node[0]]
        target_bond_type = bond_type_attributes[0][node[0]]

        for i in range(1, dimensions):
            if atom_pair_attributes[i][node[i]] != target_atom_pair:
                agree_on_atom_pair = False
                break
            if bond_type_attributes[i][node[i]] != target_bond_type:
                agree_on_bond_type = False
                break

        return agree_on_atom_pair and agree_on_bond_type

    product_graph = nx.Graph()

    if molecule:
        atom_pairs = [nx.get_node_attributes(graph, "atom_pair") for graph in L]
        bond_types = [nx.get_node_attributes(graph, "bond_type") for graph in L]

    ## list of node lists
    node_list = [sorted(list(g.nodes)) for g in L]
    n_graphs = len(L)
    ## Computes the cartesian products of all the node sets of the graphs in L.
    product_nodes = list(itertools.product(*node_list))

    ## filter product_nodes based on atom_pairs and bond_types if looking at a molecule
    if molecule:
        product_nodes = list(filter(lambda node: molecule_atom_bond_check(node, n_graphs, atom_pairs, bond_types), product_nodes))
        node_count = len(product_nodes)
    else:
        ## calculate number of potential nodes in the product graph without filtered for molecule attributes
        node_count = 1
        for i in range(n_graphs):
            node_count = node_count * len(node_list[i])
        
    ## List of vectors. Each vector 'i' specifies which edges in graph L[i] are already included in the anchor.
    anchor_vector = [ [anchor_point[i] for anchor_point in anchor_nodes] for i in range(n_graphs) ]

    ## Firstly, nodes are only added if they're in the neighbourhood of any anchor point
    ## and are added based on their edge to an anchor
    for i in range(node_count):
        node_i = product_nodes[i]

        ## check for mixed anchor
        contain_anchor_nodes = False 
        for i in range(n_graphs):
            if node_i[i] in anchor_vector[i]:
                contain_anchor_nodes = True
        
        ## Only look at nodes that consist completely of anchor nodes (i.e. not mixed w, anchor and non-anchor), 
        ## or nodes that are completely anchor-free.
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