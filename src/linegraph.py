import networkx as nx

### Authors: Tobias Klink Lehn (toleh20@student.sdu.dk) and Kasper Halkjær Beider (kbeid20@student.sdu.dk)
def line_graph(G, molecule=False):
    """
    Creates the Line Graph representation, L(G) of G.    
    L(G) is constructed in the following way: for each edge in G, make a vertex in L(G); 
    for every two edges in G that have a vertex in common, make an edge between their corresponding vertices in L(G).

        `Parameters`:
            G (Graph): A NetworkX graph.
        
        `Optional`:
            molecule (Boolean): A boolean optional to specify whether the nodes are decorated with attributes
                atom_pair: A set of atoms that a node connects e.g. {C, O}
                bond_type: The type of bond of an edge in ['s', 'd', 't', 'q']
        
        `Returns`:
            LG (Graph): A NetworkX graph in which node 'i' in LG corresponds to edge 'i' in G. 
    """

    def has_node_in_common(e1, e2):
        """
            Precondition: e1 and e2 are not identical
        """
        (x1, y1) = e1
        (x2, y2) = e2

        return x1 == x2 or x1 == y1 or y1 == x2 or y1 == y2

    LG = nx.Graph()

    G_edges = list(G.edges)
    
    if molecule:
        node_attributes = nx.get_node_attributes(G, "atom_type")
        edge_attributes = nx.get_edge_attributes(G, "bond_type")

    for i in range(len(G_edges)):
        ## Adding molecule-related labels for the nodes in the linegraph
        if molecule:
            (u, v) = G_edges[i]
            u_atom_type = node_attributes[u]
            v_atom_type = node_attributes[v]
            LG.add_node(i, atom_pair=set( [u_atom_type, v_atom_type] ), bond_type=edge_attributes[(u, v)] )
        else:
            LG.add_node(i)

    ## O(n^2)
    for i in range(len(G_edges)):
        for j in range(i + 1,len(G_edges)):
            if has_node_in_common(G_edges[i], G_edges[j]):
                LG.add_edge(i, j)

    return LG

## used in Cliques
def convert_edge_anchor_lg_list(L, edge_anchor):
    """
    Computes the node_anchor of the line graphs made from the graphs in L based on edge_anchor.

    ``Parameters``:
        L (list (Graph)): A list of networkX graphs
        edge_anchor (list: list(edge)): A list of lists of edges. In edge_anchor[l], all edges are mapped to each other.
                                        edge_anchor[l][i] is an edge from graph L[i].

    ``Returns``:
        node_map ( list (list: nodes) ): The mapping of nodes in the linegraphs. An element is thus a list of edges mapped to each other.
    """

    node_map = {}

    n_graphs = len(L)
    ## list of edge lists from all graphs in L
    edgelist = [list(L[i].edges) for i in range(n_graphs)]

    ## Transform each [(u, v), (a, b), (x, y)] into [node_i, node_j, node_k] for every l
    node_map = [ [ edgelist[i].index(l[i]) for i in range(n_graphs)] for l in edge_anchor]

    return node_map

## used in McGregor
def convert_edge_anchor(edge_anchor):
    """
        Computes the list of anchored nodes in G and the list of anchored nodes in H.
    """
    G_anchor = set()
    H_anchor = set()

    for i in range(len(edge_anchor)):
        (u_1, v_1) = edge_anchor[i][0]
        (u_2, v_2) = edge_anchor[i][1]
        G_anchor.update([u_1, v_1])
        H_anchor.update([u_2, v_2])


    return list(G_anchor), list(H_anchor)
