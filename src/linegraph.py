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

def convert_edge_anchor_lg_list(L, edge_anchor):
    """
    Computes the node_anchor of the line graphs made from the graphs in L based on edge_anchor.

    ``Parameters``:
        L (list (Graph)): A list of networkX graphs
        edge_anchor (list: list(edge)): A list of lists of edges. In edge_anchor[l], all edges are mapped to each other.
                                        edge_anchor[l][i] is an edge from graph L[i].

    ``Returns``:
        node_map ( dict: node -> list(node) ): The mapping from node in L[0] to nodes in the line graphs L[i] for i > 0.
    """

    node_map = {}

    n_graphs = len(L)
    ## list of edge lists from all graphs in L
    edgelist = [list(L[i].edges) for i in range(n_graphs)]

    ## Transform each [(u, v), (a, b), (x, y)] into [node_i, node_j, node_k] for every l
    node_map = [ [ edgelist[i].index(l[i]) for i in range(n_graphs)] for l in edge_anchor]

    return node_map

def convert_edge_anchor_lg(G, H, edge_anchor):
    """
    Assumes edge_anchor is from G_edge -> H_edge. Computes the node_anchor in line graphs of G and H by the given edge_anchor.
    
    ``Parameters``:

    ``Returns``:
        node_map ( dict: int -> int ): The mapping from node i to node j in the linegraph which corresponds to 
                                       edge i and edge j in the original graphs.
    """

    node_map = {}
    G_edges = list(G.edges)
    H_edges = list(H.edges)
    for keys in edge_anchor:
        G_edge = keys
        H_edge = edge_anchor[keys]
        LG_node = G_edges.index(G_edge)
        LH_node = H_edges.index(H_edge)
        node_map[LG_node] = LH_node

    return node_map

def convert_edge_anchor(G, H, edge_anchor):
    """
        Converts a given edge_anchor from G -> H into a node_anchor.
        That is, if edge_anchor has an entry (i, j) -> (a, b) it entails the mapping {i: a, j: b}.
        As a result, it is a precondition that the edge indices are sorted lexicographically to keep the
        integrity of the node-mapping intact when several edges share the same node.

        
        ( e.x. (1, 3) -> (a, c) entails {1: a, 3: c} but having (1, 2) -> (b, a) later results in the final mapping {1: b, 2: a, 3: c}
        when one wanted {1: a, 2: b, 3: c} )
    """

    node_anchor = {}
    for g_edges in edge_anchor:
        mapped_edge = edge_anchor[g_edges]

        node_anchor[g_edges[0]] = mapped_edge[0]
        node_anchor[g_edges[1]] = mapped_edge[1]

    return node_anchor
