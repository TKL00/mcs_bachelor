import networkx as nx

### Authors: Tobias Klink Lehn (toleh20@student.sdu.dk) and Kasper Halkjær Beider (kbeid20@student.sdu.dk)
def line_graph(G):
    """
    Creates the Line Graph representation, L(G) of G.    
    L(G) is constructed in the following way: for each edge in G, make a vertex in L(G); 
    for every two edges in G that have a vertex in common, make an edge between their corresponding vertices in L(G).

        `Parameters`:
            G (Graph): A NetworkX graph.
        
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

    for i in range(len(G_edges)):
        LG.add_node(i)

    ## O(n^2)
    for i in range(len(G_edges)):
        for j in range(i + 1,len(G_edges)):
            if has_node_in_common(G_edges[i], G_edges[j]):
                LG.add_edge(i, j)

    ## add additional decorations to nodes

    return LG

def convert_edge_anchor_lg_list(L, edge_anchor):
    """
    Computes the node_anchor of the line graphs of all graphs in L based on edge_anchor.

    ``Parameters``:
        L (list (Graph)): A list of networkX graphs
        edge_anchor (dict: edge -> [edge]): A dictionary mapping edges from L[0] to all edges in in their edge list (all of these edges are
                                            transitively and symmetrically mapped) 

    ``Returns``:
        node_map ( dict: node -> list(node) ): The mapping from node in L[0] to nodes in the line graphs L[i] for i > 0.
    """

    node_map = {}

    n_graphs = len(L)
    ## list of edge lists from all graphs in L
    edgelist = [list(L[i].edges) for i in range(n_graphs)]

    for keys in edge_anchor:
        L_0_edge = keys
        L_0_edge_index = edgelist[0].index(L_0_edge)
        mapped_nodes = edge_anchor[L_0_edge]
        mappings = []
        for graphs in range(1, n_graphs):
            ## find node index of the current edge in the mapping
            L_i_edge = mapped_nodes[graphs - 1]
            L_i_edge_index = edgelist[graphs].index(L_i_edge)
            mappings.append(L_i_edge_index)

        ## map the L[0] edge index (which is the node number in the line graph) to the list of nodes
        node_map[L_0_edge_index] = mappings

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
