import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from copy import deepcopy
from workspace import Workspace
from linegraph import line_graph as lg
from linegraph import convert_edge_anchor
from draw_graphs import draw_mcgregor_mcs_graphs
import molecules as ml

### Authors: Tobias Klink Lehn (toleh20@student.sdu.dk) and Kasper Halkjær Beider (kbeid20@student.sdu.dk)
def mcs_mcgregor(G, H, edge_anchor=[], molecule=False):
    """
    Computes the Maximum Common Subgraph using the Algorithm suggested by
    James J. McGregor in 1982.

    `Precondition`: |V_G| <= |V_H|

        `Paramters`:
            G (Graph): A NetworkX graph, nodes are integers but may be decorated with items
            H (Graph): A NetworkX graph, nodes are integers but may be decorated with items

        `Optional`:
            anchor_point (dict: int -> int): A valid one-to-one mapping from ``n`` edges in G to ``n`` edgs in H.
        
        `Returns`:
            all_mappings (list: (mapping, marcs, arcsleft)) where arcsleft is maximum:
                mapping (dict: int -> int): The node correspondence for the MCS
                marcs (np.array): The MARCS array for the MCS
                arcsleft (int): The number of arcs in G that can be mapped to arcs in H

    If an anchor point is given (i.e. a subgraph isomorphism between G and H),
    the algorithm produces a common subgraph branching out from this anchor point.
    """

    ## Initialization
    G_node_amt, H_node_amt = len(G.nodes), len(H.nodes)
    G_edge_amt, H_edge_amt = len(G.edges), len(H.edges)

    assert G_node_amt <= H_node_amt, f"The number of nodes in first input graph {G_node_amt} is larger than the number of nodes in the second input graph {H_node_amt}"

    if molecule:
        G_atom_types, H_atom_types = nx.get_node_attributes(G, "atom_type"), nx.get_node_attributes(H, "atom_type")
        G_bond_types, H_bond_types = nx.get_edge_attributes(G, "bond_type"), nx.get_edge_attributes(H, "bond_type")
    else: 
        G_atom_types, H_atom_types = {i: "" for i in G.nodes}, {i: "" for i in H.nodes}
        G_bond_types, H_bond_types = {i: "" for i in G.edges}, {i: "" for i in H.edges}

    ## Auxiliary function
    def node_to_arc_matrix(G):
        """
        Computes a |V| x |E| matrix with (v, e) = 1 if
        node `v` is incident with edge `e`.
        """
        V_size = len(G.nodes)
        A_size = len(G.edges)
        
        node_arc_matrix = np.zeros((V_size, A_size))
        edges = list(G.edges)

        ## Node pair (u, v) are both incident to the arc denoted by "Index"
        for index in range(len(G.edges)):
            (u, v) = edges[index]
            node_arc_matrix[u][index] = 1
            node_arc_matrix[v][index] = 1
        
        return node_arc_matrix
    
    def update_MARCS(MARCS, v_edges, x_edges, killed_edges, MARCS_row_ones, arcsleft):
        """
        Refines the MARCS matrix based on the edges connected to node v and node x.
        If edge 'i' is incident to node 'v' but edge 'j' is not incident with node 'x', then
        MARCS[i][j] will be set to 0. killed_edges and arcsleft will be adjusted accordingly.
        """
        for G_edge in range(G_edge_amt):
                ## v is incident to this arc.
                if v_edges[G_edge] == 1:
                    for H_edge in range(H_edge_amt):
                        ## For all edges not incident with x, v edge cannot be mapped.
                        if x_edges[H_edge] == 0:
                            if MARCS[G_edge][H_edge] != 0:
                                MARCS[G_edge][H_edge] = 0
                                ## Take note of lacking edge correspondence for future restoration.
                                killed_edges.append((G_edge, H_edge)) 
                                MARCS_row_ones[G_edge] -= 1
                                if MARCS_row_ones[G_edge] == 0:
                                    arcsleft -= 1
        return arcsleft

    def modify_MARCS_anchor(MARCS, H_mapped, G_anchor_nodes, H_anchor_nodes, current_mapping, MARCS_row_ones, arcsleft):
        for G_node in G_anchor_nodes:
            current_mapping[G_node] = "anchor"
        for H_node in H_anchor_nodes:
            H_mapped[H_node] = True

        ## refine MARCS such that no anchored edges can be mapped to anything
        for anchors in edge_anchor:
            G_edge = list(G.edges).index(anchors[0])
            H_edge = list(H.edges).index(anchors[1])

            ## The anchor edge in G can only be mapped to the anchored edge in H (G row is set to 0 except for the anchor in H)
            for i in range(H_edge_amt):
                if i != H_edge:
                    MARCS[G_edge][i] = 0
            ## Only one cell with the H-edge remains a 1
            MARCS_row_ones[G_edge] = 1

            ## H-column is set to 0 except for the anchor in G
            for i in range(G_edge_amt):
                if i != G_edge and MARCS[i][H_edge] == 1:
                    MARCS[i][H_edge] = 0
                    MARCS_row_ones[i] -= 1 ## A single cell in each row has been changed to 0
                    if MARCS_row_ones[i] == 0: arcsleft -= 1

            ## Update such that only edges incident with u, v can be mapped to incident edges to r, s
            (u, v) = anchors[0]
            G_incident = list(G.edges(u)) + list(G.edges(v))

            (r, s) = anchors[1]
            H_incident = list(H.edges(r)) + list(H.edges(s))

            ## For all incident edges to this given anchor in G, only incident edges to the anchor in H may be mapped.
            for g_edge in G_incident:
                (g1, g2) = g_edge
                try:
                    g_edge_index = list(G.edges).index((g1, g2))
                except:
                    g_edge_index = list(G.edges).index((g2, g1))
                for i in range(H_edge_amt):
                    (h1, h2) = list(H.edges)[i]
                    if (h1, h2) not in H_incident and (h2, h1) not in H_incident and MARCS[g_edge_index][i] == 1:
                        MARCS[g_edge_index][i] = 0
                        MARCS_row_ones[g_edge_index] -= 1
                        if MARCS_row_ones[g_edge_index] == 0: arcsleft -= 1
            
        return arcsleft

    def modify_MARCS_bond_type(MARCS, MARCS_row_ones, arcsleft, G_bond_types, H_bond_types):
        G_edges = list(G.edges)
        H_edges = list(H.edges)
        for i in range(G_edge_amt):
            for j in range(H_edge_amt):
                if G_bond_types[G_edges[i]] != H_bond_types[H_edges[j]]:
                    MARCS[i][j] = 0
                    MARCS_row_ones[i] -= 1
                    if MARCS_row_ones[i] == 0:
                        arcsleft -= 1
        return arcsleft

    G_node_to_arc, H_node_to_arc = node_to_arc_matrix(G), node_to_arc_matrix(H)

    G_anchor_nodes, H_anchor_nodes = convert_edge_anchor(edge_anchor)

    ## Initialize the ARC Matrix (q1 x q2) to be full of 1's (all arcs are compatible)
    MARCS = np.ones((G_edge_amt, H_edge_amt))

    ## Counting array to determine whether a row in MARCS has been reduced to all zeros
    ## used to maintain "arcsleft" instead of iterating through a row in MARCS every time a value is adjusted.
    MARCS_row_ones = [H_edge_amt for i in range(G_edge_amt)]

    ## H_mapped[i] = 1 if node "i" in H has already been assigned a vertex from G
    ## In place to ensure that a node in G in a branch isn't mapped
    ## to a node in H already occupied by a previous node in G
    H_mapped = [False for i in range(H_node_amt)]

    ## H_tried[i][j] = 1 if node "j" in H has been tried for node "i" in G
    ## In place to ensure that a node in G doesn't try to map to the same node
    ## multiple times.
    H_tried = np.zeros((G_node_amt, H_node_amt), dtype=bool)

    ## Array to save workspaces containing MARCS, arcsleft, MARCS_row_ones , and killed_edges
    workspaces = [0 for i in range(G_node_amt)]

    ## The current G -> H function, MARCSs and arcs_left
    current_mapping = {}
    all_mappings = []
    for node_index in range(G_node_amt):
        current_mapping[node_index] = ""

    ## Contains (i, j) pairs when edges (i, j) were recently unable to correspond to each other.
    killed_edges = []
    
    ## Number of total mappings found.
    counter = 0
    bestarcsleft = 0
    arcsleft = G_edge_amt

    ## Modify initial MARCS if an anchor point is present to reflect the anchored edges.
    if edge_anchor: 
        arcsleft = modify_MARCS_anchor(MARCS, H_mapped, G_anchor_nodes, H_anchor_nodes, current_mapping, MARCS_row_ones, arcsleft)
    ## If the graphs are labelled adjust MARCS such that already at the start it is not possible to map edges of differing
    ## bond types.

    if molecule:
        arcsleft = modify_MARCS_bond_type(MARCS, MARCS_row_ones, arcsleft, G_bond_types, H_bond_types)
        
    ## v in G, x in H
    v = 0
    x = None
    ## Find first node in G not a part of anchor point (lexicographic order)
    while v in G_anchor_nodes:
        v += 1

    ## In case all nodes in G are anchored, the algorithm needs not to run.
    if v >= G_node_amt:
        return current_mapping, MARCS, arcsleft
    first_non_anchor = v
    
    #####################################################################################################################
    ##                                              ALGORITHM BEGINS HERE                                              ##
    #####################################################################################################################
    while v >= 0: 
        x = None

        ## Finding a node x in H that has not already been mapped to.  
        ## Additionally, that node x in H has not been tried yet by node v in G
        for H_node in range(H_node_amt):
            if not H_tried[v][H_node] and not H_mapped[H_node] and G_atom_types[v] == H_atom_types[H_node]:
                x = H_node
                ## If v is currently mapped to a different node
                ## update said node to no longer be mapped.
                if current_mapping[v] != "":
                    H_mapped[ current_mapping[ v ] ] = False

                ## Upon mapping to a different node in H, the edges killed as a result
                ## of the previous mapping must now be restored
                for (G_edge, H_edge) in killed_edges:
                    if MARCS[G_edge][H_edge] == 0:
                        MARCS[G_edge][H_edge] = 1
                        if MARCS_row_ones[G_edge] == 0:
                            arcsleft += 1  
                        MARCS_row_ones[G_edge] += 1

                current_mapping[v] = x
                break
        ## Vertex is found
        if x is not None:
            H_tried[v][x] = True
            H_mapped[x] = True
            
            ## UPDATE MARCS
            killed_edges = []
            v_edges = G_node_to_arc[v]
            x_edges = H_node_to_arc[x]
            
            arcsleft = update_MARCS(MARCS, v_edges, x_edges, killed_edges, MARCS_row_ones, arcsleft)

            ## If the number of edges to be mapped is 'high', we either build further down the branch
            ## or save the current mapping if in a leaf node.
            if arcsleft > bestarcsleft:        ## == comes from building all MCS, even ones where arcsleft are equal
                if v == G_node_amt - 1:
                    all_mappings.append((deepcopy(current_mapping), deepcopy(MARCS), arcsleft))
                    bestarcsleft = arcsleft
                else:
                    ## Store values in workspace associated with node v
                    workspaces[v] = Workspace(MARCS, arcsleft, MARCS_row_ones, killed_edges)

                    non_anchored_node = v
                    v += 1
                    while v in G_anchor_nodes and v < G_node_amt:
                        v += 1
                    ## Skipping ahead resulted in jumping "out" of G, so we must backtrack to most
                    ## recent non-anchored node and save this mapping.
                    if v == G_node_amt:
                        all_mappings.append((deepcopy(current_mapping), deepcopy(MARCS), arcsleft))
                        bestarcsleft = arcsleft
                        v = non_anchored_node
                    ## New branch made, no edges killed and no nodes in H
                    ## have been tried in this branch.
                    else:
                        killed_edges = []
                        H_tried[v] = [False for i in range(H_node_amt)]

        ## No node in H was found - backtracking is the only option.
        else:
            ## When backtracking, ignore the tentative mapping of v in G to x in H.
            if current_mapping[v] != "":
                H_mapped[ current_mapping[ v ] ] = False
                current_mapping[v] = ""

            ## When backtracking the algorithm continues backtracking until it reaches
            ## a non anchor point.
            v -= 1
            while v in G_anchor_nodes and v > first_non_anchor:
                v -= 1
            
            ## If the algorithm has backtracked past the first non anchor point it means
            ## there are only anchor points left, therefore the algorithm stops.
            if v < first_non_anchor:
                ## return only max
                max_arcsleft = max(all_mappings, key=lambda items:items[2])[2]
                all_mappings_filtered = list(filter(lambda x: x[2] == max_arcsleft, all_mappings))
                return all_mappings
            ## Restore the saved workspace
            MARCS = workspaces[v].get_MARCS()
            arcsleft = workspaces[v].get_arcsleft()
            MARCS_row_ones = workspaces[v].get_MARCS_ones_left()
            killed_edges = workspaces[v].get_edges_killed()


def construct_cs(G, marcs):
    """
        Constructs the subgraph denoted by the mapped edges in G from MARCS.
    """
    subgraph = nx.Graph()
    G_edges = list(G.edges)

    for row in range(len(marcs)):
        for column in range(len(marcs[0])):
            if marcs[row][column] == 1:
                print(G_edges[row])
                subgraph.add_edge(G_edges[row][0], G_edges[row][1])

    return subgraph
