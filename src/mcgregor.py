### Authors: Tobias Klink Lehn (toleh20@student.sdu.dk) and Kasper Halkjær Beider (kbeid20@student.sdu.dk)
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from copy import deepcopy
from workspace import Workspace

def mcs_mcgregor(G, H, anchor_point={}):
    """
    Computes the Maximum Common Subgraph using the Algorithm suggested by
    James J. McGregor in 1982.

    `Precondition`: |V_G| <= |V_H|

        `Paramters`:
            G (Graph): A NetworkX graph, nodes are integers but may be decorated with items
            H (Graph): A NetworkX graph, nodes are integers but may be decorated with items

        `Optional`:
            anchor_point (dict: int -> int): A valid one-to-one mapping from 'n' nodes in G to 'n' nodes in H
        
        `Returns`:
            mapping (dict: int -> int): The node correspondence for the MCS
            marcs (np.array): The MARCS array for the MCS
            arcsleft (int): The number of arcs in G that can be mapped to arcs in H

    If an anchor point is given (i.e. a subgraph isomorphism between G and H),
    the algorithm produces a common subgraph branching out from this anchor point.
    """

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

    ## Initialization
    G_node_amt, H_node_amt = len(G.nodes), len(H.nodes)
    G_edge_amt, H_edge_amt = len(G.edges), len(H.edges)

    assert G_node_amt <= H_node_amt, f"The number of nodes in first input graph {G_node_amt} is larger than the number of nodes in the second input graph {H_node_amt}"

    G_node_to_arc, H_node_to_arc = node_to_arc_matrix(G), node_to_arc_matrix(H)

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

    ## The current G -> H function as well as a list of all mappings, MARCSs and arcs_left
    current_mapping = {}
    best_mapping = [""]
    for node_index in range(G_node_amt):
        current_mapping[node_index] = ""

    ## Contains (i, j) pairs when edges (i, j) were recently unable to correspond to each other.
    killed_edges = []
    
    ## Number of total mappings found.
    counter = 0
    bestarcsleft = 0
    arcsleft = G_edge_amt

    ## Modify initial MARCS if an anchor point is present to reflect the anchored edges.
    if anchor_point:
        for G_node in anchor_point:
            mapped_node = anchor_point[G_node]
            ## Take note that the mapping from G_node to mapped_node exists
            ## in H_mapped and current_mapping.
            H_mapped[mapped_node] = True
            current_mapping[G_node] = mapped_node
            v_edges = G_node_to_arc[G_node]
            x_edges = H_node_to_arc[mapped_node]
            
            arcsleft = update_MARCS(MARCS, v_edges, x_edges, killed_edges, MARCS_row_ones, arcsleft)
        
        ## Reset killed edges as these are not supposed to be restored from the anchor point perspective.
        killed_edges = []

    ## v in G, x in H
    v = 0
    x = None
    ## Find first node in G not a part of anchor point (lexicographic order)
    while v in anchor_point:
        v += 1

    ## In case all nodes in G are anchored, the algorithm needs not to run.
    if v >= G_node_amt: return current_mapping, MARCS, arcsleft
    first_non_anchor = v
    
    ##                                              ALGORITHM BEGINS HERE
    while v >= 0: 
        x = None
        ## Finding a node x in H that has not already been mapped to.  
        ## Additionally, that node x in H has not been tried yet by node v in G
        for H_node in range(H_node_amt):
            if not H_tried[v][H_node] and not H_mapped[H_node]:
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
                                        
            ## 
            if arcsleft > bestarcsleft:
                if v == G_node_amt - 1:
                    best_mapping[0] = (deepcopy(current_mapping), deepcopy(MARCS), arcsleft)
                    counter += 1
                    
                    bestarcsleft = arcsleft
                else:
                    ## Store values in workspace associated with node v
                    workspaces[v] = Workspace(MARCS, arcsleft, MARCS_row_ones, killed_edges)

                    non_anchored_node = v
                    v += 1
                    while v in anchor_point and v < G_node_amt:
                        v += 1
                    ## Skipping ahead resulted in jumping "out" of G, so we must backtrack to most
                    ## recent non-anchored node and save this mapping.
                    if v == G_node_amt:
                        best_mapping[0] = (deepcopy(current_mapping), deepcopy(MARCS), arcsleft)
                        counter += 1
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
            while v in anchor_point and v > first_non_anchor:
                v -= 1
            
            ## If the algorithm has backtracked past the first non anchor point it means
            ## there are only anchor points left, therefore the algorithm stops.
            if v < first_non_anchor:
                (return_mapping, return_marcs, return_arcsleft) = best_mapping[0]
                return return_mapping, return_marcs, return_arcsleft
            ## Restore the saved workspace
            MARCS = workspaces[v].get_MARCS()
            arcsleft = workspaces[v].get_arcsleft()
            MARCS_row_ones = workspaces[v].get_MARCS_ones_left()
            killed_edges = workspaces[v].get_edges_killed()

def pretty_printer_mapping(mapping):
    print(mapping)
    for v in mapping:
        print(v, "->", mapping[v])

def pretty_printer_graphs(G, H, mapping, marcs, anchor):
    
    ## String used for nodes mapped to each other
    label_string="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
    G_labels = {}
    H_labels = {}

    ## Assigning different labels for the nodes in the mapping
    ## a mapping u -> v from G to H assigns the same string to u and v 
    for g_nodes in mapping:
        G_labels[g_nodes] = label_string[g_nodes]
        H_labels[mapping[g_nodes]] = label_string[g_nodes]

    G_pos = nx.spring_layout(G, seed=1010203)
    H_pos = nx.spring_layout(H, seed=1010203)

    # #                                                 DRAWING/PLOTTING
    
    ## Options for node drawing
    options = {"edgecolors": "tab:gray", "node_size": 300, "alpha": 1}
    
    # H node colors
    H_node_colors = ["tab:pink" for i in range(len(H.nodes))]
    G_node_colors = ["tab:blue" for i in range(len(G.nodes))]

    ## Color anchor point in G
    for g_nodes in anchor:
        G_node_colors[g_nodes] = 'tab:orange'

    ## G DRAWING
    subax1 = plt.subplot(121)
    nx.draw_networkx_nodes(G, G_pos, nodelist=G.nodes, node_color=G_node_colors, **options)
    nx.draw_networkx_edges(G, G_pos, width=1.0, alpha=0.5)
    
    ## Determine which edges exist in the mapping and color those.
    G_edges = list(G.edges)
    G_color_edges = []
    H_edges = list(H.edges)
    H_color_edges = []
    for i in range(len(marcs)):
        for j in range(len(marcs[0])):
            if marcs[i][j] == 1:
                G_color_edges.append(G_edges[i])
                H_color_edges.append(H_edges[j])

    nx.draw_networkx_edges(
        G,
        G_pos,
        edgelist=G_color_edges,
        width=8,
        alpha=0.5,
        edge_color="tab:blue",
    )

    ## Add labels to mapped nodes from G to H
    ## such that the mapped nodes have the same label
    iterator = 0
    H_node_labels = {}
    G_node_labels = {}
    for nodes in G.nodes:
        if nodes in mapping:
            G_node_labels[nodes] = label_string[iterator]
            H_node_labels[mapping[nodes]] = label_string[iterator]
            ## mapped nodes in H will have the same color as the nodes in G
            H_node_colors[mapping[nodes]] = "tab:blue"
            iterator += 1
        ## otherwise, G nodes will be numbered from 0-n
        else:
            G_node_labels[nodes] = nodes
    ## same labeling in H with non-mapped nodes having numbers from 0-n if not already mapped.
    for nodes in H.nodes:
        if nodes not in H_node_labels:
            H_node_labels[nodes] = nodes
    
    ## Color anchor point in H AFTER mapping to overwrite color
    for g_node in anchor:
        H_node_colors[anchor[g_node]] = "tab:orange"

    nx.draw_networkx_labels(G, G_pos, G_node_labels, font_size=15, font_color="whitesmoke")


    ## H DRAWING
    subax2 = plt.subplot(122)
    nx.draw_networkx_nodes(H, nodelist=H.nodes, node_color=H_node_colors, pos=H_pos)
    nx.draw_networkx_edges(H, H_pos, width=1.0, alpha=0.5)
    nx.draw_networkx_edges(
        H,
        H_pos,
        edgelist=H_color_edges,
        width=8,
        alpha=0.5,
        edge_color="tab:blue",
    )
    nx.draw_networkx_labels(H, H_pos, H_node_labels, font_size=15, font_color="whitesmoke")

    plt.show()

if __name__ == "__main__":

    G = nx.Graph()
    G.add_edges_from([(0,1), (1,2), (2,0), (2,3)])
    anchor = dict([(0,10)])

    H = nx.complete_graph(10)
    H.add_edge(9,10)
    # H = nx.graph_atlas_g()[120]
    print("G EDGES: ", G.edges)
    print("H EDGES:", H.edges)

    mapping, marcs, arcsleft = mcs_mcgregor(G, H, anchor)
    pretty_printer_mapping(mapping)

    ## DRAWING
    pretty_printer_graphs(G, H, mapping, marcs, anchor)
