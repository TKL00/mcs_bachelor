### Authors: Tobias Klink Lehn (toleh20@student.sdu.dk) and Kasper Halkjær Beider (kbeid20@student.sdu.dk)
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from copy import deepcopy
from workspace import Workspace

def mcs_mcgregor(G, H, anchor_point={}):
    '''
    Computes the Maximum Common Subgraph using the Algorithm suggested by
    James J. McGregor in 1982.

    `Precondition`: |V_G| <= |V_H|

        `Paramters`:
            G (Graph): A NetworkX graph
            H (Graph): A NetworkX graph

        `Optional`:
            anchor_point (dict: int -> int): A valid mapping from nodes in G to nodes in H
        
        `Returns`:
            mapping (dict: int -> int): The node correspondence for the MCS
            marcs (np.array): The MARCS array for the MCS
            arcsleft (int): The number of arcs in G that can be mapped to arcs in H

    If an anchor point is given (i.e. a subgraph isomorphism between G and H),
    the algorithm produces a common subgraph branching out from this anchor point.
    '''

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

    ## Contains (i, j) pairs when edges (i, j) were recently unable to correspond to each other
    killed_edges = []
    ## Number of total mappings found
    counter = 0
    bestarcsleft = 0
    arcsleft = G_edge_amt

    ## v in G, x in H
    v = 0
    x = None 
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
        ## If a vertex is found
        if x is not None:
            H_tried[v][x] = True
            H_mapped[x] = True
            
            ## UPDATE MARCS
            killed_edges = []
            v_edges = G_node_to_arc[v]
            x_edges = H_node_to_arc[x]
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
            ## 
            if arcsleft > bestarcsleft:
                if v == G_node_amt - 1:
                    best_mapping[0] = (deepcopy(current_mapping), deepcopy(MARCS), arcsleft)
                    counter += 1
                    
                    bestarcsleft = arcsleft
                else:
                    ## Store values in workspace associated with node v
                    workspaces[v] = Workspace(MARCS, arcsleft, MARCS_row_ones, killed_edges)
                    v += 1
                    ## New branch made, no edges killed and no nodes in H
                    ## have been tried in this branch.
                    killed_edges = []
                    H_tried[v] = [False for i in range(H_node_amt)]

        ## No node in H was found - backtracking is the only option.
        else:
            ## When backtracking, ignore the mapping of v in G to x in H.
            if current_mapping[v] != "":
                H_mapped[ current_mapping[ v ] ] = False
                current_mapping[v] = ""
            v -= 1

            if v < 0:
                (return_mapping, return_marcs, return_arcsleft) = best_mapping[0]
                return return_mapping, return_marcs, return_arcsleft
            ## Restore the saved workspace
            MARCS = workspaces[v].get_MARCS()
            arcsleft = workspaces[v].get_arcsleft()
            MARCS_row_ones = workspaces[v].get_MARCS_ones_left()
            killed_edges = workspaces[v].get_edges_killed()

if __name__ == "__main__":

    G = nx.Graph()
    G.add_edges_from([(0,1), (1,2), (2,0), (2,3)])
    H = nx.graph_atlas_g()[100]
    # H = nx.Graph()
    # H.add_edges_from([(0,3), (3,2), (2,1), (1,3)])

    mapping, marcs, arcsleft = mcs_mcgregor(G, H)
    print(mapping)


    # # DRAWING
    subax1 = plt.subplot(121)
    nx.draw(G, pos=nx.spring_layout(G), with_labels=True)
    subax2 = plt.subplot(122)
    nx.draw(H, node_color='r', with_labels=True)

    print("G EDGES: ", G.edges)
    print("H EDGES:", H.edges)

    plt.show()


