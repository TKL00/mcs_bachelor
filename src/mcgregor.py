### Authors: Tobias Klink Lehn (toleh20@student.sdu.dk) and Kasper Halkjær Beider (kbeid20@student.sdu.dk)
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from workspace import Workspace

def mcs_mcgregor(G1, G2):
    '''
    Takes two NetworkX graphs as input and computes
    the Maximum Common Subraph using the Algorithm suggested by
    James J. McGregor in 1982.
    Precondition: |V_1| <= |V_2|
    '''
    

    G1_node_amt, G2_node_amt = len(G1.nodes), len(G2.nodes)
    G1_edge_amt, G2_edge_amt = len(G1.edges), len(G2.edges)

    G1_node_to_arc, G2_node_to_arc = node_to_arc_matrix(G1), node_to_arc_matrix(G2)

    ## Initialize the ARC Matrix (q1 x q2) to be full of 1's (all arcs are compatible)
    MARCS = np.ones((G1_edge_amt, G2_edge_amt))
    ## Counting array to determine whether a row in MARCS has been reduced to all zeros
    MARCS_row_ones_left = [G2_edge_amt for i in range(G1_edge_amt)]

    ## G2_mapped[i] = 1 if node "i" in G2 has already been assigned a vertex from G1
    G2_mapped = [False for i in range(G2_node_amt)]
    ## G2_tried[i][j] = 1 if node "j" in G2 has been tried for node "i" in G1
    G2_tried = np.zeros((G1_node_amt, G2_node_amt), dtype=bool)

    workspaces = [0 for i in range(G1_node_amt)]

    current_mapping = {}
    for node_index in range(G1_node_amt):
        current_mapping[node_index] = ""

    counter = 0
    bestarcsleft = 0
    arcsleft = G1_edge_amt
    v_i = 0
    x_i = None 
    while v_i >= 0: 

        # if v_i == 3:
        #     print("Current Mapping:", current_mapping)
        #     print("ARCS LEFT:", arcsleft, "\t BEST ARCS LEFT", bestarcsleft)
        #     print("MARCS:")
        #     print(MARCS)
        #     print()
            
        x_i = None
        for index in range(G2_node_amt):
            if not G2_tried[v_i][index] and not G2_mapped[index]:
                # Map v_i to x_i
                x_i = index
                if current_mapping[v_i] != "":
                    G2_mapped[ current_mapping[ v_i ] ] = False
                
                ## Reset MARCS rows for v_i
                v_i_edges = G1_node_to_arc[v_i]
                for v_i_edge_index in range(G1_edge_amt):
                    if v_i_edges[v_i_edge_index] == 1:
                        if MARCS_row_ones_left[v_i_edge_index] == 0:
                            arcsleft += 1
                        MARCS[v_i_edge_index] = [1 for i in range(G2_edge_amt)]
                        MARCS_row_ones_left[v_i_edge_index] = G2_edge_amt

                current_mapping[v_i] = x_i
                break
        ## If a vertex is found
        if x_i is not None:
            G2_tried[v_i][x_i] = True
            G2_mapped[x_i] = True

            ## UPDATE MARCS
            v_i_edges = G1_node_to_arc[v_i]
            x_i_edges = G2_node_to_arc[x_i]
            for v_i_edge_index in range(G1_edge_amt):
                ## v_i is incident to this arc
                if v_i_edges[v_i_edge_index] == 1:
                    for x_i_edge_index in range(G2_edge_amt):
                        ## For all edges not incident with x_i, v_i arcs cannot be mapped
                        if x_i_edges[x_i_edge_index] == 0:
                            if MARCS[v_i_edge_index][x_i_edge_index] != 0:
                                MARCS[v_i_edge_index][x_i_edge_index] = 0
                                MARCS_row_ones_left[v_i_edge_index] -= 1
                                if MARCS_row_ones_left[v_i_edge_index] == 0:
                                    arcsleft -= 1                            
            
            if arcsleft >= bestarcsleft:
                if v_i == G1_node_amt - 1:
                    print("CURRENT MAPPING:")
                    print(current_mapping)

                    print("MARCS:")
                    print(MARCS)
                    counter += 1
                    
                    bestarcsleft = arcsleft
                else:
                    ## store copy
                    workspaces[v_i] = Workspace(MARCS, arcsleft, MARCS_row_ones_left)
                    v_i += 1
                    G2_tried[v_i] = [False for i in range(G2_node_amt)]

        else:
            ## When backtracking after a mapping has occurred between v_i and x_i, undo the "map" status of x_i
            if current_mapping[v_i] != "":
                G2_mapped[ current_mapping[ v_i ] ] = False
                current_mapping[v_i] = ""
            v_i -= 1
            if v_i < 0:
                print("Counter:\t", counter)
                return
            MARCS, arcsleft, MARCS_row_ones_left = workspaces[v_i].MARCS, workspaces[v_i].arcsleft, workspaces[v_i].MARCS_ones_left

def node_to_arc_matrix(G):
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



if __name__ == "__main__":


    # G = nx.Graph()
    # H = nx.complete_graph(5)

    # G.add_nodes_from([i for i in range(10)])

    # a = [i for i in range(9, -1, -1)]
    # b = [i for i in range(0, 10)]
    # tuples = list(zip(a,b))
    # G.add_edges_from(tuples)
    # G.add_edges_from([(1,2), (3,4)])

    # G = nx.Graph()
    # G.add_edges_from([(0,1), (1,2), (2,3), (3,0)])
    # H = nx.complete_graph(3)

    H = nx.Graph()
    H.add_edges_from([(0,1), (1,2), (2,3), (3,0)])
    G = nx.complete_graph(5)

    # DRAWING
    subax1 = plt.subplot(121)
    nx.draw(G)
    subax2 = plt.subplot(122)
    nx.draw(H, node_color='r')

    plt.show()

    print("H EDGES:\t", H.edges)
    print("G EDGES:\t", G.edges)
    mcs_mcgregor(H, G)


