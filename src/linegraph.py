import networkx as nx
import matplotlib.pyplot as plt

### Authors: Tobias Klink Lehn (toleh20@student.sdu.dk) and Kasper Halkjær Beider (kbeid20@st
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

    ## O(n^2)
    for i in range(len(G_edges)):
        for j in range(i + 1,len(G_edges)):
            if has_node_in_common(G_edges[i], G_edges[j]):
                LG.add_edge(i, j)

    return LG


# if __name__ == "__main__":

    
#     G = nx.Graph()
#     G.add_edges_from([(0,1), (1,2), (1,3), (1,7), (2,3), (3,4), (3,6), (4,5), (6,7), (7, 8)])
#     G_pos = nx.spring_layout(G, seed=1010203)

#     LG = line_graph(G)
#     LG_pos = nx.spring_layout(LG, seed=1010203)

#     subax1 = plt.subplot(121)
#     nx.draw(G, G_pos, with_labels=True)

#     subax2 = plt.subplot(122)
#     nx.draw(LG, LG_pos, with_labels=True)

#     plt.show()

    