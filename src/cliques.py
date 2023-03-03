from productgraph import product_graph as pg
from linegraph import line_graph as lg
from linegraph import convert_edge_anchor_lg
import networkx as nx
from draw_graphs import draw_product_graph
from draw_graphs import draw_two_graphs

### Authors: Tobias Klink Lehn (toleh20@student.sdu.dk) and Kasper Halkjær Beider (kbeid20@student.sdu.dk)
def mcs_leviBarrowBurstall(G, H, node_anchor={}, edge_anchor={}):
    """
        Computes the Maximum Common Subgraph using the Algorithm suggested by
        G. Levi and H.G. Barrow + R.M. Burstall in 1973 and 1975 respectively.
        
        `Paramters`:
            G (Graph): A NetworkX graph, nodes are integers but may be decorated with items

            H (Graph): A NetworkX graph, nodes are integers but may be decorated with items

        `Optional`:
            node_anchor (dict: int -> int): A valid one-to-one mapping from 'n' nodes in G to 'n' nodes in H 

            edge_anchor (dict: int -> int): A valid one-to-one mapping from 'n' edges in G to 'n' edges in H, overrules a given node_anchor
        
        `Returns`:
            all_mappings (list( dict: node -> node) or list( dict: edge -> edge)): 
            
            A valid mapping from nodes in G to nodes in H if a node_anchor was given.
            Otherwise, a mapping from edges in G to edges in H.
            
    """

    ## Use Line-Graph implementation
    if edge_anchor:
        LG = lg(G)
        LH = lg(H)

        ## Mapping: v in LG -> u in LH
        computed_node_anchor = convert_edge_anchor_lg(G, H, edge_anchor)
        
        ## Compute product graph.
        mod_product_graph = pg(LG, LH)

        ## Anchor points in product graph correspond to (x, y) nodes where 
        ## x is a node in LG and y is a node in LH.
        anchor_points = [ (g, computed_node_anchor[g]) for g in computed_node_anchor ]
    
    ## No use of line graphs
    else:
        mod_product_graph = pg(G, H)

        anchor_points = [ (g, node_anchor[g]) for g in node_anchor]
    
    ## Finds all cliques containing the anchor points
    anchored_cliques = list(nx.find_cliques(mod_product_graph, anchor_points))

    ## Determines the size of the largest cliques for each mapping
    if not edge_anchor:
        max_clique_size = max(dict(nx.node_clique_number(mod_product_graph, cliques=anchored_cliques)).values())
    else:
        max_clique_size = max(dict(nx.node_clique_number(mod_product_graph, nodes=anchor_points, cliques=anchored_cliques)).values())
    
    all_max_cliques = [cliques for cliques in anchored_cliques if len(cliques) == max_clique_size]

    ## Convert the lists of cliques into lists of dictionaries
    if edge_anchor:
        ## In each clique, a node pair (i, j) corresponds to the i'th node
        ## in LG (i'th edge in G) being mapped to the j'th node in LH (j'th
        ## edge in H). 
        G_edges = list(G.edges)
        H_edges = list(H.edges)
        all_mappings = []
        for cliques in all_max_cliques:
            new_mapping = {}
            for points in cliques:
                new_mapping[G_edges[points[0]]] = H_edges[points[1]]
            all_mappings.append(new_mapping)
        
        return all_mappings
    
    else:
        ## In each clique, a node pair (i, j) corresponds to the i'th node
        ## in G being mapped to the j'th node in H. 
        all_mappings = []
        for cliques in all_max_cliques:
            new_mapping = {}
            for points in cliques:
                new_mapping[points[0]] = points[1]
            all_mappings.append(new_mapping)
        
        return all_mappings

if __name__ == "__main__":
    G = nx.Graph()
    G.add_edges_from([(0,1), (0,2), (1,2), (2,3), (3,4)])

    H = nx.Graph()
    H.add_edges_from([(0,1), (1,2), (1,3), (1,7), (2,3), (3,4), (3,6), (6, 7), (7, 8), (4,5)])

    edge_anchor = {
            (3, 4): (4, 5),
            (2, 3): (3, 4)
    }

    node_anchor = {
        4: 5
    }

    pgh = pg(G, H)

    print("RESULT FROM EDGE ANCHOR")
    subgraphs = mcs_leviBarrowBurstall(G, H, edge_anchor=edge_anchor)
    print(subgraphs)

    print("RESULT FROM NODE ANCHOR")
    subgraphs = mcs_leviBarrowBurstall(G, H, node_anchor=node_anchor)
    print(subgraphs)
    draw_two_graphs(G, H)
    


