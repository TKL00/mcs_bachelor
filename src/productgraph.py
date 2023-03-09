import networkx as nx
import itertools

### Authors: Tobias Klink Lehn (toleh20@student.sdu.dk) and Kasper Halkjær Beider (kbeid20@student.sdu.dk)
def product_graph(G, H):
    """
        Computes the modular product of two NetworkX graphs, G and H.
    """

    product_graph = nx.Graph()

    ## Computes the cartesian products of the two node sets (sorted)
    product_nodes = itertools.product(sorted(G.nodes), sorted(H.nodes))
    product_graph.add_nodes_from(product_nodes)
    node_count = len(product_graph.nodes)

    product_graph_nodes = list(product_graph.nodes)

    for i in range(node_count):
        ## "u" and "x" from G, "v" and "y" from H
        (u, v) = product_graph_nodes[i]
        for j in range(i + 1, node_count):
            (x, y) = product_graph_nodes[j]

            if u != x and v != y:
                u_neighbours, v_neighbours = list(G.adj[u]), list(H.adj[v])

                if x in u_neighbours and y in v_neighbours:
                    product_graph.add_edge( (u, v), (x, y), color="blue" )
                
                elif x not in u_neighbours and y not in v_neighbours:
                    product_graph.add_edge( (u, v), (x, y), color="red" )

    return product_graph 