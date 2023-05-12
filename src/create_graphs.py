import networkx as nx
from draw_graphs import draw_one_graph
import random

def inspect_graph(path):

    g = nx.read_adjlist(path)
    draw_one_graph(g)


if __name__ == "__main__":

    # inspect_graph("../unlabelled_anchored_graphs/0.txt")

    # GRAPH_AMT = 10
    # SIZE_CLASSES = [5, 10, 20, 50, 100]

    
    # for size_class in SIZE_CLASSES:
    #     graphs_written = 0
    #     while graphs_written != GRAPH_AMT:

    #         random_number_edges = random.randint(int(size_class*1.5), size_class * 3)
    #         g = nx.dense_gnm_random_graph(size_class, random_number_edges)
    #         if nx.is_connected(g):
    #             path = "../unlabelled_graphs/" + str(size_class) + "_" + str(graphs_written) + ".txt"
    #             nx.write_adjlist(g, path)
    #             graphs_written += 1
                
            


    
    graphs_written = 0
    while graphs_written < 1:

        
        random_number_nodes = random.randint(8, 10)
        random_number_edges = random.randint(15, 20)
        graph = nx.dense_gnm_random_graph(random_number_nodes, random_number_edges)
        
        if nx.is_connected(graph):
            draw_one_graph(graph) 
            nx.write_adjlist(graph, "../unlabelled_anchored_graphs/" + "4" + ".txt")
        
            graphs_written += 1
            