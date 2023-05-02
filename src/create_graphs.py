import networkx as nx
from draw_graphs import draw_one_graph
import random

def inspect_graph(path):

    g = nx.read_adjlist(path)
    draw_one_graph(g)


if __name__ == "__main__":

    inspect_graph("../unlabelled_graphs/0.txt")

    # GRAPH_AMT = 10

    # graphs_written = 0
    # for i in range(GRAPH_AMT):

    #     random_number_nodes = random.randint(10, 19)
    #     random_number_edges = random.randint(random_number_nodes, random_number_nodes * 3)
    #     g = nx.dense_gnm_random_graph(random_number_nodes, random_number_edges)
    #     if nx.is_connected(g):
    #         path = "../unlabelled_graphs/" + str(graphs_written) + ".txt"
    #         nx.write_adjlist(g, path)
    #         draw_one_graph(g)
    #         graphs_written += 1


    


    # graph = nx.dense_gnm_random_graph(15, 20)
    # draw_one_graph(graph) 
    # nx.write_adjlist(graph, "../unlabelled_graphs/test.txt")

    # g = nx.read_adjlist("../unlabelled_graphs/test.txt")
    

    # with open("../unlabelled_graphs/test.txt", "w") as f:
        