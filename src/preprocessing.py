import networkx as nx
from queue import Queue
from draw_graphs import draw_one_graph
import copy 

def BFS_w_distance(G, anchored_nodes):

    source = anchored_nodes[0]
    node_distances = {}
    color_dict = {i: "w" for i in G.nodes}
    pred_dict = {i: None for i in G.nodes}

    node_distances[source] = 0
    color_dict[source] = "b"

    Q = Queue()
    Q.put(source)
    
    while not Q.empty():
        u = Q.get()
        neighbours = G.adj[u]
        for node in neighbours:
            if color_dict[node] == "w":
                color_dict[node] == "g"
                ## all anchor nodes have 0 distance in the graph
                if node in anchored_nodes:
                    node_distances[node] = 0
                else:
                    node_distances[node] = node_distances[u] + 1
                pred_dict[node] = u
                Q.put(node)
        color_dict[u] = "b"
    
    return node_distances

def anchor_reach(L, A):
    """
        L is a list of graphs
        A is a list of anchored edges in each graph. E.g. A[0] all anchored edges in L[0] etc.
        
        returns distance maps and shortest distance
    """

    ## anchored_nodes[0] is a list of anchored nodes from L[0]
    ## for every edge in A[i], iterate through each node and insert it
    anchored_nodes = { i: list(set([ node for edge in A[i] for node in edge])) for i in range(len(L)) }
    
    distance_map = {}
    
    for i in range(len(L)):
        graph = L[i]
        graph_anchors = anchored_nodes[i]
        ## Use first anchor node as source (irrelevant)
        node_distance = BFS_w_distance(graph, graph_anchors)
        distance_map[i] = node_distance

    max_distances = []
    for graph in distance_map:
        node_distances = distance_map[graph]
        max_distance = max(node_distances.values())
        max_distances.append(max_distance)
    
    shortest_distance = min(max_distances)

    return distance_map, shortest_distance
        
def shrink_graphs(L, shortest_distance, distance_map):
    """
        L is a list of graphs
        shortest_distance is the shortest max distance from any anchor in L
        distance_map is a dictionary s.t distance_map[0][i] is the distance from node i to the anchor in graph 0
    """
    ## Protection against introducing unintended side effects
    copied_graphs = [copy.deepcopy(graph) for graph in L]

    for i in range(len(copied_graphs)):
        graph = copied_graphs[i]
        distances = distance_map[i]
        
        ## Remove all nodes further away than the max distance
        for nodes in distances:
            if distances[nodes] > shortest_distance:
                graph.remove_node(nodes)

    return copied_graphs
    
        
    

    