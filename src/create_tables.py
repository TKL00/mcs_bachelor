from mcgregor import mcs_mcgregor
from draw_graphs import draw_mcgregor_mcs_graphs, draw_graphs, draw_one_graph
from graph_format import compute_anchor
from cliques import iterative_approach
from graph_format import convert_graph_file
import networkx as nx
import multiprocessing
import os
import time
import itertools as it
from preprocessing import shrink_graphs, anchor_reach

# test_signals.py
# only available in unix systems
import functools
import signal
import time
def timeout(max_timeout, default=None):
    """Timeout decorator, parameter in seconds."""
    def timeout_decorator(func):
        """Wrap the original function."""
        @functools.wraps(func)
        def func_wrapper(*args, **kwargs):
            """Timeout using signal."""
            class MyTimeoutError(Exception):
                pass
            def handler(signum, frame):
                raise MyTimeoutError(
                    "{0} - Timeout after {1} seconds".format(func.__name__, max_timeout)
                )
            # set the timeout handler
            signal.signal(signal.SIGALRM, handler)
            signal.alarm(max_timeout)
            result = default
            try:
                result = func(*args, **kwargs)
            except MyTimeoutError as exc:
                # handle the timeout
                print("", end="")
            finally:
                # Cancel the timer
                signal.alarm(0)
            return result
        return func_wrapper
    return timeout_decorator


MC_GREGOR_TIMEOUT = 10
CLIQUES_TIMEOUT = 30

@timeout(MC_GREGOR_TIMEOUT)
def call_mcgregor(g1, g2):
    return mcs_mcgregor(g1, g2)

@timeout(CLIQUES_TIMEOUT)
def call_cliques(graphs, anchor, bool1, bool2):
    return iterative_approach(graphs, anchor, bool1, bool2)

def mcgregor_same_class(list_of_graphs):
    for i in range(len(list_of_graphs)):
        g1 = list_of_graphs[i]
        for j in range(i + 1, len(list_of_graphs)):
            g2 = list_of_graphs[j]
            time_before = time.time()
            res = call_mcgregor(g1, g2)
            time_after = time.time()
            if res:
                print(f"{len(g1.nodes)} \t {len(g2.nodes)} \t {time_after-time_before}")
            else:
                print(f"{len(g1.nodes)} \t {len(g2.nodes)} \t timed out after {MC_GREGOR_TIMEOUT} seconds")
            

def mcgregor_across_class(list_of_graphs_c1, list_of_graphs_c2):
    for i in range(len(list_of_graphs_c1)):
        g1 = list_of_graphs_c1[i]
        for j in range(len(list_of_graphs_c2)):
            g2 = list_of_graphs_c2[j]
            time_before = time.time()
            res = call_mcgregor(g1, g2)
            time_after = time.time()
            if res:
                print(f"{len(g1.nodes)} \t {len(g2.nodes)} \t {time_after-time_before}")
            else:
                print(f"{len(g1.nodes)} \t {len(g2.nodes)} \t timed out after {MC_GREGOR_TIMEOUT} seconds")

def fix_edge_ordering(G):
    edge_set = list(G.edges)
    sorted_edge_set = [sorted(i) for i in edge_set]
    new_G = nx.Graph()
    new_G.add_nodes_from(sorted([i for i in G.nodes]))
    new_G.add_edges_from(sorted_edge_set)

    return new_G


def recursive_unlabelled_anchor(indices, graphs, anchored_edges):
    anchor = []
    for i in range(len(anchored_edges[0])):
        new_anchor = []
        for index in indices:
            new_anchor.append(anchored_edges[index][i])
        anchor.append(new_anchor)
    
    graphs_to_input = [graphs[i] for i in indices]

    time_before = time.time()
    res = call_cliques(graphs_to_input, anchor, True, False)
    time_after = time.time()
    if res:
        max_size = max([len(i) for i in res])
        print(*indices, sep=" ",end="")
        print(f"\t{round(time_after-time_before, ndigits=2)}\t{max_size}")
    else:
        print(*indices, sep=" ",end="")
        print(f"\t timed out after {CLIQUES_TIMEOUT} seconds")

## MCGREGOR
def table1():
    
    path = "../unlabelled_graphs"
    all_graphs = []
    for (root, dirs, files) in os.walk(path):
        
        for file_name in files:
            full_file_path = os.path.join(path, file_name)
            g = nx.read_adjlist(full_file_path, nodetype=int)
            all_graphs.append(g)
    
    five = list(filter(lambda g: len(g.nodes) == 5, all_graphs))
    ten = list(filter(lambda g: len(g.nodes) == 10, all_graphs))
    twenty = list(filter(lambda g: len(g.nodes) == 20, all_graphs))

    print(f"g1 n\tg2 n\ttime (s)")
    mcgregor_same_class(five)
    mcgregor_same_class(ten)
    mcgregor_same_class(twenty)

    mcgregor_across_class(five, ten)
    mcgregor_across_class(ten, twenty)
            
    return None


## CLIQUES
def table3():

    path = "../unlabelled_anchored_graphs"

    all_graphs = []
    anchored_edges = []
    for (root, dirs, files) in os.walk(path):
        
        anchor_files = sorted(list(filter(lambda name: "anchor" in name, files)))
        graph_files = sorted(list(filter(lambda name: name.endswith(".txt") and "anchor" not in name, files)))
        for graph_file in graph_files:
            full_file_path = os.path.join(path, graph_file)
            if full_file_path.endswith(".txt"):
                g = nx.read_adjlist(full_file_path, nodetype=int)
                all_graphs.append(g)
        
        for anchor_file in anchor_files:
            anchor_path = os.path.join(path, anchor_file)
            with open(anchor_path) as f:
                graph_anchored_edges = []
                for lines in f:
                    indices = lines.split(",")
                    graph_anchored_edges.append((int(indices[0]), int(indices[1])))
                anchored_edges.append(graph_anchored_edges)
    
    
    fixed_graphs = [fix_edge_ordering(G) for G in all_graphs]

    indices = [i for i in range(len(all_graphs))]

    res = []
    for i in range(2, len(indices)):
        res.append(list(it.permutations(indices, i)))

    all_orders = list(it.chain(*res))

    print(f"graph seq\ttime (s)\tmax extension")
    for order in all_orders:
        recursive_unlabelled_anchor(list(order), fixed_graphs, anchored_edges)

    return None 


## Gathering data from the graphs supplied by Daniel
def table4():
    print(f"filename\tnumber of graphs\tlargest # nodes\tgraph seq\tmax distance\tmax extension\ttime(s)")
    
    clique_seq = [
           [2, 1, 0, 3],
           [2, 1, 0, 3],
           [2, 1, 0, 3],
           [2, 1, 0, 3],
           [0, 2, 1],
           [0, 3, 1, 2],
           [0, 2, 1, 3],
           [0, 2, 1, 3],
           [0, 2, 1],
           [0, 3, 2, 1],
           [0, 1, 2, 3]
    ]
    
    path = "../labelled_graphs"
    all_graphs = []
    all_anchors = []
    file_names = []
    for (root, dirs, files) in os.walk(path):
        for file_name in files:
            file_names.append(file_name)
            
            full_file_path = os.path.join(path, file_name)
            graphs, anchors = convert_graph_file(full_file_path)
            all_graphs.append(graphs)
            all_anchors.append(anchors)
    
    distance_classes = [6, 7, 8, 9]
    for i in range(len(all_graphs) - 1):
        
        ## print without shrinkage
        file_name = file_names[i]
        Gs = all_graphs[i]
        max_size = max([len(g.nodes) for g in Gs])
        n_graphs = len(Gs)
        As = all_anchors[i]
        seq = clique_seq[i]
        graph_seq = [Gs[i] for i in seq]
        anchor_seq = [As[i] for i in seq]
        
        graph_anchor = compute_anchor(graph_seq, anchor_seq, True)
        chosen_anchor = graph_anchor[0]
        

        time_before = time.time()
        res = call_cliques(graph_seq, chosen_anchor, True, True)
        time_after = time.time()
        if res:
            max_length = max([len(i) for i in res])
            print(f"{file_name}\t{n_graphs}\t{max_size}\t{seq}\tinfinity\t{max_length}\t{time_after-time_before} ")
        else:
            print(f"{file_name}\t{n_graphs}\t{max_size}\t{seq}\tinfinity\t-\ttimed out after {CLIQUES_TIMEOUT} seconds ")
        
        
        dist_map, shortest_distance = anchor_reach(graph_seq, chosen_anchor)
        
        ## print distance class results
        for dist_class in distance_classes:
            shrunk_graphs = shrink_graphs(graph_seq, dist_class, dist_map)

            if "ethanol_backward.txt" in file_name and dist_class == 7:
                print("In if statement create_tables")
                print(shrunk_graphs)
                for i in shrunk_graphs:
                    print("In for loop create_tables")
                    draw_one_graph(i)
            time_before = time.time()
            try:
                res = call_cliques(shrunk_graphs, chosen_anchor, True, True)
                time_after = time.time()
                if res:
                    max_length = max([len(i) for i in res])
                    print(f"{file_name}\t{n_graphs}\t{max_size}\t{seq}\t{dist_class}\t{max_length}\t{time_after-time_before} ")
                else:
                    print(f"{file_name}\t{n_graphs}\t{max_size}\t{seq}\t{dist_class}\t-\ttimed out after {CLIQUES_TIMEOUT} seconds ")
            except:
                print("ERROR")

            
        
    return None 


if __name__ == "__main__":
    # table1()
    # table3()
    table4()