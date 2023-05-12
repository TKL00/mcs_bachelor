from productgraph import product_graph_no_limit as pgnl
from productgraph import product_graph_limit as pgl
from linegraph import line_graph as lg
from linegraph import convert_edge_anchor_lg_list
from itertools import chain
import graph_format
import numpy as np
import networkx as nx
import networkx.algorithms.isomorphism as iso
from queue import Queue
import copy
from draw_graphs import draw_product_graph, draw_graphs, draw_blue_connected_components, draw_molecules, draw_one_graph
from preprocessing import anchor_reach, shrink_graphs
from molecules import *
import time


### Authors: Tobias Klink Lehn (toleh20@student.sdu.dk) and Kasper Halkjær Beider (kbeid20@student.sdu.dk)
def mcs_list_leviBarrowBurstall(L, edge_anchor, limit_pg=True, molecule=False):
    """
        Computes the Maximum Common Subgraph using the Algorithm suggested by
        G. Levi and H.G. Barrow + R.M. Burstall in 1973 and 1975 respectively.
        
        `Paramters`:
            L (list(Graph)): A list of networkX graphs.

            edge_anchor (list: list(edge)): A valid one-to-one between edges in the graphs in L. An element X in the edge_anchor 
                                            is thus a list of edges, where X[i] is an edge in L[i]. All edges in X are mapped to each other.

        `Optional`:
            limit_pg (boolean): Indicates whether the product graph should be limited to the neighbourhood of anchors or not, default to true

            molecule (boolean): Indicates whether the graphs in L are decorated molecules. If true, it is expected that each graph has
                                attribute "atom_type" on nodes and "bond_type" on edges. This further limits the tuples in the product graph.
                                Default to false.
        
        `Returns`:
            all_mappings (list( list (list(edge))): Lists of mapping. That is, each element in all_mappings is a list of edge mappings between the graphs
                                                    following the same format as edge_anchor.
            
            A valid mapping from nodes in G to nodes in H if a node_anchor was given.
            Otherwise, a mapping from edges in G to edges in H.
            
    """

    def blue_component_filter(PG, N, A, edge_color_lookup):
        """
            Given the "blue" neighbourhood N of the anchored nodes in A in the Product Graph PG, computes a list of
            disjoint blue connected components reachable from the anchor point via blue edges. Blue connected components
            not reachable from A via blue edges are dismissed.

            Precondition: All nodes in A are connected to all nodes in N either by blue or red edges.

            `Parameters`:

                PG (Graph): A NetworkX graph, nodes are tuples of integers (u, v) but may be decorated with items
                
                N (list(node)): A list of nodes in the blue neighbourhood of A. That is, if a node is in N, it is adjacent to all nodes in A.

                A (list(node)): A list of anchored nodes in PG.

                edge_color_lookup (dict(edge -> str)): A dictionary mapping edges in PG to a color. The value may either be "red" or "blue"
            
            `Returns`:

                all_components (list( list(node) )): A list of disjoint blue connected components all reachable by blue edges from the anchor, A.
        """
        ## BFS to find all nodes that are connected to A through blue edges
        ## all nodes in the product graph are white in the beginning
        ## NOTE: optimization for not whitening all the nodes in PG, but only the ones in A and N.
        
        ## Avoiding all nodes not in A and N by giving them a scary color
        node_color_lookup = { node: "p" for node in PG.nodes}                               ## SNAK OM DET HER PÅ ET TIDSPUNKT
        for nodes in A:
            node_color_lookup[nodes] = "w"
        for nodes in N:
            node_color_lookup[nodes] = "w"
        
        ## Starting from a node in the anchor and running BFS from that, 
        ## modified to only look at the neighbours that can be reached via a blue edge
        source = A[0]
        Q = Queue()
        Q.put(source)
        while not Q.empty():
            u = Q.get()
            neighbours = PG.adj[u]
            for v in neighbours:
                if node_color_lookup[v] == "w":
                    edge_color = edge_color_lookup[(u, v)] if (u, v) in edge_color_lookup else edge_color_lookup[(v, u)]                    
                    ## Don't consider nodes that are reached by red edges
                    if edge_color == "blue":
                        node_color_lookup[v] = "g"
                        Q.put(v)
            node_color_lookup[u] = "b"

        ## all nodes with a "black" color can be reached from A via blue edges. The reached nodes are either in A or N.
        marked_nodes = [nodes for nodes in node_color_lookup if node_color_lookup[nodes] == "b"]
        
        ## disregard nodes in the neighbourhood not reachable by blue edges
        filtered_N = [nodes for nodes in N if nodes in marked_nodes]
        
        ## Only the anchor
        if len(filtered_N) == 0:
            return []

        ## BFS on filtered nodes to find their connected components
        ## reset color of blue nodes
        for blue_nodes in filtered_N:
            node_color_lookup[blue_nodes] = "w"

        all_components = []
        ## continue BFS until all components are found (i.e. when all nodes of the filtered_N have been re-colored)
        while len(filtered_N) > 0:
            blue_source = filtered_N.pop()
            Q.put(blue_source)
            ## New component consists of at least this one node from the blue neighbourhood.
            current_component = [blue_source]        
            while not Q.empty():
                u = Q.get()
                neighbours = PG.adj[u]
                for v in neighbours:
                    if node_color_lookup[v] == "w":
                        edge_color = edge_color_lookup[(u, v)] if (u, v) in edge_color_lookup else edge_color_lookup[(v, u)]
                        
                        ## Don't consider nodes that are reached by red edges
                        if edge_color == "blue":
                            node_color_lookup[v] = "g"
                            current_component.append(v)
                            Q.put(v)
                node_color_lookup[u] = "b"
            ## remaining blue neighbours who haven't built their own component yet. 
            filtered_N = [nodes for nodes in filtered_N if node_color_lookup[nodes] == "w"]
            all_components.append(current_component)
            current_component = []

        return all_components

    def connected_MCS(listN, PG, A, edge_color_lookup):
        """
            Implementation of Algorithm 1 as suggested by Akbar Davoodi.

            Finds all maximal cliques in all blue connected components and creates
            the union of these with the anchor point.
        """
        MCSs = []
        ## Create node induced subgraph of the modular product graph containing all blue connected components
        blue_component_graph = nx.Graph(nx.induced_subgraph(PG, chain(*listN)))
        
        component_cliques = list(nx.find_cliques(blue_component_graph))

        ## For each clique, create the induced subgraph concatenated with the anchor
        ## do BFS on this new graph, removing all nodes that are not reachable by a blue edge from the anchor
        for comp_clique in component_cliques:
            bfs_graph = nx.Graph(nx.induced_subgraph(PG, A + comp_clique))
            node_color_lookup = {node: 'w' for node in bfs_graph.nodes}
            edge_color_lookup = nx.get_edge_attributes(bfs_graph, "color")
            source = A[0]
            Q = Queue()
            Q.put(source)
            while not Q.empty():
                u = Q.get()
                neighbours = bfs_graph.adj[u]
                for v in neighbours:
                    if node_color_lookup[v] == "w":
                        edge_color = edge_color_lookup[(u, v)] if (u, v) in edge_color_lookup else edge_color_lookup[(v, u)]                    
                        ## Don't consider nodes that are reached by red edges
                        if edge_color == "blue":
                            node_color_lookup[v] = "g"
                            Q.put(v)
                node_color_lookup[u] = "b"
            
            ## 
            reachable_nodes = [node for node in comp_clique if node_color_lookup[node] == 'b']
            MCSs.append(reachable_nodes + A)


        return MCSs

    n_graphs = len(L)
    edge_lists = [list(L[i].edges) for i in range(n_graphs)]

    ## Use Line-Graph implementation

    linegraphs = [lg(L[i], molecule=molecule) for i in range(n_graphs)]

    ## List of anchor nodes in the linegraphs
    computed_node_anchor = convert_edge_anchor_lg_list(L, edge_anchor)
    ## Unfold anchor nodes in the modular product graph
    anchor_nodes =  [ tuple(v) for v in computed_node_anchor]
    
    ## Compute product graph, either constraining it to only include A and N, or include all possible nodes.
    mod_product_graph = pgl(linegraphs, anchor_nodes, molecule=molecule) if limit_pg else pgnl(linegraphs)

    ## If no nodes are added, |anchor| = 1 and N = Ø. If product graph only contains 
    ## anchor nodes (|anchor| >= 2), then N = Ø.
    if not mod_product_graph.nodes or anchor_nodes == mod_product_graph.nodes:
        return [edge_anchor]

    ## Mapping edges in the product graph to their color
    color_dictionary = nx.get_edge_attributes(mod_product_graph, "color")

    ## computing N by intersecting all neighbours among the anchor points (NOTE: possible optimization when computing the intersection)
    common_neighbours_N = list(mod_product_graph.adj[anchor_nodes[0]].keys())
    for nodes in anchor_nodes:
        common_neighbours_N = [value for value in common_neighbours_N if value in mod_product_graph.adj[nodes].keys()]

    listN = blue_component_filter(mod_product_graph, common_neighbours_N, anchor_nodes, color_dictionary)

    ## If no components exist, the anchor is the MCS
    if len(listN) == 0:
        return [edge_anchor]
    else:
        MCSs = connected_MCS(listN, mod_product_graph, anchor_nodes, color_dictionary)

        all_mappings = []
        ## list of clique extensions (i.e. list of list of tuples)
        for mappings in MCSs:
            current_mapping = []
            for tuples in mappings:
                ## (i.e. one tuple of the form (a, b, c, d, ..., z)) up to the number of graphs
                ## tranformed into their corresponding edge in the graph
                mapped_edges = [edge_lists[i][tuples[i]] for i in range(n_graphs)]
                current_mapping.append(mapped_edges)
            all_mappings.append(current_mapping)
        return all_mappings

def iterative_approach(L, edge_anchor, limit_pg=True, molecule=False):
    """
        Computes the maximum common subgraph of all graphs in L w.r.t the anchors in edge_anchor.
        
        `Parameters`
            L (list: Graph): List of NetworkX graphs. Graphs may be decorated with labels.

            edge_anchor (list( list(Edge))): A list of edge anchors. Each edge anchor is a list of edges,
                                            each edge belongs to the graph in the order they are listed.
                                            The order is always edges from L[0], L[1] and so on.
                                            
        `Optional`
            limit_pg: Indicates whether the product graph should be limited to the neighbourhood of anchors or not, default to true.

            molecule: Indicates whether the graphs in L are decorated molecules. If true, it is expected that each graph has
                      attribute "atom_type" on nodes and "bond_type" on edges. This further limits the tuples in the product graph.
                      Default to false.
    """

    def create_induced_graph(mapping, graph_extract_attributes):
        """
            
            Creates the edge induced subgraph
            
            The attributes are inhertied from graph_extract_attributes (which should include all edges in mapping)
        """
        induced_graph = nx.Graph()
        ## Create edge-induced graph based on the given mapping
        ## circumvent the issue of out-of-order edges
        nodes_to_add = set()
        for edge_list in mapping:
            (u, v) = edge_list[0]
            nodes_to_add.add(u)
            nodes_to_add.add(v)
        nodes_to_add = sorted(list(nodes_to_add))
        induced_graph.add_nodes_from(nodes_to_add)
        induced_graph.add_edges_from(sorted([edge_list[0] for edge_list in mapping]))
        if molecule:
            atom_types = nx.get_node_attributes(graph_extract_attributes, "atom_type")
            bond_types = nx.get_edge_attributes(graph_extract_attributes, "bond_type")

            nx.set_node_attributes(induced_graph, atom_types, "atom_type")
            nx.set_edge_attributes(induced_graph, bond_types, "bond_type")
        
        return induced_graph

    def find_unique_graphs(all_mappings, graph_to_induce):
        """
            all_mappings: a list containing the mappings returned from mcs_list_leviBarrowBurstall
            graph_to_induce: The graph to induce based on each mapping
            
            Filters out isomorphic graphs, as it is unnecessary to consider all branches 
            of isomorphic graphs, only needs to look at one. 

            Returns :  
                unique_graphs: a list containing all unique graphs
                unique_mappings: a dict ontaining the mappings corresponding to each unique graph in unique_graphs
        """
        index_counter = 0
        ## Saving graphs that are not isomporphic to already seen graphs
        unique_graphs = []
        ## mapping each index in unique_mappings to their mapping
        unique_mappings = {}
        for mappings in all_mappings:
            
            ## Creates the induced graph from the current mapping
            found_subgraph = create_induced_graph(mappings, graph_to_induce)

            found_isomorph = False
            ## The found_subgraph needs to be checked against each graph already added as a unique graph
            for graph in unique_graphs:
                if molecule:
                    node_match = iso.categorical_node_match("atom_type", "")
                    edge_match = iso.categorical_edge_match("bond_type", "")
                    if nx.is_isomorphic(found_subgraph, graph, node_match, edge_match):
                        found_isomorph = True
                        break
                else:
                    if nx.is_isomorphic(found_subgraph, graph):
                        found_isomorph = True
                        break
                
            ## Only add graphs to the list if it is not isomorphic to any existing graphs
            if not found_isomorph:
                unique_graphs.append(found_subgraph)
                unique_mappings[index_counter] = mappings
                index_counter += 1
        
        return unique_graphs, unique_mappings

    def _iterative_approach_rec(L, current_mcs_graph, to_mcs_graph, all_mappings, current_mapping, anchor_bound, anchor, graph_amt, limit_pg=True, molecule=False):
        ## If end of L is reached, add the current mapping to the global list of mappings
        if to_mcs_graph == graph_amt:
            all_mappings.append(current_mapping)
            return

        graph_one = current_mcs_graph
        graph_two = L[to_mcs_graph]

        ## Map edges from current best graph to the upcoming "to_mcs_graph"
        new_anchor = [ [lists[0], lists[to_mcs_graph] ] for lists in anchor ]

        mcs = mcs_list_leviBarrowBurstall([graph_one, graph_two], new_anchor, limit_pg, molecule)
        
        ## Filter duplicates, no need to branch multiple times for identical mappings
        filtered_mcs = []
        for l in mcs:
            if sorted(l) not in filtered_mcs:
                filtered_mcs.append(sorted(l))
        
        unique_graphs, unique_mappings = find_unique_graphs(filtered_mcs, graph_one)  
        
        for i in range(len(unique_graphs)):
            graph_to_recurse = unique_graphs[i]
            mapping_to_recurse = unique_mappings[i]
            ## Only add extensions of the anchor
            if len(mapping_to_recurse) > anchor_bound:  
                
                ## If no mapping has been created yet, this is the first mapping
                if not current_mapping:
                    continue_mapping = mapping_to_recurse
                
                ## otherwise, update the current mapping with the new edge in this graph
                else:
                    new_current_mapping = copy.deepcopy(current_mapping)
                    continue_mapping = []
                    for i in range(len(mapping_to_recurse)):
                        edge_list = mapping_to_recurse[i]
                        ## Update the current mapping to include this found mapping, and update this mapping to track previous mappings including itself
                        for j in range(len(current_mapping)):
                            current_map = current_mapping[j]
                            ## update those lists that include already mapped edges
                            if edge_list[0] == current_map[0]: 
                                ## add newly mapped member to old list of mapped members
                                new_current_mapping[j].append(edge_list[1])
                                ## only move those mapped edges forward with newly mapped edges
                                continue_mapping.append(new_current_mapping[j])
                ## Continue recursively
                _iterative_approach_rec(L, graph_to_recurse, to_mcs_graph + 1, all_mappings, continue_mapping, anchor_bound, anchor, graph_amt, limit_pg, molecule)

    ## Use anchor_size as guard in the recursive step, terminating branches that reach this length
    anchor_size = len(edge_anchor)
    mapping_list = []

    ## first recursive step is between graph 0 and graph 1. 
    ## the mapping list is updated for each recursive call that ends up with 
    ## an actual extension of the anchor.
    _iterative_approach_rec(L, L[0], 1, mapping_list, [], anchor_size, edge_anchor, len(L), limit_pg, molecule)

    if not mapping_list:
        mapping_list = [edge_anchor]

    return mapping_list

def all_products(L, edge_anchor, limit_pg=True, molecule=False):
    """
        L: List of graphs
    """

    subgraphs = mcs_list_leviBarrowBurstall(L, edge_anchor, limit_pg, molecule)

    return subgraphs


def test_graphs(Gs, As, seq, molecules=False):

    graph_seq = [Gs[i] for i in seq]
    anchor_seq = [As[i] for i in seq]

    test_anchors = graph_format.compute_anchor(graph_seq, anchor_seq, molecule=molecules)
    print(f"Number of computed anchors: {len(test_anchors)}")

    global_maximum = 0
    for anchor in test_anchors:
        time_before = time.time()
        res_iterative = iterative_approach(graph_seq, anchor, molecule=True)
        time_after = time.time() 
        
        map_lengths = max([len(i) for i in res_iterative])
        if map_lengths > global_maximum: global_maximum = map_lengths

        max_mapping = list(filter(lambda x: len(x) == map_lengths, res_iterative))

        # if map_lengths >= 8:
        #     draw_molecules(graph_seq, [max_mapping[0]], anchor)

        print(f"Max extension: {map_lengths}")
        print(f"Number of extensions: {len(res_iterative)}")
        print(f"Number of extensions of max size: {len(max_mapping)}")
        print(f"time spent: {time_after-time_before} seconds")
        print()
    
    print(global_maximum)


# if __name__ == "__main__":          
#     graphs, anchors = graph_format.convert_graph_file("../labelled_graphs/alcohol_dehydrogenase_ethanol_backward.txt")
#     dist_map, shortest_distance = anchor_reach(graphs, anchors)

#     shrinked_graphs = shrink_graphs(graphs, 7, dist_map)

#     test_graphs(graphs, anchors, [0, 2, 1], True)

    
        
        

        