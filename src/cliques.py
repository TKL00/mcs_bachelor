from productgraph import product_graph_no_limit as pgnl
from productgraph import product_graph_limit as pgl
from linegraph import line_graph as lg
from linegraph import convert_edge_anchor_lg_list
from itertools import chain
import graph_format
import networkx as nx
from queue import Queue
import copy
from draw_graphs import draw_product_graph, draw_two_graphs, draw_blue_connected_components, draw_molecules
from molecules import *
import time


### Authors: Tobias Klink Lehn (toleh20@student.sdu.dk) and Kasper Halkjær Beider (kbeid20@student.sdu.dk)
def mcs_list_leviBarrowBurstall(L, edge_anchor, limit_pg=True, molecule=False):
    """
        Computes the Maximum Common Subgraph using the Algorithm suggested by
        G. Levi and H.G. Barrow + R.M. Burstall in 1973 and 1975 respectively.
        
        `Paramters`:
            G (Graph): A NetworkX graph, nodes are integers but may be decorated with items

            H (Graph): A NetworkX graph, nodes are integers but may be decorated with items

            edge_anchor (list: list(edge)): A valid one-to-one mapping from 'n' edges in G to 'n' edges in H, overrules a given node_anchor

        `Optional`:
            limit_pg (boolean): Indicates whether the product graph should be limited to the neighbourhood of anchors or not, default to true

            molecule (boolean): Indicates whether the graphs in L are decorated molecules. If true, it is expected that each graph has
                                attribute "atom_type" on nodes and "bond_type" on edges. 
        
        `Returns`:
            all_mappings (list( dict: node -> node) or list( dict: edge -> edge)): 
            
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
        L: List of graphs
    """

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

        for found_mapping in mcs:
            ## Only add extensions of the anchor
            if len(found_mapping) > anchor_bound:  
                ## Create copy for each mapping to recurse on
                new_current_mapping = copy.deepcopy(current_mapping) 

                ## Create induced graph to potentially to recurse on
                current_mcs_graph = nx.Graph()   
                nodes_to_add = set()
                for edge_list in found_mapping:
                    (u, v) = edge_list[0]
                    nodes_to_add.add(u)
                    nodes_to_add.add(v)
                nodes_to_add = sorted(list(nodes_to_add))
                current_mcs_graph.add_nodes_from(nodes_to_add)    
                current_mcs_graph.add_edges_from(sorted([edge_list[0] for edge_list in found_mapping]))

                # Inherit labels from current mcs
                if molecule:
                    atom_types = nx.get_node_attributes(graph_one, "atom_type")
                    bond_types = nx.get_edge_attributes(graph_one, "bond_type")
                    nx.set_node_attributes(current_mcs_graph, atom_types, "atom_type")
                    nx.set_edge_attributes(current_mcs_graph, bond_types, "bond_type")
                ## Update the found mapping with the previously found mapped edges

                continue_mapping = []
                for i in range(len(found_mapping)):
                    edge_list = found_mapping[i]
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
                _iterative_approach_rec(L, current_mcs_graph, to_mcs_graph + 1, all_mappings, continue_mapping, anchor_bound, anchor, graph_amt, limit_pg, molecule)

    ## Use anchor_size as guard in the recursive step, terminating branches that reach this length
    anchor_size = len(edge_anchor)
    mapping_list = []
    ## Map L[0] edges to L[1] edges as initial anchor
    start_anchor = [ [lists[0], lists[1]] for lists in edge_anchor]
    graph_one = L[0]
    graph_two = L[1]

    mcs = mcs_list_leviBarrowBurstall([graph_one, graph_two], start_anchor, limit_pg, molecule)
    for mappings in mcs:
        current_mcs_graph = nx.Graph()
        ## Create edge-induced graph based on the given mapping
        ## circumvent the issue of out-of-order edges
        nodes_to_add = set()
        for edge_list in mappings:
            (u, v) = edge_list[0]
            nodes_to_add.add(u)
            nodes_to_add.add(v)
        nodes_to_add = sorted(list(nodes_to_add))
        current_mcs_graph.add_nodes_from(nodes_to_add)
        current_mcs_graph.add_edges_from(sorted([edge_list[0] for edge_list in mappings]))
        if molecule:
            atom_types = nx.get_node_attributes(graph_one, "atom_type")
            bond_types = nx.get_edge_attributes(graph_one, "bond_type")

            nx.set_node_attributes(current_mcs_graph, atom_types, "atom_type")
            nx.set_edge_attributes(current_mcs_graph, bond_types, "bond_type")

        _iterative_approach_rec(L, current_mcs_graph, 2, mapping_list, copy.deepcopy(mappings), anchor_size, edge_anchor, len(L), limit_pg, molecule)

    ## If nothing was added to the mapping list along the way, the anchor is the only common subgraph
    if not mapping_list:
        mapping_list = [edge_anchor]

    return mapping_list

def all_products(L, edge_anchor, limit_pg=True, molecule=False):
    """
        L: List of graphs
    """

    subgraphs = mcs_list_leviBarrowBurstall(L, edge_anchor, limit_pg, molecule)

    return subgraphs

if __name__ == "__main__":

    graphs, anchors = graph_format.convert_graph_file("../out.txt")
    all_anchors = graph_format.compute_anchor(graphs, anchors, molecule=True)

    graph_one_atoms = nx.get_node_attributes(graphs[0], "atom_type")
    graph_one_bonds = nx.get_edge_attributes(graphs[0], "bond_type")

    graph_two_atoms = nx.get_node_attributes(graphs[1], "atom_type")
    graph_two_bonds = nx.get_edge_attributes(graphs[1], "bond_type")
    
    # print("graph one")
    # for i in graph_one_atoms: print(i, graph_one_atoms[i])
    # for i in graph_one_bonds: print(i, graph_one_bonds[i])
    # print("graph two")
    # for i in graph_two_atoms: print(i, graph_two_atoms[i])
    # for i in graph_two_bonds: print(i, graph_two_bonds[i])
    trial_graphs = [graphs[1], graphs[2], graphs[0], graphs[3]]
    trial_anchors = [anchors[1], anchors[2], anchors[0], anchors[3]]
    # fail_anchor = [[(53, 54), (10, 14), (0, 1), (5, 14)], [(54, 123), (5, 14), (1, 10), (6, 14)], [(123, 124), (6, 10), (10, 11), (6, 7)], [(53, 124), (5, 6), (0, 11), (5, 7)]]
    trial_all_anchors = graph_format.compute_anchor(trial_graphs, trial_anchors, molecule=True)

    # iterative_approach(trial_graphs, edge_anchor=fail_anchor, molecule=True)

    for trial_anchor in trial_all_anchors:
        print(f"anchor {trial_anchor}")
        time_before = time.time()
        res_iterative = iterative_approach(trial_graphs, edge_anchor=trial_anchor, molecule=True)
        time_after = time.time()
        print(res_iterative)
        print(f"Time taken: {time_after - time_before}")
        
        
    print(len(all_anchors))