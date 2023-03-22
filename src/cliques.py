from productgraph import product_graph_no_limit as pgnl
from productgraph import product_graph_limit as pgl
from linegraph import line_graph as lg
from linegraph import convert_edge_anchor_lg_list
import networkx as nx
from queue import Queue
import copy
from draw_graphs import draw_product_graph, draw_two_graphs, draw_blue_connected_components, draw_molecules


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

            edge_anchor (dict: int -> [int]): A valid one-to-one mapping from 'n' edges in G to 'n' edges in H, overrules a given node_anchor
        
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

                all_components (list( list(node) )): A list of blue connected components all reachable by blue edges from the anchor, A.
        """
        ## BFS to find all nodes that are connected to A through blue edges
        ## all nodes in the product graph are white in the beginning
        ## NOTE: optimization for not whitening all the nodes in PG, but only the ones in A and N.
        
        ## Avoiding all nodes not in A and N by giving them a scary color
        node_color_lookup = { node: "p" for node in PG.nodes}
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
            for v in PG.adj[u]:
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
                for v in PG.adj[u]:
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

        for component in listN:
            ## create vertex induced subgraph
            comp_subgraph = nx.Graph(nx.induced_subgraph(PG, component))
            ## compute L
            cliques = list(nx.find_cliques(comp_subgraph))
            for clique in cliques:
                blue_edge_found = False
                for node in clique:
                    neighbs = PG.adj[node]
                    for neighbour in neighbs:
                        if neighbour in A:
                            ## determine edge color between the node and its anchor neighbour (if one exists)
                            edge_color = edge_color_lookup[(neighbour, node)] if (neighbour, node) in edge_color_lookup else edge_color_lookup[(node, neighbour)]
                            if edge_color == "blue":
                                ## Upon blue edge discovery, add entire clique and anchor as a MCS (clique is extension of A)
                                MCSs.append(clique + A)
                                blue_edge_found = True
                                break
                    ## pursue next clique if this one already has a blue edge
                    if blue_edge_found: break
        return MCSs

    G_edges = list(G.edges)
    H_edges = list(H.edges)

    ## Use Line-Graph implementation
    if edge_anchor:
        LG = lg(G)
        LH = lg(H)

        ## Mapping: v in LG -> u in LH
        computed_node_anchor = convert_edge_anchor_lg_list([G, H], edge_anchor)
        
        ## Compute product graph.
        mod_product_graph = pgl([LG, LH], computed_node_anchor)
        ## Mapping edges in the "LG x LH" to their color
        color_dictionary = nx.get_edge_attributes(mod_product_graph, "color")

        ## Anchor points in product graph correspond to (x, y) nodes where 
        ## x is a node in LG and y is a node in LH. Adjusted
        anchor_points = [ (g, computed_node_anchor[g][0]) for g in computed_node_anchor ]

        ## computing N by intersecting all neighbours among the anchor points (NOTE: possible optimization when computing the intersection)
        common_neighbours_N = list(mod_product_graph.adj[anchor_points[0]].keys())
        for nodes in anchor_points:
            common_neighbours_N = [value for value in common_neighbours_N if value in mod_product_graph.adj[nodes].keys()]
            
        listN = blue_component_filter(mod_product_graph, common_neighbours_N, anchor_points, color_dictionary)

        ## If no components exist, the anchor is the MCS
        if len(listN) == 0:
            return edge_anchor
        else:
            MCSs = connected_MCS(listN, mod_product_graph, anchor_points, color_dictionary)
            all_mappings = []
            for mappings in MCSs:
                current_mapping = {}
                for tuples in mappings:
                    (u, v) = tuples
                    current_mapping[G_edges[u]] = H_edges[v]
                all_mappings.append(current_mapping)
            return all_mappings
            # draw_blue_connected_components(mod_product_graph, anchor_points, listN, color_dictionary)

    ## NOTE: Ask for "permission" to clean this up such that we only focus on edge-anchored version.
    ## No use of line graphs
    else:
        mod_product_graph = pgnl(G, H)

        anchor_points = [ (g, node_anchor[g]) for g in node_anchor]

                                                ### FIND CLIQUES CONTAINING ANCHOR AND RETURN MAP
        ## Finds all cliques containing the anchor points
        anchored_cliques = list(nx.find_cliques(mod_product_graph, anchor_points))
        
        ## In each clique, a node pair (i, j) corresponds to the i'th node
        ## in G being mapped to the j'th node in H. 
        all_mappings = []
        for cliques in anchored_cliques:
            new_mapping = {}
            for points in cliques:
                new_mapping[points[0]] = points[1]
            all_mappings.append(new_mapping)
        
        return all_mappings

def mcs_list_leviBarrowBurstall(L, edge_anchor, limit_pg=True, molecule=False):
    """
        Computes the Maximum Common Subgraph using the Algorithm suggested by
        G. Levi and H.G. Barrow + R.M. Burstall in 1973 and 1975 respectively.
        
        `Paramters`:
            G (Graph): A NetworkX graph, nodes are integers but may be decorated with items

            H (Graph): A NetworkX graph, nodes are integers but may be decorated with items

            edge_anchor (dict: int -> int): A valid one-to-one mapping from 'n' edges in G to 'n' edges in H, overrules a given node_anchor

        `Optional`:
            limit_pg (boolean): A boolean to indicate whether the product graph should be limited to the neighbourhood of anchors or not, default to true
        
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
            for v in PG.adj[u]:
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
                for v in PG.adj[u]:
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

        for component in listN:
            ## create vertex induced subgraph
            comp_subgraph = nx.Graph(nx.induced_subgraph(PG, component))
            ## compute L
            cliques = list(nx.find_cliques(comp_subgraph))
            for clique in cliques:
                blue_edge_found = False
                for node in clique:
                    neighbs = PG.adj[node]
                    for neighbour in neighbs:
                        if neighbour in A:
                            ## determine edge color between the node and its anchor neighbour (if one exists)
                            edge_color = edge_color_lookup[(neighbour, node)] if (neighbour, node) in edge_color_lookup else edge_color_lookup[(node, neighbour)]
                            if edge_color == "blue":
                                ## Upon blue edge discovery, add entire clique and anchor as a MCS (clique is extension of A)
                                MCSs.append(clique + A)
                                blue_edge_found = True
                                break
                    ## pursue next clique if this one already has a blue edge
                    if blue_edge_found: break
        return MCSs

    n_graphs = len(L)
    edge_lists = [list(L[i].edges) for i in range(n_graphs)]

    ## Use Line-Graph implementation

    linegraphs = [lg(L[i], molecule=molecule) for i in range(n_graphs)]

    ## Mapping: v in L[0] -> u in [L[1] for i > 0]
    computed_node_anchor = convert_edge_anchor_lg_list(L, edge_anchor)
    
    ## Compute product graph, either constraining it to only include A and N, or include all possible nodes.
    mod_product_graph = pgl(linegraphs, computed_node_anchor, molecule=molecule) if limit_pg else pgnl(linegraphs)

    ## Mapping edges in the product graph to their color
    color_dictionary = nx.get_edge_attributes(mod_product_graph, "color")

    ## Anchor points in product graph correspond to (v_1, v_2, ..., v_n) nodes where 
    ## v_0 is a node in L[0] and v_i is a node in L[i]. Folding out the list of mappings.
    anchor_points = [ (l_zero, *computed_node_anchor[l_zero]) for l_zero in computed_node_anchor ]

    ## computing N by intersecting all neighbours among the anchor points (NOTE: possible optimization when computing the intersection,)
    common_neighbours_N = list(mod_product_graph.adj[anchor_points[0]].keys())
    for nodes in anchor_points:
        common_neighbours_N = [value for value in common_neighbours_N if value in mod_product_graph.adj[nodes].keys()]

    listN = blue_component_filter(mod_product_graph, common_neighbours_N, anchor_points, color_dictionary)
    ## If no components exist, the anchor is the MCS
    if len(listN) == 0:
        return edge_anchor
    else:
        MCSs = connected_MCS(listN, mod_product_graph, anchor_points, color_dictionary)

        all_mappings = []
        ## list of clique extensions
        for mappings in MCSs:
            current_mapping = {}
            for tuples in mappings:
                L_0 = edge_lists[0][tuples[0]]
                mapped_edges = []
                for i in range(1, n_graphs):
                    mapped_edges.append(edge_lists[i][tuples[i]])
                ## v -> [x_2, x_3, ..., x_n]
                current_mapping[L_0] = mapped_edges
            all_mappings.append(current_mapping)
        return all_mappings


def iterative_approach(L, anchored_edges):
    """
        L: List of graphs
    """

    return None

def all_products(L, anchored_edges):
    """
        L: List of graphs
    """

    subgraphs = mcs_list_leviBarrowBurstall(L, anchored_edges)

    return None


if __name__ == "__main__":
    G = nx.Graph()
    G.add_edges_from([(0,1), (0,2), (1,2), (2,3), (3,4)])
    G_node_attributes = {
                         0: {"atom_type":"C"},
                         1: {"atom_type":"N"},
                         2: {"atom_type":"C"},
                         3: {"atom_type":"O"},
                         4: {"atom_type":"C"}
                        }
    G_edge_attributes = {
                        (0, 1): {"bond_type": "s"},
                        (0, 2): {"bond_type": "s"},
                        (1, 2): {"bond_type": "s"},
                        (2, 3): {"bond_type": "s"},
                        (3, 4): {"bond_type": "s"},
                        }
    nx.set_edge_attributes(G, G_edge_attributes)
    nx.set_node_attributes(G, G_node_attributes)
    
    lG = lg(G, molecule=True)

    H = nx.Graph()
    H.add_edges_from([(0,1), (1,2), (1,3), (1,7), (2,3), (3,4), (3,6), (6, 7), (7, 8), (4,5)])
    H_node_attributes = {0: {"atom_type":"C"},
                         1: {"atom_type":"C"},
                         2: {"atom_type":"C"},
                         3: {"atom_type":"C"},
                         4: {"atom_type":"O"},
                         5: {"atom_type":"C"},
                         6: {"atom_type":"C"},
                         7: {"atom_type":"C"},
                         8: {"atom_type":"C"},
                        }
    H_edge_attributes = {
                        (0, 1): {"bond_type": "s"},
                        (1, 2): {"bond_type": "s"},
                        (1, 3): {"bond_type": "s"},
                        (1, 7): {"bond_type": "s"},
                        (2, 3): {"bond_type": "s"},
                        (3, 4): {"bond_type": "s"},
                        (3, 6): {"bond_type": "s"},
                        (6, 7): {"bond_type": "s"},
                        (7, 8): {"bond_type": "s"},
                        (4, 5): {"bond_type": "s"},
                        }
    nx.set_edge_attributes(H, H_edge_attributes)
    nx.set_node_attributes(H, H_node_attributes)
    lH = lg(H, molecule=True)

    I = nx.Graph()
    I.add_edges_from([(0,1),(0,3),(1, 2),(2, 3)])
    I_node_attributes = {
                         0: {"atom_type":"C"},
                         1: {"atom_type":"O"},
                         2: {"atom_type":"C"},
                         3: {"atom_type":"C"},
                        }
    I_edge_attributes = {
                        (0, 1): {"bond_type": "s"},
                        (0, 3): {"bond_type": "s"},
                        (1, 2): {"bond_type": "s"},
                        (2, 3): {"bond_type": "s"},
                        }
    nx.set_edge_attributes(I, I_edge_attributes)
    nx.set_node_attributes(I, I_node_attributes)

    lI = lg(I, molecule=True)

    lG_atom_attributes = nx.get_node_attributes(lG, "atom_pair")
    lH_atom_attributes = nx.get_node_attributes(lH, "atom_pair")
    lI_atom_attributes = nx.get_node_attributes(lI, "atom_pair")

    lG_bond_attributes = nx.get_node_attributes(lG, "bond_type")
    lH_bond_attributes = nx.get_node_attributes(lH, "bond_type")
    lI_bond_attributes = nx.get_node_attributes(lI, "bond_type")

    list_edge_anchor = {
        (3, 4): [(4, 5), (0, 1)],
        (2, 3): [(3, 4), (1, 2)]
    }

    pg = pgl([lG, lH, lI], lg_node_anchor=convert_edge_anchor_lg_list([G, H, I], list_edge_anchor), molecule=True)

    # subgraphs_list_three = mcs_list_leviBarrowBurstall([G, H, I], list_edge_anchor)
    # print(subgraphs_list_three)
          
    # subgraphs_list_three_molecules = mcs_list_leviBarrowBurstall([G, H, I], list_edge_anchor, molecule=True)
    # print(subgraphs_list_three_molecules)


    propanic_acid = nx.Graph()
    propanic_acid.add_edges_from([(0, 1), (1, 2), (2, 3), (2, 4), (4, 5), (0, 8), (0, 9), (0, 10), (1, 6), (1, 7)])
    propanic_acid_node_attributes = {
                         0: {"atom_type":"C"},
                         1: {"atom_type":"C"},
                         2: {"atom_type":"C"},
                         3: {"atom_type":"O"},
                         4: {"atom_type":"O"},
                         5: {"atom_type": "H"},
                         6: {"atom_type": "H"},
                         7: {"atom_type":"C"},
                         8: {"atom_type": "H"},
                         9: {"atom_type": "H"},
                         10: {"atom_type": "H"}
                        }
    propanic_acid_edge_attributes = {
                        (0, 1): {"bond_type": "s"}, 
                        (1, 2): {"bond_type": "s"},
                        (2, 3): {"bond_type": "d"}, 
                        (2, 4): {"bond_type": "s"}, 
                        (4, 5): {"bond_type": "s"}, 
                        (0, 8): {"bond_type": "s"}, 
                        (0, 9): {"bond_type": "s"}, 
                        (0, 10): {"bond_type": "s"}, 
                        (1, 6): {"bond_type": "s"}, 
                        (1, 7): {"bond_type": "s"}
                        }

    nx.set_node_attributes(propanic_acid, propanic_acid_node_attributes)
    nx.set_edge_attributes(propanic_acid, propanic_acid_edge_attributes)
    

    methane_acid = nx.Graph()
    methane_acid.add_edges_from([(0, 1), (0, 2), (0, 4), (2, 3)])
    methane_acid_node_attributes = {
                        0: {"atom_type":"C"},
                        1: {"atom_type":"O"},
                        2: {"atom_type":"O"},
                        3: {"atom_type":"H"},
                        4: {"atom_type":"H"},
    }
    methane_acid_edge_attributes = {
                        (0, 1): {"bond_type": "d"},
                        (0, 2): {"bond_type": "s"}, 
                        (0, 4): {"bond_type": "s"}, 
                        (2, 3): {"bond_type": "s"}
    }
    nx.set_node_attributes(methane_acid, methane_acid_node_attributes)
    nx.set_edge_attributes(methane_acid, methane_acid_edge_attributes)

    methanol = nx.Graph()
    methanol.add_edges_from([(0, 1), (0, 3), (0, 4), (0,5), (1, 2)])
    methanol_node_attributes = {
                         0: {"atom_type":"C"},
                         1: {"atom_type":"O"},
                         2: {"atom_type":"H"},
                         3: {"atom_type":"H"},
                         4: {"atom_type":"H"},
                         5: {"atom_type": "H"},
                        }
    methanol_edge_attributes = {
                        (0, 1): {"bond_type": "s"}, 
                        (0, 3): {"bond_type": "s"}, 
                        (0, 4): {"bond_type": "s"}, 
                        (0, 5): {"bond_type": "s"}, 
                        (1, 2): {"bond_type": "s"}
    }
    nx.set_node_attributes(methanol, methanol_node_attributes)
    nx.set_edge_attributes(methanol, methanol_edge_attributes)

    molecule_edge_anchor = {
        (2, 4): [(0, 2), (0, 1)]
    }

    molecule_subgraph = mcs_list_leviBarrowBurstall([propanic_acid, methane_acid, methanol], molecule_edge_anchor, molecule=True)
    
    draw_molecules([propanic_acid, methane_acid, methanol], molecule_subgraph,molecule_edge_anchor)
