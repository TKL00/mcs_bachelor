import networkx as nx
import itertools

def switch_bond_type(input):
    """Switch function to convert the bond type from given data to our format
    Args:
        input (String): String specifying the bond type from the input
    Returns:
        String: specifying the bond type in our format
    """
    if input == "-":
        return "s"
    elif input == "=":
        return "d"
    elif input == ":":
        return "a"
    elif input == "==":
        return "t"
    elif input == "-=":
        return "s/d"
    elif input == "=-":
        return "d/s"
    else:
        return "q"

def convert_graph_file(path):
    

    ## The file containing our data is currently hardcodedd, could potentially be changed
    with open(path, "r") as file:
        
        ## Reading the first line from the file to get past the first "---New Instance---" line
        file.readline()
        reading_nodes = True
        
        ## Lists to save the graphs and the anchor edges
        all_graphs = []
        all_anchors = []
        
        ## Creating a new graph and dictionaries to contain the attributes, and a list to save the anchors
        G = nx.Graph()
        node_attribute_dict = {}
        edge_attribute_dict = {}
        anchor_edges = []
        
        for string in file:
            ## If the string is "New Instance" the current graph is completed and needs to be saved
            ## and preparing to create next graph
            if(string == "---New Instance---\n"):
                ## Save current graph
                nx.set_node_attributes(G, node_attribute_dict)
                nx.set_edge_attributes(G, edge_attribute_dict)
                
                all_graphs.append(G)
                all_anchors.append(anchor_edges)
                
                
                G = nx.Graph()
                node_attribute_dict = {}
                edge_attribute_dict = {}
                anchor_edges = []
                reading_nodes = True
            ## If string "###" is read, the nodes for the current graph has been read and the edges comes after.
            elif(string == "###\n"):
                reading_nodes = False
            
            else: 
                ## Adding nodes
                if reading_nodes:
                    string_split = string.split(" ")
                    node = int(string_split[0])
                    G.add_node(node)
                    node_attribute_dict[node] = {"atom_type": string_split[1].strip()}
                ## Adding edges
                else:
                    string_split = string.split(" ")
                    edge = (int(string_split[0]), int(string_split[1]))
                    G.add_edge(int(string_split[0]), int(string_split[1]))
                    ## Checking if edge is an anchor edge
                    if(string_split[2] == "anchor"):
                        if(len(string_split) == 8):
                            bond_1 = string_split[4].strip()
                            bond_2 = string_split[6].strip()    
                            bond = bond_1 + bond_2
                        else:
                            if(string_split[4] != ","):
                                bond = string_split[4].strip()
                            else:
                                bond = string_split[5].strip()
                        anchor_edges.append(edge)
                    else:
                        bond = string_split[2].strip()
                    ## Using switch function to convert to our format of bond type
                    correct_bond = switch_bond_type(bond)
                    edge_attribute_dict[edge] = {"bond_type": correct_bond}


    #Saving the last graph
    nx.set_node_attributes(G, node_attribute_dict)
    nx.set_edge_attributes(G, edge_attribute_dict)
    all_graphs.append(G)
    all_anchors.append(anchor_edges)
    
    return all_graphs, all_anchors

def compute_anchor(Gs, AEs, molecule=False):
    """
        Gs: List of graphs
        AEs: List of anchored edges (list of edges). AEs[i] -> anchored edges in Gs[i]

        Returns: A list of anchors (list of lists of edges mapped together)        

    """

    computed_anchors = []
    n_graphs = len(Gs)
    n_anchored_edges = len(AEs[0])
    ## g_edge_types[i] is a dictionary from edge types in anchor to
    ## the anchor edges of this type in g_i
    g_edge_types = []

    ## combinations are limited based on edge type
    if molecule:
        ## Compute g_edge_types
        for i in range(n_graphs):
            g = Gs[i]
            anchors = AEs[i]
            ## Mapping edge types to list of edges
            g_edge_type_map = {}
            g_atom_type = nx.get_node_attributes(g, "atom_type")
            g_bond_type = nx.get_edge_attributes(g, "bond_type") 
            for j in range(n_anchored_edges):
                ## An anchor edge_j in G_i
                (u, v) = anchors[j]
                atom_pair = tuple(sorted((g_atom_type[u], g_atom_type[v])))
                ## Ignore networkX edge ordering problems
                try:
                    bond_type = g_bond_type[(u, v)]
                except:
                    bond_type = g_bond_type[(v, u)]
                
                ## edges are identified by unique atom_pair and bond type
                edge_type = (atom_pair, bond_type)
                ## If edge type not discovered yet, add the key and value
                if edge_type not in g_edge_type_map:
                    g_edge_type_map[edge_type] = [(u, v)]
                ## otherwise add to list
                else:
                    g_edge_type_map[edge_type].append((u, v))
            g_edge_types.append(g_edge_type_map)
        
        ## g_edge_types now contains a list of dictionaries, one for each g in Gs

        ## Save the first dictionary for graph 0 to retrieve attributes
        init_edge_type_dict = g_edge_types[0]

        ## List of all possible edge types e.g. (('O', 'P'), 's'), (('H', 'O), 'd')
        possible_edge_types = [key for key in init_edge_type_dict]
        ## Number of edges of the different types
        n_edge_type_edges = {key: len(init_edge_type_dict[key]) for key in init_edge_type_dict}
        n_edge_types = len(possible_edge_types)
        edge_type_options = {}

        ## Compute all combinations of edges of the same edge type across all graphs
        for edge_type in possible_edge_types:
            n_edges_of_type = n_edge_type_edges[edge_type]          
            indices = [i for i in range(n_edges_of_type)]
            
            permutations_of_indices = []
            for j in range(n_graphs):
                permutations_of_indices.append(list(itertools.permutations(indices)))
            
            ## list of tuples, each tuple being a mapping of an edge type ((x, y), (a, b), (u, v)) where x -> a -> u and y -> b -> v, no one is mapped to the same node.
            mappings_of_edge_type = list(itertools.product(*permutations_of_indices))
            edge_type_options[edge_type] = mappings_of_edge_type
            
        all_edge_type_map_combinations = list(itertools.product(*edge_type_options.values()))

        ## Each combination tuple corresponds to (x, y, ..) where x is a mapping of edge type 0 in all graphs, y is a mapping of all edges of type 1 and so on.
        ## That is, each tuple is an anchor option.
        for combinations in all_edge_type_map_combinations:
            new_anchor = []
            ## Each combination has tuples for each edge type in the order of edges in possible_edge_types
            for i in range(n_edge_types):
                current_edge_type = possible_edge_types[i]
                edge_type_tuple = combinations[i]

                ## Each contains the number of edges for the given edge type
                for j in range(n_edge_type_edges[current_edge_type]):
                    edges_mapped_together = []
                    
                    ## Each edge must, for all graphs, be mapped to an edge of that type in that graph
                    for k in range(n_graphs):
                        current_graph_tuple = edge_type_tuple[k]
                        g_edge = g_edge_types[k][current_edge_type][current_graph_tuple[j]]
                        edges_mapped_together.append(g_edge)
                    
                    new_anchor.append(edges_mapped_together)
                    
            computed_anchors.append(new_anchor)
        
        return computed_anchors
        
    ## all combinations
    else:
        return []
