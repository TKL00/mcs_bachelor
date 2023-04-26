import networkx as nx

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
    else:
        return "q"

## The file containing our data is currently hardcodedd, could potentially be changed
with open("../out.txt", "r") as file:
    
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
            print(f"Inside the new instance")
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
            print("in ###")
        
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
                #print(f"adding edge: {edge}")
                G.add_edge(int(string_split[0]), int(string_split[1]))
                ## Checking if edge is an anchor edge
                if(string_split[2] == "anchor"):
                    ## Finding out if the bond type is for 
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

print(nx.get_node_attributes(all_graphs[0], "atom_type"))