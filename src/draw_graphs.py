import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import math
import copy

### Authors: Tobias Klink Lehn (toleh20@student.sdu.dk) and Kasper Halkjær Beider (kbeid20@student.sdu.dk)

GLOBAL_SEED = 1010203

def draw_product_graph(G, H, PGH):
    """
        Draws two NetworkX graphs next to each other, along with their
        associated product graph.
    """
    G_pos = nx.spring_layout(G, seed=101)
    H_pos = nx.spring_layout(H, seed=101)

    PGH_pos = nx.spring_layout(PGH, seed=101)
    G_node_labels = {}
    H_node_labels = {}
    PGH_node_labels = {}


    for nodes in G.nodes:
        G_node_labels[nodes] = nodes

    for nodes in H.nodes:
        H_node_labels[nodes] = nodes
    
    for nodes in PGH.nodes:
        PGH_node_labels[nodes] = nodes

    subax1 = plt.subplot(221)
    subax1.set_title("G")
    nx.draw_networkx_nodes(G, G_pos, nodelist=G.nodes)
    nx.draw_networkx_labels(G, G_pos, G_node_labels, font_size=15, font_color="whitesmoke")
    nx.draw_networkx_edges(G, G_pos, edgelist=G.edges)
    
    subax2 = plt.subplot(222)
    subax2.set_title("H")
    nx.draw_networkx_nodes(H, H_pos, nodelist=H.nodes)
    nx.draw_networkx_labels(H, H_pos, H_node_labels, font_size=15, font_color="whitesmoke")
    nx.draw_networkx_edges(H, H_pos, edgelist=H.edges)
    
    ## Uses one column in the second row
    subax3 = plt.subplot(212)
    subax3.set_title("Product Graph")
    nx.draw_networkx_nodes(PGH, PGH_pos, nodelist=PGH.nodes)
    nx.draw_networkx_labels(PGH, PGH_pos, PGH_node_labels, font_size=7, font_color="whitesmoke")
    nx.draw_networkx_edges(PGH, PGH_pos, edgelist=PGH.edges)


    plt.show()

def draw_graphs(L, mappings, edge_anchor):
    """
        Draws graphs in L in a grid with two columns and appropriate amount of rows.
        Highlights edge_anchor by grey edges, labels the edges based on the mappings.
    """

    label_string = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"

    n_graphs = len(L)
    n_rows = math.ceil(n_graphs / 2)

    for i in range(n_graphs):
        graph = L[i]

    ## Draw each mapping separately
    for mapping in mappings:
        fig, axs = plt.subplots(ncols=2, nrows=n_rows)
        graph_index = 0

        for ax in axs.flat:
            ## Remove border
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)
            if graph_index == n_graphs: break
            ax.set_title(f"Graph {graph_index}")

            
            graph = L[graph_index]
            node_labels = {i: i for i in graph.nodes}
            position = nx.spring_layout(graph, seed=101)

            ## Highlight anchored edges
            for lists in edge_anchor:
                edge_to_draw = [lists[graph_index]]
                nx.draw_networkx_edges(graph, position, edge_to_draw, width=8, edge_color="gray", alpha=0.3, ax=ax)
            
            ## Draw each edge individually, color is based on their bond type
            for edge in graph.edges:
                nx.draw_networkx_edges(graph, position, [edge], width=1.0, edge_color="black", ax=ax)
            
            ## Draw mapped edge labels
            label_iterator = 0
            for lists in mapping:
                label = {lists[graph_index]: label_string[label_iterator]}
                nx.draw_networkx_edge_labels(graph, position, label, ax=ax, font_color="lightseagreen", font_family="calibri", font_size=15)
                label_iterator += 1

            
            # ## Draw nodes and atom labels
            nx.draw_networkx_nodes(graph, position, graph.nodes, node_size=175, node_color="tab:blue", ax=ax)
            nx.draw_networkx_labels(graph, position, labels=node_labels, ax=ax)
            graph_index += 1

        ## Remove x- and y-axis
        plt.axis('off')
        plt.show()
    
def draw_blue_connected_components(PG, anchor, blue_components, color_map):

    PG_pos = nx.spring_layout(PG, seed=101)

    nodes_to_draw = copy.deepcopy(anchor)
    for component in blue_components:
        nodes_to_draw.extend(component)

    blue_component_graphs = []
    component_labels = {}
    for component in blue_components:
        nodes_to_add = copy.deepcopy(component)
        nodes_to_add.extend(anchor)
        subgraph = nx.Graph(nx.induced_subgraph(PG, nodes_to_add))
        blue_component_graphs.append(subgraph)
        for nodes in component:
            component_labels[nodes] = nodes

    anchor_labels = {point: point for point in anchor}
    
    subax1 = plt.subplot(111)

    ## draw each subgraph individually
    for blue_subgraphs in blue_component_graphs:

        nx.draw_networkx_nodes(blue_subgraphs, PG_pos, blue_subgraphs.nodes, node_color="tab:blue", node_size=400)
        nx.draw_networkx_labels(blue_subgraphs, PG_pos, component_labels, font_size=8)
        nx.draw_networkx_edges(PG, PG_pos, nodelist=blue_subgraphs.nodes, width=0.5, alpha=0.1)
        
        for edge in blue_subgraphs.edges:
            (u, v) = edge
            color = color_map[(u, v)] if (u, v) in color_map else color_map[(v, u)]
            nx.draw_networkx_edges(blue_subgraphs, PG_pos, edgelist=[edge], width=0.5, alpha=0.5, edge_color=color)


    for i in range(len(nodes_to_draw)):
        node_i = nodes_to_draw[i]
        for j in range(i + 1, len(nodes_to_draw)):
            node_j = nodes_to_draw[j]

            edge = (node_i, node_j) if (node_i, node_j) in color_map else (node_j, node_i)
            if edge in color_map:
                color = color_map[edge] if edge in color_map else color_map[edge]
                nx.draw_networkx_edges(PG, PG_pos, edgelist=[edge], width=0.5, alpha=0.5, edge_color=color)

    ## draw the anchor lastly
    nx.draw_networkx_nodes(PG, PG_pos, anchor, node_color="tab:orange", node_size=400)
    nx.draw_networkx_labels(PG, PG_pos, anchor_labels, font_size=8)

    plt.show()

    return None

def draw_molecules(L, mappings, edge_anchor):
    """
        Draws "molecule" graphs in L in a grid with two columns and appropriate amount of rows.
        Highlights edge_anchor by grey edges, labels the edges based on the mappings.
    """

    ## List of atom_type and bond_type maps for each graph
    atom_types = [nx.get_node_attributes(graph, "atom_type") for graph in L]
    bond_types = [nx.get_edge_attributes(graph, "bond_type") for graph in L]

    label_string = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"

    n_graphs = len(L)
    n_rows = math.ceil(n_graphs / 2)

    for i in range(n_graphs):
        graph = L[i]

    ## Draw each mapping separately
    for mapping in mappings:
        fig, axs = plt.subplots(ncols=2, nrows=n_rows)
        graph_index = 0

        for ax in axs.flat:
            ## Remove border
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)
            if graph_index == n_graphs: break
            ax.set_title(f"Graph {graph_index}")

            
            graph = L[graph_index]
            position = nx.spring_layout(graph, seed=101)

            ## Highlight anchored edges
            for lists in edge_anchor:
                edge_to_draw = [lists[graph_index]]
                nx.draw_networkx_edges(graph, position, edge_to_draw, width=8, edge_color="gray", alpha=0.3, ax=ax)
            
            ## Draw each edge individually, color is based on their bond type
            for edge in graph.edges:
                edge_bond_type = bond_types[graph_index][edge]

                if edge_bond_type == 's':
                    nx.draw_networkx_edges(graph, position, [edge], width=1.0, edge_color="black", ax=ax)

                elif edge_bond_type == 'd':
                    nx.draw_networkx_edges(graph, position, [edge], width=3.0, edge_color="black", ax=ax)
                    nx.draw_networkx_edges(graph, position, [edge], width=1.3, edge_color="white", ax=ax)

                elif edge_bond_type == 't':
                    nx.draw_networkx_edges(graph, position, [edge], width=7.0, edge_color="black", ax=ax)
                    nx.draw_networkx_edges(graph, position, [edge], width=4.25, edge_color="white", ax=ax)
                    nx.draw_networkx_edges(graph, position, [edge], width=1.6, edge_color="black", ax=ax)

                elif edge_bond_type == "q":
                    nx.draw_networkx_edges(graph, position, [edge], width=10, edge_color="black", ax=ax)
                    nx.draw_networkx_edges(graph, position, [edge], width=7, edge_color="white", ax=ax)
                    nx.draw_networkx_edges(graph, position, [edge], width=4, edge_color="black", ax=ax)
                    nx.draw_networkx_edges(graph, position, [edge], width=1.00, edge_color="white", ax=ax)
                
                ## aromatic bond
                else:
                    nx.draw_networkx_edges(graph, position, [edge], width=3, edge_color="black", ax=ax, style="--")
            
            ## Draw mapped edge labels
            label_iterator = 0
            for lists in mapping:
                label = {lists[graph_index]: label_string[label_iterator]}
                nx.draw_networkx_edge_labels(graph, position, label, ax=ax, font_color="lightseagreen", font_family="calibri", font_size=15)
                label_iterator += 1

            
            # ## Draw nodes and atom labels
            nx.draw_networkx_nodes(graph, position, graph.nodes, node_size=175, node_color="white", ax=ax)
            nx.draw_networkx_labels(graph, position, labels=atom_types[graph_index], ax=ax)
            graph_index += 1

        ## Remove x- and y-axis
        plt.axis('off')
        plt.show()




## MCGREGOR DRAWINGS
def draw_mcgregor_mcs_graphs(G, H, mapping, marcs, anchor={}):
    """
    Provides a graphical representation of graphs G and H, highlighting
    their subgraph denoted by the 'mapping', 'marcs' and an optional anchor mapping.

        `Parameters`:
            G (Graph): A NetworkX graph, nodes are integers but may be decorated with items
            H (Graph): A NetworkX graph, nodes are integers but may be decorated with items
            mapping (dict: int -> int): The node correspondence for the MCS
            marcs (np.array): The MARCS array for the MCS

        `Optional`:
            anchor (dict: int -> int): The initial anchor point. Can be declared with dict([(a, b), (b, c), ...])
    """
    
    
    ## String used for nodes mapped to each other
    label_string="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"

    G_pos = nx.spring_layout(G, seed=101)
    H_pos = nx.spring_layout(H, seed=101)

    # #                                                 DRAWING/PLOTTING
    
    ## Options for node drawing
    options = {"edgecolors": "tab:gray", "node_size": 300, "alpha": 1}
    
    # H node colors
    H_node_colors = ["tab:pink" for i in range(len(H.nodes))]
    G_node_colors = ["tab:blue" for i in range(len(G.nodes))]

    ## Color anchor point in G
    for g_nodes in anchor:
        G_node_colors[g_nodes] = 'tab:orange'
    
    ## Determine which edges exist in the mapping and color those.
    G_edges = list(G.edges)
    G_color_edges = []
    H_edges = list(H.edges)
    H_color_edges = []
    for i in range(len(marcs)):
        for j in range(len(marcs[0])):
            if marcs[i][j] == 1:
                G_color_edges.append(G_edges[i])
                H_color_edges.append(H_edges[j])

    ## Add labels to mapped nodes from G to H
    ## such that the mapped nodes have the same label
    iterator = 0
    H_node_labels = {}
    G_node_labels = {}
    for nodes in G.nodes:
        if nodes in mapping:
            G_node_labels[nodes] = str(nodes) + "," + label_string[iterator]
            H_node_labels[mapping[nodes]] = str(mapping[nodes]) + "," + label_string[iterator]
            ## mapped nodes in H will have the same color as the nodes in G
            H_node_colors[mapping[nodes]] = "tab:blue"
            iterator += 1
        ## otherwise, the nodes in G will have the label according to itself (0 - n in principel)
        else:
            G_node_labels[nodes] = nodes
    ## same labeling in H with non-mapped nodes having numbers from 0-n if not already mapped.
    for nodes in H.nodes:
        if nodes not in H_node_labels:
            H_node_labels[nodes] = nodes
    
    ## Color anchor point in H AFTER mapping to overwrite color
    for g_node in anchor:
        H_node_colors[anchor[g_node]] = "tab:orange"

    ## G DRAWING
    subax1 = plt.subplot(121)
    subax1.set_title("G")
    nx.draw_networkx_nodes(G, G_pos, nodelist=sorted(G.nodes), node_color=G_node_colors, **options)
    nx.draw_networkx_edges(G, G_pos, width=1.0, alpha=0.5)
    nx.draw_networkx_edges(
        G,
        G_pos,
        edgelist=G_color_edges,
        width=8,
        alpha=0.5,
        edge_color="tab:blue",
    )
    nx.draw_networkx_labels(G, G_pos, G_node_labels, font_size=10, font_color="whitesmoke")


    ## H DRAWING
    subax2 = plt.subplot(122)
    subax2.set_title("H")
    nx.draw_networkx_nodes(H, nodelist=sorted(H.nodes), node_color=H_node_colors, pos=H_pos)
    nx.draw_networkx_edges(H, H_pos, width=1.0, alpha=0.5)
    nx.draw_networkx_edges(
        H,
        H_pos,
        edgelist=H_color_edges,
        width=8,
        alpha=0.5,
        edge_color="tab:blue",
    )
    nx.draw_networkx_labels(H, H_pos, H_node_labels, font_size=10, font_color="whitesmoke")

    plt.show()


def draw_one_graph(G):

    G_pos = nx.spring_layout(G, seed=101)
    nx.draw_networkx_nodes(G, G_pos, nodelist=G.nodes, node_color="tab:blue")
    nx.draw_networkx_labels(G, G_pos, {i:i for i in G.nodes}, font_size=15, font_color="whitesmoke")
    nx.draw_networkx_edges(G, G_pos, width=1.0, alpha=0.75)
    nx.draw_networkx_edge_labels(G, G_pos, {i:i for i in G.edges}, font_color="lightseagreen")

    plt.show()