import networkx as nx
import matplotlib.pyplot as plt
from linegraph import convert_edge_anchor_lg

### Authors: Tobias Klink Lehn (toleh20@student.sdu.dk) and Kasper Halkjær Beider (kbeid20@student.sdu.dk)
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

    PGH_nodes = list(PGH.nodes)

    for i in range(len(G.nodes)):
        G_node_labels[i] = i

    for i in range(len(H.nodes)):
        H_node_labels[i] = i
    
    for i in range(len(PGH.nodes)):
        PGH_node_labels[PGH_nodes[i]] = PGH_nodes[i]

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

def draw_two_graphs(G, H):
    """
        Draws two NetworkX graphs next to each other.
    """
    G_pos = nx.spring_layout(G, seed=101)
    H_pos = nx.spring_layout(H, seed=101)
    G_node_labels = {}
    H_node_labels = {}

    for i in range(len(G.nodes)):
        G_node_labels[i] = i

    for i in range(len(H.nodes)):
        H_node_labels[i] = i


    subax1 = plt.subplot(121)
    subax1.set_title("G")
    nx.draw_networkx_nodes(G, G_pos, nodelist=G.nodes)
    nx.draw_networkx_labels(G, G_pos, G_node_labels, font_size=15, font_color="whitesmoke")
    nx.draw_networkx_edges(G, G_pos, edgelist=G.edges)
    
    subax2 = plt.subplot(122)
    subax2.set_title("H")
    nx.draw_networkx_nodes(H, H_pos, nodelist=H.nodes)
    nx.draw_networkx_labels(H, H_pos, H_node_labels, font_size=15, font_color="whitesmoke")
    nx.draw_networkx_edges(H, H_pos, edgelist=H.edges)
    
    plt.show()
    

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
    G_labels = {}
    H_labels = {}

    ## Assigning different labels for the nodes in the mapping
    ## a mapping u -> v from G to H assigns the same string to u and v 
    for g_nodes in mapping:
        G_labels[g_nodes] = label_string[g_nodes]
        H_labels[mapping[g_nodes]] = label_string[g_nodes]

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

    ## G DRAWING
    subax1 = plt.subplot(121)
    nx.draw_networkx_nodes(G, G_pos, nodelist=sorted(G.nodes), node_color=G_node_colors, **options)
    nx.draw_networkx_edges(G, G_pos, width=1.0, alpha=0.5)
    
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

    nx.draw_networkx_edges(
        G,
        G_pos,
        edgelist=G_color_edges,
        width=8,
        alpha=0.5,
        edge_color="tab:blue",
    )

    ## Add labels to mapped nodes from G to H
    ## such that the mapped nodes have the same label
    iterator = 0
    H_node_labels = {}
    G_node_labels = {}
    for nodes in G.nodes:
        if nodes in mapping:
            G_node_labels[nodes] = label_string[iterator]
            H_node_labels[mapping[nodes]] = label_string[iterator]
            ## mapped nodes in H will have the same color as the nodes in G
            H_node_colors[mapping[nodes]] = "tab:blue"
            iterator += 1
        ## otherwise, G nodes will be numbered from 0-n
        else:
            G_node_labels[nodes] = nodes
    ## same labeling in H with non-mapped nodes having numbers from 0-n if not already mapped.
    for nodes in H.nodes:
        if nodes not in H_node_labels:
            H_node_labels[nodes] = nodes
    
    ## Color anchor point in H AFTER mapping to overwrite color
    for g_node in anchor:
        H_node_colors[anchor[g_node]] = "tab:orange"

    nx.draw_networkx_labels(G, G_pos, G_node_labels, font_size=15, font_color="whitesmoke")


    ## H DRAWING
    subax2 = plt.subplot(122)
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
    nx.draw_networkx_labels(H, H_pos, H_node_labels, font_size=15, font_color="whitesmoke")

    plt.show()


def draw_mcgregor_mcs_lgs(G, H, LG, LH, mapping, marcs, edge_anchor={}):

    """
        Provides a graphical representation of graphs G, H and their corresponding linegraphs LG, LH, highlighting
        their subgraph denoted by the 'mapping', 'marcs' and an optional edge anchor mapping between G and H.
        The anchored edges between G and H are highlighted, and their corresponding nodes in the line graphs are
        highlighted accordingly.

        `Parameters`:
            G (Graph): A NetworkX graph, nodes are integers but may be decorated with items
            H (Graph): A NetworkX graph, nodes are integers but may be decorated with items
            LG (Graph): The linegraph of G, nodes are integers but may be decorated with items
            LH (Graph): The linegraph of H, nodes are integers but may be decorated with items
            mapping (dict: int -> int): The node correspondence for the MCS in LG and LH.
            marcs (np.array): The MARCS array for the MCS

        `Optional`:
            edge_anchor (dict: (int, int) -> (int, int)): The initial anchored edges between G and H.
    """

    label_string="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
    line_node_anchor = convert_edge_anchor_lg(G, H, edge_anchor)
    G_edge_labels, H_edge_labels, G_node_labels, H_node_labels, LG_node_labels, LH_node_labels = {}, {}, {}, {}, {}, {}

    G_edges = list(G.edges)
    H_edges = list(H.edges)
    LG_edges = list(LG.edges)
    LH_edges = list(LH.edges)

    G_pos = nx.spring_layout(G, seed=1010203)
    H_pos = nx.spring_layout(H, seed=1010203)
    LG_pos = nx.spring_layout(LG, seed=1010203)
    LH_pos = nx.spring_layout(LH, seed=1010203)

    ## Node colours
    H_node_colors = ["tab:pink" for i in range(len(H.nodes))]
    G_node_colors = ["tab:blue" for i in range(len(G.nodes))]
    LH_node_colors = ["tab:pink" for i in range(len(LH.nodes))]
    LG_node_colors = ["tab:blue" for i in range(len(LG.nodes))]

    # Edge colors in line graph mapping
    LG_color_edges = []
    LH_color_edges = []
    for i in range(len(marcs)):
        for j in range(len(marcs[0])):
            if marcs[i][j] == 1:
                LG_color_edges.append(LG_edges[i])
                LH_color_edges.append(LH_edges[j])

    ## edge color for anchored edges
    G_anchor_color_edges = [keys for keys in edge_anchor]
    H_anchor_color_edges = [edge_anchor[keys] for keys in edge_anchor]

    ## node color for anchored nodes in line graphs
    LG_anchor_color_nodes = [keys for keys in line_node_anchor]
    LH_anchor_color_nodes = [line_node_anchor[keys] for keys in line_node_anchor]

    ## 1-n numbering of the nodes
    for i in range(len(G.nodes)):
        G_node_labels[i] = i
    
    for i in range(len(H.nodes)):
        H_node_labels[i] = i

    for i in range(len(LG.nodes)):
        LG_node_labels[i] = i
    
    for i in range(len(LH.nodes)):
        LH_node_labels[i] = i

    ## Compute labels for the edges and nodes in the mapping
    for LG_node in mapping:
        LH_node = mapping[LG_node]

        LG_node_labels[LG_node] = label_string[LG_node]
        LH_node_labels[LH_node] = label_string[LG_node]
        
        G_edge_labels[G_edges[LG_node]] = label_string[LG_node]
        H_edge_labels[H_edges[LH_node]] = label_string[LG_node]
    
    ###                 DRAWING

    ## G
    subax1 = plt.subplot(221)
    subax1.set_title("G")
    nx.draw_networkx_nodes(G, G_pos, nodelist=sorted(G.nodes), node_color=G_node_colors)
    nx.draw_networkx_edges(G, G_pos, width=1.0, alpha=0.5)
    ## color anchored edges
    nx.draw_networkx_edges(G, G_pos, width=8, alpha=0.5, edgelist=G_anchor_color_edges, edge_color="tab:orange")
    nx.draw_networkx_labels(G, G_pos, G_node_labels, font_size=15, font_color="whitesmoke")
    nx.draw_networkx_edge_labels(G, G_pos, G_edge_labels)

    ## H
    subax2 = plt.subplot(222)
    subax2.set_title("H")
    nx.draw_networkx_nodes(H, H_pos, nodelist=sorted(H.nodes), node_color=H_node_colors)
    nx.draw_networkx_edges(H, H_pos, width=1.0, alpha=0.5)
    ## color anchored edges
    nx.draw_networkx_edges(H, H_pos, width=8, alpha=0.5, edgelist=H_anchor_color_edges, edge_color="tab:orange")
    nx.draw_networkx_labels(H, H_pos, H_node_labels, font_size=15, font_color="whitesmoke")
    nx.draw_networkx_edge_labels(H, H_pos, H_edge_labels)

    ## LG
    subax3 = plt.subplot(223)
    subax3.set_title("Linegraph of G")
    nx.draw_networkx_nodes(LG, LG_pos, nodelist=sorted(LG.nodes), node_color=LG_node_colors)
    nx.draw_networkx_nodes(LG, LG_pos, nodelist=LG_anchor_color_nodes, node_color="tab:orange")
    nx.draw_networkx_edges(LG, LG_pos, width=1.0, alpha=0.5)
    nx.draw_networkx_edges(LG, LG_pos, edgelist=LG_color_edges, width=5, alpha=0.5, edge_color="tab:red")
    nx.draw_networkx_labels(LG, LG_pos, LG_node_labels, font_size=15, font_color="whitesmoke")

    ## LH
    subax4 = plt.subplot(224)
    subax4.set_title("Linegraph of H")
    nx.draw_networkx_nodes(LH, LH_pos, nodelist=sorted(LH.nodes), node_color=LH_node_colors)
    nx.draw_networkx_nodes(LH, LH_pos, nodelist=LH_anchor_color_nodes, node_color="tab:orange")
    nx.draw_networkx_edges(LH, LH_pos, width=1.0, alpha=0.5)
    nx.draw_networkx_edges(LH, LH_pos, edgelist=LH_color_edges, width=5, alpha=0.5, edge_color="tab:red")
    nx.draw_networkx_labels(LH, LH_pos, LH_node_labels, font_size=15, font_color="whitesmoke")
    

    plt.show()