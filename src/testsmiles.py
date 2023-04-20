from pysmiles import read_smiles
import networkx as nx
from matplotlib import pyplot as plt
from cliques import mcs_list_leviBarrowBurstall, iterative_approach
from draw_graphs import draw_molecules

def merge_graphs_from_smiles(smiles_list):

    if not smiles_list: return None
    
    graph_list = [read_smiles(smile, explicit_hydrogen=True) for smile in smiles_list]

    n_graphs = len(graph_list)
    to_union_1 = graph_list[0]
    for i in range(1, n_graphs):
        to_union_2 = graph_list[i]
        relabel_offset = len(to_union_1.nodes)
        to_union_2_copy = nx.relabel_nodes(to_union_2, {i: relabel_offset + i for i in to_union_2.nodes})

        to_union_1 = nx.union(to_union_1, to_union_2_copy)
    
    return to_union_1

def convert_labels(graph):
    re_map_edges = {
        1: "s",
        2: "d",
        3: "t",
        4: "q",
        1.5: "a"
    }

    edge_order = nx.get_edge_attributes(graph, "order")
    node_elements = nx.get_node_attributes(graph, "element")
    node_charges = nx.get_node_attributes(graph, "charge")

    atom_labels = { node: { "atom_type": str(node_elements[node]) + str(node_charges[node])} if node_charges[node] != 0 else {"atom_type": str(node_elements[node])} for node in graph.nodes}
    bond_labels = { edge: { "bond_type": re_map_edges[edge_order[edge]]} for edge in graph.edges}

    nx.set_node_attributes(graph, atom_labels)
    nx.set_edge_attributes(graph, bond_labels)





if __name__ == "__main__":

    nad = "C1=CC(=C[N+](=C1)C2C(C(C(O2)COP(=O)([O-])OP(=O)([O-])OCC3C(C(C(O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O)C(=O)N"
    mol_1 = read_smiles(nad, explicit_hydrogen=True)


    nadph = "C1C=CN(C=C1C(=O)N)C2C(C(C(O2)COP(=O)([O-])OP(=O)([O-])OCC3C(C(C(O3)N4C=NC5=C(N=CN=C54)N)OP(=O)([O-])[O-])O)O)O"
    mol_2 = read_smiles(nadph, explicit_hydrogen=True)

    twopg = "C(C(C(=O)[O-])OP(=O)([O-])[O-])O"
    mol_3 = read_smiles(twopg, explicit_hydrogen=True)

    water = "O"
    mol_4 = read_smiles(water, explicit_hydrogen=True)

    fdp = "C(C1C(C(C(O1)(COP(=O)([O-])[O-])O)O)O)OP(=O)([O-])[O-]"
    mol_5 = read_smiles(fdp)

    mol_3.add_edges_from([(5, 14), (6, 10)])
    nx.set_edge_attributes(mol_3,
        {
            (5, 14): {"order": 1},
            (6, 10): {"order": 1},
            (14, 5): {"order": 1},
            (10, 6): {"order": 1}
        }
    )
    nad_nadph_list = [nad, nadph]
    union_graph = merge_graphs_from_smiles(nad_nadph_list)
    union_graph.add_edges_from([(64, 110), (37, 111)])   # Adding the edges from H that are not in G
    nx.set_edge_attributes(union_graph,
        {
            (64,110): {"order": 1},
            (37, 111): {"order": 1},
            (110, 64): {"order": 1},
            (111, 37): {"order": 1}
        }
    )

    water_fdp_list = [water, fdp]
    water_fdp = merge_graphs_from_smiles(water_fdp_list)
    water_fdp.add_edges_from([(1, 10), (11, 0)])
    nx.set_edge_attributes(water_fdp,
                           {
                                (1, 10): {"order": 1},
                                (0, 11): {"order": 1}
                           })
    convert_labels(water_fdp)    

        # node_labels = {node: str(node) + nx.get_node_attributes(water_fdp, "atom_type")[node] for node in water_fdp.nodes}

        # pos = nx.spring_layout(water_fdp, seed=101)
        # nx.draw_networkx_edges(water_fdp, pos, water_fdp.edges)
        # nx.draw_networkx_nodes(water_fdp, pos, water_fdp.nodes)
        # nx.draw_networkx_labels(water_fdp, pos, labels=node_labels)
        # plt.show()


    convert_labels(union_graph)
    convert_labels(mol_3)


    edge_anchor = [
        [(110, 111), (6, 10), (10, 11)],
        [(37, 111), (5, 6), (0, 11)],
        [(37, 64), (5, 14), (0, 1)],
        [(64, 110), (10, 14), (1, 10)]
    ]


    res_iterative = iterative_approach([union_graph, mol_3, water_fdp], edge_anchor=edge_anchor, molecule=True)
    print("Iterative done")
    # res_mass = mcs_list_leviBarrowBurstall([union_graph, mol_3, water_fdp], edge_anchor=edge_anchor, molecule=True)
    # print("Mass, done")

    filtered_res_iterative = []
    for val in res_iterative:
        if val not in filtered_res_iterative:
            filtered_res_iterative.append(val)

    print(len(filtered_res_iterative))

    # filtered_res_mass = []
    # for val in res_mass:
    #     if val not in filtered_res_mass:
    #         filtered_res_mass.append(val)
    
    # print(f"Number of mappings for mass: {len(filtered_res_mass)}, number of mappings for iterative: {len(filtered_res_iterative)}")


    # iterator = 0
    # for val in filtered_res_iterative:
    #     print(f"val: {val}")
    #     draw_molecules([union_graph, mol_3, water_fdp], [val], edge_anchor)
    #     iterator += 1
    #     if iterator == 10: break