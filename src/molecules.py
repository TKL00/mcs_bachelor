import networkx as nx


def propanic_acid():
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
                         7: {"atom_type":"H"},
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

    return propanic_acid

def methane_acid():
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

    return methane_acid



def methanol():
    methanol = nx.Graph()
    methanol.add_edges_from([(0, 1), (0, 3), (0, 4), (0,5), (1, 2)])
    methanol_node_attributes = {
                         0: {"atom_type":"C"},
                         1: {"atom_type":"O"},
                         2: {"atom_type":"H"},
                         3: {"atom_type":"H"},
                         4: {"atom_type":"H"},
                         5: {"atom_type":"H"},
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

    return methanol

def glucose():
    glucose = nx.Graph()
    glucose.add_edges_from([(0,1), (0,5), (0, 7), (0, 12), (1, 2), (1, 13), (1, 14), (2, 3), (2, 16), (2, 18), (3 , 4), (3, 6), 
                            (3, 19), (4, 5), (4, 21), (4, 22), (7, 8), (7, 9), (7, 10), (10, 11), (14, 15), (16, 17), (19, 20), (22, 23)])
    glucose_node_attributes = {
                        0: {"atom_type": "C"},
                        1: {"atom_type": "C"},
                        2: {"atom_type": "C"},
                        3: {"atom_type": "C"},
                        4: {"atom_type": "C"},
                        5: {"atom_type": "O"},
                        6: {"atom_type": "H"},
                        7: {"atom_type": "C"},
                        8: {"atom_type": "H"},
                        9: {"atom_type": "H"},
                        10: {"atom_type": "O"},
                        11: {"atom_type": "H"},
                        12: {"atom_type": "H"},
                        13: {"atom_type": "H"},
                        14: {"atom_type": "O"},
                        15: {"atom_type": "H"},
                        16: {"atom_type": "O"},
                        17: {"atom_type": "H"},
                        18: {"atom_type": "H"},
                        19: {"atom_type": "O"},
                        20: {"atom_type": "H"},
                        21: {"atom_type": "H"},
                        22: {"atom_type": "O"},
                        23: {"atom_type": "H"}        
    }
    glucose_edge_attributes = {
                        (0,1): {"bond_type":"s"}, 
                        (0,5): {"bond_type":"s"}, 
                        (0, 7): {"bond_type":"s"}, 
                        (0, 12): {"bond_type":"s"}, 
                        (1, 2): {"bond_type":"s"}, 
                        (1, 13): {"bond_type":"s"}, 
                        (1, 14): {"bond_type":"s"}, 
                        (2, 3): {"bond_type":"s"}, 
                        (2, 16): {"bond_type":"s"}, 
                        (2, 18): {"bond_type":"s"}, 
                        (3 , 4): {"bond_type":"s"}, 
                        (3, 6): {"bond_type":"s"}, 
                        (3, 19): {"bond_type":"s"}, 
                        (4, 5): {"bond_type":"s"}, 
                        (4, 21): {"bond_type":"s"}, 
                        (4, 22): {"bond_type":"s"}, 
                        (7, 8): {"bond_type":"s"}, 
                        (7, 9): {"bond_type":"s"}, 
                        (7, 10): {"bond_type":"s"}, 
                        (10, 11): {"bond_type":"s"}, 
                        (14, 15): {"bond_type":"s"}, 
                        (16, 17): {"bond_type":"s"}, 
                        (19, 20): {"bond_type":"s"}, 
                        (22, 23): {"bond_type":"s"}
    }
    nx.set_node_attributes(glucose, glucose_node_attributes)
    nx.set_edge_attributes(glucose, glucose_edge_attributes)

    return glucose
    
def caffeine():
    g = nx.Graph()
    g.add_edges_from([
        (0, 1),
        (0, 2),
        (0, 23),
        (2, 3),
        (2, 7),
        (3, 4),
        (3, 5),
        (3, 6),
        (7, 8),
        (7, 9),
        (9, 10),
        (9, 14),
        (10, 11),
        (10, 12),
        (10, 13),
        (14, 15),
        (14, 23),
        (15, 16),
        (16, 17),
        (16, 18),
        (18, 19),
        (18, 23),
        (19, 20),
        (19, 21),
        (19, 22)
    ])
    g_node_attr = {
        0: {"atom_type": "C"},
        1: {"atom_type": "O"},
        2: {"atom_type": "N"},
        3: {"atom_type": "C"},
        4: {"atom_type": "H"},
        5: {"atom_type": "H"},
        6: {"atom_type": "H"},
        7: {"atom_type": "C"},
        8: {"atom_type": "O"},
        9: {"atom_type": "N"},
        10: {"atom_type": "C"},
        11: {"atom_type": "H"},
        12: {"atom_type": "H"},
        13: {"atom_type": "H"},
        14: {"atom_type": "C"},
        15: {"atom_type": "N"},
        16: {"atom_type": "C"},
        17: {"atom_type": "H"},
        18: {"atom_type": "N"},
        19: {"atom_type": "C"},
        20: {"atom_type": "H"},
        21: {"atom_type": "H"} ,
        22: {"atom_type": "H"},
        23: {"atom_type": "C"},
    }
    g_edge_attr = {
        (0, 1): {"bond_type":"d"},
        (0, 2): {"bond_type":"s"},
        (0, 23): {"bond_type":"s"},
        (2, 3): {"bond_type":"s"},
        (2, 7): {"bond_type":"s"},
        (3, 4): {"bond_type":"s"},
        (3, 5): {"bond_type":"s"},
        (3, 6): {"bond_type":"s"},
        (7, 8): {"bond_type":"d"},
        (7, 9): {"bond_type":"s"},
        (9, 10): {"bond_type":"s"},
        (9, 14): {"bond_type":"s"},
        (10, 11): {"bond_type":"s"},
        (10, 12): {"bond_type":"s"},
        (10, 13): {"bond_type":"s"},
        (14, 15): {"bond_type":"s"},
        (14, 23): {"bond_type":"d"},
        (15, 16): {"bond_type":"d"},
        (16, 17): {"bond_type":"s"},
        (16, 18): {"bond_type":"s"},
        (18, 19): {"bond_type":"s"},
        (18, 23): {"bond_type":"s"},
        (19, 20): {"bond_type":"s"},
        (19, 21): {"bond_type":"s"},
        (19, 22): {"bond_type":"s"}
    }

    nx.set_node_attributes(g, g_node_attr)
    nx.set_edge_attributes(g, g_edge_attr)

    return g

def methane():

    methane = nx.Graph()
    methane.add_edges_from([(0, 1), (0, 2), (0, 3), (0, 4)])
    methane_node_attributes = {
        0: {"atom_type": "C"},
        1: {"atom_type": "H"},
        2: {"atom_type": "H"},
        3: {"atom_type": "H"},
        4: {"atom_type": "H"}
    }
    methane_edge_attributes = {
        (0, 1): {"bond_type": "s"},
        (0, 2): {"bond_type": "s"},
        (0, 3): {"bond_type": "s"},
        (0, 4): {"bond_type": "s"}
    }

    nx.set_node_attributes(methane, methane_node_attributes)
    nx.set_edge_attributes(methane, methane_edge_attributes)

    return methane