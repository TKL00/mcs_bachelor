import copy

class Workspace:
    """
    Workspace class for the McGregor algorithm. 


    `MARCS`: A |E1| x |E2| matrix with a 1 in position (r,s) if edge `r` can be mapped to edge `s`

    `arcsleft`: An integer value that denotes the number of rows in MARCS that are not all 0s.

    `MARCS_ones_left`: An array of length (rows_in_MARCS) where MARCS_ones_left[i] denotes the number of 1s left in row `i`

    `edges_killed`: An array of integer tuples (r, s) such that all tuples correspond to 'killed' edges that lead
    to the current state of MARCS.
    """

    def __init__(self, MARCS, arcsleft, MARCS_ones_left, edges_killed):
        self.MARCS = copy.deepcopy(MARCS)
        self.arcsleft = arcsleft
        self.MARCS_ones_left = copy.deepcopy(MARCS_ones_left)
        self.edges_killed = edges_killed
    
    def get_MARCS(self):
        return self.MARCS
    
    def get_arcsleft(self):
        return self.arcsleft

    def get_MARCS_ones_left(self):
        return self.MARCS_ones_left

    def get_edges_killed(self):
        return self.edges_killed