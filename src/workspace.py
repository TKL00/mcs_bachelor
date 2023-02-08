import copy

class Workspace:

    def __init__(self, MARCS, arcsleft, MARCS_ones_left):
        self.MARCS = copy.deepcopy(MARCS)
        self.arcsleft = arcsleft
        self.MARCS_ones_left = copy.deepcopy(MARCS_ones_left)
    
    def get_MARCS(self):
        return self.MARCS
    
    def get_arcsleft(self):
        return self.arcsleft