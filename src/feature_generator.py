import pandas as pd
from utils import di_nucleosides, tri_nucleosides, tetra_nucleosides


class FeatureGenerator():
    def __init__(self):
        
        self.task = "Feature Generator for Promoters".title()
    
        
    def create_features(self, dataset = None, type = None):
        if isinstance(dataset, pd.pandas.DataFrame):
            if type == "di":
                self.nucleosides = di_nucleosides
            
                for index, in range(dataset.shape[0]):
                    self.value = []
                    self.sequence = dataset[index, "sequence"]
                    
                    for idx in (len(self.sequence) - 1 - 1): # (1: GAP, 1: Di -->Di)
                        for _, nucleosides in enumerate(self.nucleosides):
                            seq = self.sequence[idx] + self.sequence[idx + 2] # A_A
                                
                            if seq == nucleosides:
                                self.value.append(1)
                                
                            else:
                                self.value.append(0)
                                
                    for idx, val in enumerate(self.value):
                        dataset[index, f"1_GAP_KMer_di_{idx}"]=val
                        
                for index, in range(dataset.shape[0]):
                    self.value = []
                    self.sequence = dataset[index, "sequence"]
                    
                    for idx in (len(self.sequence) - 2 - 1): # (2: GAP, 1: Di -->Di)
                        for _, nucleosides in enumerate(self.nucleosides):
                            seq = self.sequence[idx] + self.sequence[idx + 3] +  self.sequence[idx + 4] # A__AA
                                
                            if seq == nucleosides:
                                self.value.append(1)
                                
                            else:
                                self.value.append(0)
                                
                    for idx, val in enumerate(self.value):
                        dataset[index, f"2_GAP_KMer_di_{idx}"]=val
                        
                for index, in range(dataset.shape[0]):
                    self.value = []
                    self.sequence = dataset[index, "sequence"]
                    
                    for idx in (len(self.sequence) - 3 - 1): # (3: GAP, 1: Di -->Di)
                        for _, nucleosides in enumerate(self.nucleosides):
                            seq = self.sequence[idx] + self.sequence[idx + 4] + self.sequence[idx + 5] + + self.sequence[idx + 6] # A___AAA
                                
                            if seq == nucleosides:
                                self.value.append(1)
                                
                            else:
                                self.value.append(0)
                                
                    for idx, val in enumerate(self.value):
                        dataset[index, f"3_GAP_KMer_di_{idx}"]=val
                        
                
            elif type == "tri":
                self.nucleosides = tri_nucleosides
                pass
            
            elif type == "tetra":
                self.helper = 3
                self.Gap = 3
                self.nucleosides = tetra_nucleosides