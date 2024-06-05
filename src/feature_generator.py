import pandas as pd
from tqdm import tqdm
from utils import di_nucleosides, tri_nucleosides, tetra_nucleosides


import warnings
warnings.filterwarnings('ignore')


class FeatureGenerator():
    def __init__(self):
        
        self.task = "Feature Generator for Promoters".title()
    
        
    def generate_features(self, dataset = None, type = None):
        if isinstance(dataset, pd.pandas.DataFrame):
            if type == "di":
                self.nucleosides = di_nucleosides
            
                for index in tqdm(range(dataset.shape[0])):
                    self.value = []
                    self.sequence = dataset.loc[index, "sequence"]
                    
                    for idx in range(len(self.sequence) - 1 - 1): # (1: GAP, 1: Di -->Di)
                        for _, nucleosides in enumerate(self.nucleosides):
                            seq = self.sequence[idx] + self.sequence[idx + 2] # A_A
                                
                            if seq == nucleosides:
                                self.value.append(1)
                                
                            else:
                                self.value.append(0)               
                                
                    for idx, val in enumerate(self.value):
                        dataset.loc[index, f"1_GAP_KMer_di_{idx}"]=val
                        
                for index in tqdm(range(dataset.shape[0])):
                    self.value = []
                    self.sequence = dataset.loc[index, "sequence"]
                    
                    for idx in range(len(self.sequence) - 2 - 1): # (2: GAP, 1: Di -->Di)
                        for _, nucleosides in enumerate(self.nucleosides):
                            seq = self.sequence[idx] + self.sequence[idx + 3] # A__A
                                
                            if seq == nucleosides:
                                self.value.append(1)
                                
                            else:
                                self.value.append(0)
                                
                    for idx, val in enumerate(self.value):
                        dataset.loc[index, f"2_GAP_KMer_di_{idx}"]=val
                        
                for index in tqdm(range(dataset.shape[0])):
                    self.value = []
                    self.sequence = dataset.loc[index, "sequence"]
                    
                    for idx in range(len(self.sequence) - 3 - 1): # (3: GAP, 1: Di -->Di)
                        for _, nucleosides in enumerate(self.nucleosides):
                            seq = self.sequence[idx] + self.sequence[idx + 4] # A___A
                                
                            if seq == nucleosides:
                                self.value.append(1)
                                
                            else:
                                self.value.append(0)
                                
                    for idx, val in enumerate(self.value):
                        dataset.loc[index, f"3_GAP_KMer_di_{idx}"]=val
                
                print("di nucleosides features generation is completed.".title())        
                        
                return dataset
                        
                
            elif type == "tri":
                self.nucleosides = tri_nucleosides
                pass
            
            elif type == "tetra":
                self.helper = 3
                self.Gap = 3
                self.nucleosides = tetra_nucleosides