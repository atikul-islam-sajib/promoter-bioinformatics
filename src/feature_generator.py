import os
import pandas as pd
from utils import config
from tqdm import tqdm
from utils import di_nucleosides, tri_nucleosides, tetra_nucleosides, single_nucleosides


import warnings
warnings.filterwarnings('ignore')


class FeatureGenerator():
    def __init__(self):
        
        self.task = "Feature Generator for Promoters".title()
    
        
    def generate_features(self, dataset = None, type = None):
        if isinstance(dataset, pd.pandas.DataFrame):
            if type == "single":
                for nucleosides in tqdm(single_nucleosides):
                    for i in range(len(dataset.loc[0, "sequence"])):
                        feature_name = 'Index' + str(i) + '_' + nucleosides
                        
                        values = []
                        
                        for j in range(dataset.shape[0]):
                            sequence = dataset.loc[j, "sequence"]
                            if sequence[i] == nucleosides:
                                values.append(1)
                            else:
                                values.append(0)
                                
                        dataset[feature_name] = values
                        
                print("single nucleosides features generation is completed...".title())
                        
                dataset.to_csv(os.path.join(config()["path"]["PROCESSED_DATA_PATH"], "single.csv"))
                        
                return dataset
            
            elif type == "di":
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
                
                print("di nucleosides features generation is completed...".title())
                
                dataset.to_csv(os.path.join(config()["path"]["PROCESSED_DATA_PATH"], "single_di.csv"))        
                        
                return dataset
                        
                
            elif type == "tri":
                self.nucleosides = tri_nucleosides
                
                for index in tqdm(range(dataset.shape[0])):
                    self.value = []
                    self.sequence = dataset.loc[index, "sequence"]
                    
                    for idx in range(len(self.sequence) - 1 - 2): # (1: GAP, 1: tri -->tri)
                        for _, nucleosides in enumerate(self.nucleosides):
                            seq = self.sequence[idx] + self.sequence[idx + 2] + self.sequence[idx + 3] # A_AA
                                
                            if seq == nucleosides:
                                self.value.append(1)
                                
                            else:
                                self.value.append(0)               
                                
                    for idx, val in enumerate(self.value):
                        dataset.loc[index, f"1_GAP_KMer_tri_{idx}"]=val
                        
                        
                for index in tqdm(range(dataset.shape[0])):
                    self.value = []
                    self.sequence = dataset.loc[index, "sequence"]
                    
                    for idx in range(len(self.sequence) - 2 - 2): # (2: GAP, 1: tri -->tri)
                        for _, nucleosides in enumerate(self.nucleosides):
                            seq = self.sequence[idx] + self.sequence[idx + 3] + self.sequence[idx + 4] # A__AA
                                
                            if seq == nucleosides:
                                self.value.append(1)
                                
                            else:
                                self.value.append(0)               
                                
                    for idx, val in enumerate(self.value):
                        dataset.loc[index, f"2_GAP_KMer_tri_{idx}"]=val
                        
                for index in tqdm(range(dataset.shape[0])):
                    self.value = []
                    self.sequence = dataset.loc[index, "sequence"]
                    
                    for idx in range(len(self.sequence) - 3 - 2): # (3: GAP, 1: tri -->tri)
                        for _, nucleosides in enumerate(self.nucleosides):
                            seq = self.sequence[idx] + self.sequence[idx + 4] + self.sequence[idx + 5] # A___AA
                                
                            if seq == nucleosides:
                                self.value.append(1)
                                
                            else:
                                self.value.append(0)               
                                
                    for idx, val in enumerate(self.value):
                        dataset.loc[index, f"3_GAP_KMer_tri_{idx}"]=val
                        
                print("tri nucleosides features generation is completed...".title())
                
                dataset.to_csv(os.path.join(config()["path"]["PROCESSED_DATA_PATH"], "single_di_tri.csv"))
                
                return dataset
            
            elif type == "tetra":
                
                self.nucleosides = tetra_nucleosides
                
                for index in tqdm(range(dataset.shape[0])):
                    self.value = []
                    self.sequence = dataset.loc[index, "sequence"]
                    
                    for idx in range(len(self.sequence) - 1 - 3): # (1: GAP, 3: tetra -->tetra)
                        for _, nucleosides in enumerate(self.nucleosides):
                            seq = self.sequence[idx] + self.sequence[idx + 2] + self.sequence[idx + 3] + self.sequence[idx + 4]# A_AAA
                                
                            if seq == nucleosides:
                                self.value.append(1)
                                
                            else:
                                self.value.append(0)               
                                
                    for idx, val in enumerate(self.value):
                        dataset.loc[index, f"1_GAP_KMer_tetra_{idx}"]=val
                        
                for index in tqdm(range(dataset.shape[0])):
                    self.value = []
                    self.sequence = dataset.loc[index, "sequence"]
                    
                    for idx in range(len(self.sequence) - 2 - 3): # (2: GAP, 3: tetra -->tetra)
                        for _, nucleosides in enumerate(self.nucleosides):
                            seq = self.sequence[idx] + self.sequence[idx + 3] + self.sequence[idx + 4] + self.sequence[idx + 5]# A__AAA
                                
                            if seq == nucleosides:
                                self.value.append(1)
                                
                            else:
                                self.value.append(0)               
                                
                    for idx, val in enumerate(self.value):
                        dataset.loc[index, f"2_GAP_KMer_tetra_{idx}"]=val
                        
                for index in tqdm(range(dataset.shape[0])):
                    self.value = []
                    self.sequence = dataset.loc[index, "sequence"]
                    
                    for idx in range(len(self.sequence) - 3 - 3): # (3: GAP, 3: tetra -->tetra)
                        for _, nucleosides in enumerate(self.nucleosides):
                            seq = self.sequence[idx] + self.sequence[idx + 4] + self.sequence[idx + 5] + self.sequence[idx + 6]# A___AAA
                                
                            if seq == nucleosides:
                                self.value.append(1)
                                
                            else:
                                self.value.append(0)               
                                
                    for idx, val in enumerate(self.value):
                        dataset.loc[index, f"3_GAP_KMer_tetra_{idx}"]=val
                        
                print("tetra nucleosides features generation is completed...".title())
                
                dataset.to_csv(os.path.join(config()["path"]["PROCESSED_DATA_PATH"], "single_di_tri_tetra.csv"))
                        
                return dataset
            
            elif type == "GC":
                self.nucleosides = single_nucleosides
                
                self.GC_Content = []
                self.value = []
                
                for _, sequence in tqdm(enumerate(dataset.loc[:, "sequence"])):
                    G = sequence.count("G".upper())
                    C = sequence.count("C".upper())
                    A = sequence.count("A".upper())
                    T = sequence.count("T".upper())
                    
                    self.GC_Content.append((G + C)/(A + C + G + T))
                    
                    self.GC_Count = (G + C)
                    
                    if self.GC_Count > 10:
                        self.value.append(1)
                        
                    else:
                        self.value.append(0)
                    
                dataset["GC_Content"] = self.GC_Content
                dataset["GC_Content>10"] = self.value
                
                print("GC Content features generation is completed...".title())
                
                dataset.to_csv(os.path.join(config()["path"]["PROCESSED_DATA_PATH"], "final_dataset.csv"))
                        
                return dataset                
            
            else:
                raise TypeError("Please select the type to create the Features for Promoters.".capitalize())
            
            
            
if __name__ == "__main__":
    feature_generator = FeatureGenerator()
    
    for type in ["single", "di", "tri", "tetra", "GC"]:
        _ = feature_generator.generate_features(
            type=type, dataset="./data/processed/dataset.csv"
        )