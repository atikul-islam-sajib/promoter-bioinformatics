import os
import pandas as pd
from utils import dump, load, config
from feature_generator import FeatureGenerator

import warnings
warnings.filterwarnings('ignore')

class Loader():
    def __init__(self, dataset = None, test_size = 0.20, shuffle = True, random_state = 42, preprocess=False):
        self.dataset = dataset
        self.test_size = test_size
        self.random_state = random_state
        self.shuffle = shuffle
        self.preprocess = preprocess
        
        
        self.config = config()
        
    def load_dataset(self):
        dataset = pd.read_csv(self.dataset)
        
        if isinstance(dataset, pd.pandas.DataFrame):
            
            if not self.preprocess:
            
                dataset = dataset.rename(columns={dataset.columns[-1]: "sequence" , dataset.columns[0]: "labels"})
                
                if "S10" in dataset.columns:
                    dataset.drop(["S10"], axis = 1, inplace=True)
                else:
                    Exception("S10 Feature is not available in the dataset".capitalize())
                
            dataset["labels"] = dataset["labels"].map({"+": 1, "-": 0}).astype("int")
            dataset["sequence"] = dataset["sequence"].str.replace("\t", "").str.upper()
            
            dataset.to_csv(os.path.join(self.config["path"]["PROCESSED_DATA_PATH"], "dataset.csv"))
            
            dump(
                value=dataset, filename=os.path.join(self.config["path"]["PROCESSED_DATA_PATH"], "dataset.pkl")
            )
            
            return dataset
        
        else:
            ValueError("Dataset should be in the format of pandas dataFrame".capitalize())
            
    def create_features(self, dataset=None):
        
        for type in ["single", "di", "tri", "tetra"]:
            _ = FeatureGenerator().generate_features(dataset = dataset, type=type)
    
    
if __name__ == "__main__":
    loader = Loader(
        dataset="./data/raw/promoters.data",
        test_size=0.25,
        shuffle=True,
        random_state=42,
        preprocess=False
    )
    
    dataset = loader.load_dataset()
    
    loader.create_features(dataset=dataset)