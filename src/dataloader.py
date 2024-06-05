import os
import argparse
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
            
    def split_dataset(self, **kwargs):
        pass
            
    def create_features(self, dataset=None, features_type=None):
        
        if isinstance(dataset, pd.pandas.DataFrame) and isinstance(features_type, list):
        
            for type in features_type:
                _ = FeatureGenerator().generate_features(dataset = dataset, type=type)
                
            print(
                "Features Generation is completed and saved in the path # {}".format(self.config["path"]["PROCESSED_DATA_PATH"])
            )
                
        else:
            raise ValueError(
                "Features cannot be generated: Features type should be in list format and dataFrame should be in pandas format".capitalize()
            )
            
    @staticmethod
    def dataset_details():
        pass
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Dataloader for Promoters".title())
    
    parser.add_argument(
        "--dataset", type=str, default=config()["dataloader"]["dataset"], help="Define the dataset, CSV file".capitalize()
    )
    parser.add_argument(
        "--test_size", type=float, default=config()["dataloader"]["test_size"], help="Define the train and test split of the dataset".capitalize()
    )
    parser.add_argument(
        "--shuffle", type=bool, default=config()["dataloader"]["shuffle"], help="Define whether dataset needs to be shuffle or not".capitalize()
    )
    parser.add_argument(
        "--random_state", type=int, default=config()["dataloader"]["random_state"], help="Define the random state of the dataset".capitalize()
    )
    parser.add_argument(
        "--preprocess", type=bool, default=config()["dataloader"]["preprocess"], help="Define the dataset needs to preprocess or not".capitalize()
    )
    parser.add_argument(
        "--type", type=list, default=config()["dataloader"]["features_type"], help="Define the types of features - how will it be created".capitalize()
    )
    
    args = parser.parse_args()
    
    loader = Loader(
        dataset=args.dataset,
        test_size=args.test_size,
        shuffle=args.shuffle,
        random_state=args.random_state,
        preprocess=args.preprocess
    )
    
    dataset = loader.load_dataset()
    
    loader.create_features(dataset=dataset, features_type=args.type)