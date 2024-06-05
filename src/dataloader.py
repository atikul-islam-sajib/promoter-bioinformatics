import pandas as pd
from feature_generator import FeatureGenerator

class Loader():
    def __init__(self, dataset = None, test_size = 0.20, shuffle = True, random_state = 42, preprocess=False):
        self.dataset = dataset
        self.test_size = test_size
        self.random_state = random_state
        self.shuffle = shuffle
        self.preprocess = preprocess
        
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
            
            return dataset
        
        else:
            ValueError("Dataset should be in the format of pandas dataFrame".capitalize())
            
    def create_features(self, dataset=None):
        
        di_nucleosides_features = FeatureGenerator().generate_features(
            dataset = dataset,
            type="di"
        )
        
        tri_nucleosides_features = FeatureGenerator().generate_features(
            dataset = dataset,
            type="tri"
        )
        
        tetra_nucleosides_features = FeatureGenerator().generate_features(
            dataset = dataset,
            type="tetra"
        )
        
        print(di_nucleosides_features.shape)
        print(di_nucleosides_features.isnull().sum().sum())
        
        print(di_nucleosides_features.head())
        
        print(dataset.shape)
        print(dataset.isnull().sum().sum())
        print(dataset.head())
        
        print(tri_nucleosides_features.head())
        print(tri_nucleosides_features.shape)
        
        print(tetra_nucleosides_features.head())
        print(tetra_nucleosides_features.shape)
        
    
    
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