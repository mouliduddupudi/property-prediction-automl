'''
Scaffold split is based on the scaffold of the molecules so that train/val/test set is more structurally different.
'''

from random import Random
from collections import defaultdict
from rdkit import Chem
from tqdm import tqdm

from rdkit.Chem import Scaffolds
from rdkit.Chem.Scaffolds import MurckoScaffold
def scaffold_split(df, seed = 2021, frac = [0.8,0.0,0.2], entity:str = "SMILES"):
    """
    reference: https://github.com/chemprop/chemprop/blob/master/chemprop/data/scaffold.py
    Split the dataset on basis of different fragments or scaffolds
    @params:
        frac: split ratios in dataset
        df:   Input dataset DataFrame 
    @return:
        train, validation, test data in dictionary
    """
    random = Random(seed)

    s = df[entity].values
    scaffolds = defaultdict(set)
    idx2mol = dict(zip(list(range(len(s))) ,s))

    error_smiles = 0
    for i, smiles in tqdm(enumerate(s), total=len(s)):
        try:
            scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol = Chem.MolFromSmiles(smiles), includeChirality = False)
            scaffolds[scaffold].add(i)
        except:
            print(smiles + ' returns RDKit error and is thus omitted...')
            error_smiles += 1

    train, val, test = [], [], []
    train_size = int((len(df) - error_smiles) * frac[0])
    val_size = int((len(df) - error_smiles) * frac[1])
    test_size = (len(df) - error_smiles) - train_size - val_size
    train_scaffold_count, val_scaffold_count, test_scaffold_count = 0, 0, 0

    index_sets = list(scaffolds.values())
    big_index_sets = []
    small_index_sets = []
    for index_set in index_sets:
        if len(index_set) > val_size / 2 or len(index_set) > test_size / 2:
            big_index_sets.append(index_set)
        else:
            small_index_sets.append(index_set)
    random.seed(seed)
    random.shuffle(big_index_sets)
    random.shuffle(small_index_sets)
    index_sets = big_index_sets + small_index_sets

    if frac[2] == 0:
        for index_set in index_sets:
            if len(train) + len(index_set) <= train_size:
                train += index_set
                train_scaffold_count += 1
            else:
                val += index_set
                val_scaffold_count += 1
    else:
        for index_set in index_sets:
            if len(train) + len(index_set) <= train_size:
                train += index_set
                train_scaffold_count += 1
            elif len(val) + len(index_set) <= val_size:
                val += index_set
                val_scaffold_count += 1
            else:
                test += index_set
                test_scaffold_count += 1

    
    return {'train': df.iloc[train].reset_index(drop = True),
            'valid': df.iloc[val].reset_index(drop = True),
            'test': df.iloc[test].reset_index(drop = True)}
  
  
# the dataset is df and column with smiles data is SMILES
#splitted_data  = scaffold_split(df, seed = 2021, frac = [0.8,0.0,0.2], entity = "SMILES") 
