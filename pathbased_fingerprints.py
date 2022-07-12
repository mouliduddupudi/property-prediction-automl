from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
from tqdm import tqdm
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
import numpy as np
from rdkit.Chem import MACCSkeys

def fingerprint_mols_process_image(x):
    mol=Chem.MolFromSmiles(x)
    if mol:
        fp = rdmolops.RDKFingerprint(mol, fpSize=2048, minPath=1, maxPath=7)
        arr=fps_to_arr([fp])
        arr=list(arr[0])
        arr.insert(0,x)
            
        return arr
    else:
        return None

def fps_to_arr(fps):
        fp_dim = len(fps[0])
        arr = np.zeros((len(fps), fp_dim), dtype=np.bool)
        for i, fp in enumerate(fps):
            onbits = list(fp.GetOnBits())
            arr[i, onbits] = 1
        return arr
      
      
fps=[fingerprint_mols_process_image(x) for x in data['SMILES']]

fp_columns=[str(i) for i in range(int(len(fps[0])-1))]
fp_columns.insert(0,'SMILES')

dataset = pd.DataFrame(fps,columns=fp_columns)

dataset.replace({False: 0, True: 1}, inplace=True)
dataset.reset_index(drop=True, inplace=True)
