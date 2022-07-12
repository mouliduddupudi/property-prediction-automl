'''
The following code represents the smiles data in descriptors form.
'''
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
from tqdm import tqdm
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors

def descriptors(data):
    print("desp_calculation")
    descriptors_list = [x[0] for x in Descriptors._descList]
    desp_df=pd.DataFrame(columns=descriptors_list)
    calc = MoleculeDescriptors.MolecularDescriptorCalculator([x[0] for x in Descriptors._descList])
    for i,r in tqdm(enumerate(data)):
        try:
            ds = list(calc.CalcDescriptors(Chem.MolFromSmiles(r)))
            desp_df.loc[i]=ds
        except:
            continue
    
    return desp_df
