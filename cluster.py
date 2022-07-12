'''
Creating clusters for smiles data using rdkit module. The initial dataset will be modified with a label columns specified to which cluster the molecule\\
belongs. Based on this, different clusters can be divided and trained for one's specific purpose.
'''


import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from statistics import mean
from rdkit import DataStructs

from rdkit.ML.Cluster import Butina
#Define clustering setup
def ClusterFps(fps,cutoff=0.2):
    # first generate the distance matrix:
    dists = []
    nfps = len(fps)
    for i in range(1,nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i])
        #print(sims)
        dists.extend([1-x for x in sims])

    # now cluster the data:
    print(mean(dists))
    cs = Butina.ClusterData(dists,nfps,cutoff,isDistData=True,reordering=True)#cutoff==distance cutoff
    return cs

rationale=list(df['SMILES'])     #--> takes smiles data from input dataset
ms = [Chem.MolFromSmiles(x) for x in rationale]
fps = [AllChem.GetMorganFingerprintAsBitVect(x,2,1024) for x in ms]
clusters=ClusterFps(fps,cutoff=0.8)
print("clusters len: ", len(clusters))

index=[]
mols=[]
df['labels']=''           #--> Adds another column in the original dataset, df, with cluster lables
j=1
for i in clusters:

    for inx,r in df.iterrows():
        if inx in i:
            df.loc[inx,'labels']=j
    j=j+1
    
num_of_clusters = []
for i in clusters:
    num_of_clusters.append(len(i))
print(num_of_clusters)   #--> prints number of molecules in each cluster in descending order
