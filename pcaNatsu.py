#!/Users/siida/.pyenv/shims/python
import sys
from tqdm import tqdm 
import numpy as np
from scipy.spatial import distance
from MDAnalysis import Universe 

def argParse():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("-i" , "--itraj"     , required=True)
    p.add_argument("-s1", "--selection1", required=True)
    p.add_argument("-s2", "--selection2", required=True)
    p.add_argument("-ir", "--ref")
    args = p.parse_args()
    return args.itraj, args.selection1, args.selection2, args.ref

def distanceMatrix(trajFile, sele1, sele2, ref = None):

    if ref == None: 
        u = Universe(trajFile)

    else:
        print("* Reference is given.")
        u = Universe(ref, trajFile)

    s1, s2 = u.select_atoms(sele1), u.select_atoms(sele2)
    print("* pair 1: \n    ", s1)
    print("* pair 2: \n    ", s2)

    distances = []
    for itraj in tqdm(u.trajectory):
        pos1, pos2 = s1.positions, s2.positions

        if sele1 == sele2: #if self-distance pair calculation
            dist = np.triu(distance.cdist(pos1, pos2, metric='euclidean')) #symmetrical matrix if self-distances 
            dist = dist[dist !=0] # remove 0 elements 

        else: #if differnt distance pair calculation 
            dist = distance.cdist(pos1, pos2, metric='euclidean').flatten()

        distances.append(dist)
    return distances

def run(fileNameTraj, s1, s2, reference = None):
   """
   Args:
       fileNameTraj : an input trajectory whose format is PDB.
       reference    : a reference pdb file 
       s1, s2       : selection 1 and selection 2 ; the selection rule is based on MDAnalysis 
   Returns:
       x: coordinates on PC space (numpy array), pca: PCA object of sklearn.decompositio.PCA that generates x
   """
   from sklearn.decomposition import PCA

   if reference == None:
       distances = distanceMatrix(fileNameTraj, s1, s2)

   else:
       distances = distanceMatrix(fileNameTraj, s1, s2, reference)
   
   pca = PCA()
   #pca = PCA(n_components=5)
   x   = pca.fit_transform(distances)
   return x, pca

def write(x, pca):
   """
   Args:
       x: a numpy array storing coordinates on PC space. This is genereated by the function "run".
       pca: PCA object of sklearn.decompositio.PCA that generates x
   Returns: 
       Nothing

   """
   with open("pca.dat", "w") as fout:
   
       fout.write("#id, pc1, pc2, pc3, pc4, pc5 \n")
       for iratio in pca.explained_variance_ratio_: 
           fout.write("#{} %\n".format(iratio * 100.0))
   
       for i, cods in enumerate(x, 1): 
#           cods_str = ["{:5.3f}".format(icod) for icod in cods ] 
#           fout.write("{:s} {} \n".format("MODEL"+str(i), cods_str))
           fout.write("{:s} {:5.3f} {:5.3f} {:5.3f} {:5.3f} {:5.3f}\n".format(
               "MODEL"+str(i), cods[0], cods[1], cods[2], cods[3], cods[4]))
   
       print("Shape of output: {}".format(x.shape))
 
if __name__ == "__main__":
    fileNameTraj, s1, s2, ref = argParse()

    if ref == None:
        x, pca = run(fileNameTraj, s1, s2)
    
    else: 
        x, pca = run(fileNameTraj, s1, s2, ref)

    write(x, pca)
