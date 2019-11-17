#!/Users/siida/anaconda3/bin/python
import numpy as np
from sys import exit
from glob import glob
from MDAnalysis import Universe
import MDAnalysis
from tqdm import tqdm
import scipy.spatial.distance as distance

def isViolation(univ, cutoff=6.0, lower_v=10, upper_v=20):
    """
    args:
        univ: Universe, which is an object in MDAnalysis
    returns:
        is_Violation: bool
    """
    distances = distance.cdist(univ.atoms.positions, univ.atoms.positions, metric='euclidean')
    #print(distances)
    nViolations = len(distances[distances<=cutoff])
    print(nViolations)

    if nViolations > lower_v and nViolations < upper_v:
        return False
    else:
        return True

if __name__ == "__main__":
    inps = glob("/Volumes/HD-siida/gtail_b1_sys/analysis/aligned_complexes/aligned_*")
    for i in inps:
        u = Universe(i)
#        u = Universe("aligned_0-8.pdb")
        if isViolation(u) == False:
            print(i)
