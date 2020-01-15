#!/Users/siida/anaconda3/bin/python
import numpy as np
from sys import exit
import os
from glob import glob
import argparse
from MDAnalysis import Universe
import MDAnalysis
from tqdm import tqdm
import scipy.spatial.distance as distance
from pprint import pprint
#import math

def main():
    u = Universe('/Volumes/HD-siida/gtail_b1_sys/analysis/merged_aligned_complexes.pdb')
    #u = Universe('complex_models.pdb')
    print(u.atoms.segids)
    ca_integrinAB = u.select_atoms('segid A B and name CA')
    ca_lamininE8  = u.select_atoms('segid C D E and name CA')
    #pprint([''.join(resn+str(resi)) for resn,resi in zip(ca_lamininE8.resnames, ca_lamininE8.resids)])

    lower, upper = 6.0, 10.0
    #cutoff = 2.75 # water molecule diameter
    with open('model_no.out','w') as fout:
        #fout.write(f'#MODEL NO, nViolations (if distances<={cutoff}) \n')
        fout.write(f'#MODEL NO, nViolations (if r<{lower}), nContacts ({lower}<=r<={upper}) \n')

        for i, frame in enumerate(tqdm(u.trajectory)):
            distances = distance.cdist(ca_integrinAB.positions,ca_lamininE8.positions, metric='euclidean')
            #nViolations = len(distances[distances<=cutoff])
            nViolations = len(distances[distances<lower])
            nContacts   = len(distances[(distances<=upper) & (distances>=lower)])
#            score = nContacts -nViolations
            score = np.log(nContacts/nViolations)
            fout.write(f'{i}, {nViolations}, {nContacts}, {score}\n')
            #if nViolations >= 50:
            #    fout.write(f'# {i}; Too much violations   = {nViolations}\n')
            #else:
            #    fout.write(f'{i} # Acceptable violations = {nViolations}\n')

if __name__ == "__main__":
    main()
