#!/Users/siida/anaconda3/bin/python
import numpy as np
from glob import glob
import argparse
from MDAnalysis import Universe
import MDAnalysis
from tqdm import tqdm
import scipy.spatial.distance as distance

def main():
    u = Universe('/Volumes/HD-siida/gtail_b1_sys/analysis/merged_aligned_complexes.pdb')
    #u = Universe('complex_models.pdb')
    print(u.atoms.segids)
    ca_integrinAB = u.select_atoms('segid A B and name CA')
    ca_lamininE8  = u.select_atoms('segid C D E and name CA')

    lower, upper = 6.0, 10.0
    with open('model_no.out','w') as fout:
        fout.write(f'#MODEL NO, nViolations (if r<{lower}), nContacts ({lower}<=r<={upper}) \n')

        for i, frame in enumerate(tqdm(u.trajectory),1): # Note that i starts with 1.
            distances = distance.cdist(ca_integrinAB.positions,ca_lamininE8.positions, metric='euclidean')
            #nViolations = len(distances[distances<=cutoff])
            nViolations = len(distances[distances<lower])
            nContacts   = len(distances[(distances<=upper) & (distances>=lower)])
#            score = nContacts -nViolations
            if nViolations != 0:
                score = -0.59 * np.log(nContacts/nViolations)
            else:
                score = np.nan
            fout.write(f'{i}, {nViolations}, {nContacts}, {score}\n')

if __name__ == "__main__":
    main()
