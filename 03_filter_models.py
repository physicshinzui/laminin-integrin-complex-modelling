#!/Users/siida/anaconda3/bin/python
import sys
import os
from glob import glob
import random
from MDAnalysis import Universe
import MDAnalysis
from tqdm import tqdm

def reada_vio(violationFile):
    with open(violationFile,'r') as fin:
        scores = []
        for line in fin:
            if not line.startswith('#'):
                row = line.split(',')
                model_number  = int(row[0].strip())
                nViolations   = int(row[1].strip())
                nContacts     = int(row[2].strip())
                score         = float(row[3].strip())
                scores.append([model_number, score])
    return scores

def pymol_process(model_numbers):
    for number in model_numbers:
        cmd.load(f'complexes/{number}.pdb')

def main():
#    model_numbers, nViolations_of_models = reada_vio('model_no.out')
    model_no_w_scores = reada_vio('model_no.out')
    lower_bound, upper_bound = 50, 60

#    filtered_model_numbers = [ model_number for model_number, n_violations in zip(model_numbers, nViolations_of_models)
#                               if n_violations <= upper_bound and n_violations >= lower_bound]
#    print(filtered_model_numbers, len(filtered_model_numbers))
#    print(model_no_w_scores)
    top_models = sorted(model_no_w_scores, key=lambda x:x[1], reverse=True)[0:10]
    top_model_numbers = [i[0] for i in top_models]
    print(top_models)
    print(top_model_numbers)

#    sys.exit()
    u = Universe("/Volumes/HD-siida/gtail_b1_sys/analysis/merged_aligned_complexes.pdb")
#    u = Universe("complex_models.pdb")
    models = u.select_atoms('all')
    with MDAnalysis.Writer("complex_models.pdb", models.n_atoms) as W:
        for i, frame in tqdm(enumerate(u.trajectory)):
            if i in top_model_numbers:
                W.write(models)

main()
