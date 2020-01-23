#!/Users/siida/anaconda3/bin/python
import sys
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
    model_no_w_scores = reada_vio('model_no.out')

    top_models = sorted(model_no_w_scores, key=lambda x:x[1], reverse=False)[0:100]
    top_model_numbers = [i[0] for i in top_models] # MODEL NO like a NMR model
    with open('rank.log','w') as flog:
        flog.write('#Rank, MODEL NO, Score\n')
        for i, elem in enumerate(top_models,1):
            flog.write(f'{i}, {elem[0]}, {elem[1]}\n')

    u = Universe("/Volumes/HD-siida/gtail_b1_sys/analysis/merged_aligned_complexes.pdb")
    models = u.select_atoms('all')
    with MDAnalysis.Writer("complex_models.pdb", models.n_atoms) as W:
        for model_no in tqdm(top_model_numbers):
            u.trajectory[model_no-1] # to position the uni to the address specified by model_no - 1
            #print(model_no -1)
            W.write(models) # the positioned uni was yielded here

main()
