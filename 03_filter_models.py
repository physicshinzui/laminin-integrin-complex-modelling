#!/Users/siida/anaconda3/bin/python
from sys import exit
import os
from glob import glob
import random

def reada_vio(violationFile):
    with open(violationFile,'r') as fin:
        model_numbers, nViolations_of_models = [], []
        for line in fin:
            if not line.startswith('#'):
                row = line.split(',')
                model_numbers.append(int(row[0].strip()))
                nViolations_of_models.append(int(row[1].strip()))
    return model_numbers, nViolations_of_models

def pymol_process(model_numbers):
    for number in model_numbers:
        cmd.load(f'complexes/{number}.pdb')
                
def main():
    model_numbers, nViolations_of_models = reada_vio('model_no.out')
    upper_bound, lower_bound = 20, 10
    filtered_model_numbers = [ model_number for model_number, n_violations in zip(model_numbers, nViolations_of_models) 
                               if n_violations <= upper_bound and n_violations >= lower_bound] 
    print(filtered_model_numbers)

    pymol_process(random.sample(filtered_model_numbers, 10))

main()
