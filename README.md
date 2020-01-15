# Integrin-laminin Complex modelling

## Outline
1. `01_gen_compl_candidates.py` generates los of PDBs, in which superimposition is carried out according to the C-alpha positions of the gamma tail in both laminin system and integrin-gtail one.

2. `02_calc_violation.py` contains a function that judge if an PDB of the PDBs generated is violated by considering C-alpha distances

FixME:
A problem is the computational expence.
if it is done, it will take 27 days!!!

# 9.1.2020
the number of conformations reduced to thoese in each cluster c1, c2, ....
at this time, the computatinal time was 2 hours.


# 10.1.2020
When I peformed pca, i did  
python pcaNatsu.py -i complex_models.pdb -s1 'segid C and resid 2981 3197 2832' -s2 '(segid A and resid 150 224) or (segid B and resid 166)'

maybe i should reduce the range of the nubmer of violaions...
    try: the number of violation (now) 10 to 50 -> 20 to 30
