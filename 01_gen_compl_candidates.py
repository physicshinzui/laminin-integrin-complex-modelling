#!/Users/siida/anaconda3/bin/python
import numpy as np
from sys import exit
import os
from glob import glob
import MDAnalysis
from MDAnalysis import Universe
from MDAnalysis.analysis import align
from tqdm import tqdm

"""
This program modells laminin-integrin complexes, based on MD trajectories
"""
int_gtail_ref_file  = './templates_traj/integrin_gtail.pdb'
int_gtail_xtc_file  = '../../data_to_be_analysed/traj/xtc_files/n_all_skipped_100.xtc'
lm_ref_file = './templates_traj/laminin.pdb'
lm_xtc_file = '../../../laminin/data_to_be_analysed/traj/all_skipped_10.xtc'

u_int_gtail = Universe(int_gtail_ref_file, int_gtail_xtc_file)
u_lm        = Universe(lm_ref_file, lm_xtc_file)

# Create selected objects for gtail in both trajectories.
gtail_on_integ = u_int_gtail.select_atoms('segid I and resid 1605-1609 and name CA')
gtail_in_lm    = u_lm.select_atoms('segid C and resid 1605-1609 and name CA')

#gtail_on_integ = u_int_gtail.select_atoms('segid E and resid 1605 and name CA')
#gtain_in_lm    = u_lm.select_atoms('segid C and resid 1605 and name CA')

for frame_int_g in tqdm(u_int_gtail.trajectory[0:1]):
    for frame_lm in u_lm.trajectory[0:1]:
        mobile = gtail_in_lm.positions    -    gtail_in_lm.atoms.center_of_mass()
        ref    = gtail_on_integ.positions - gtail_on_integ.atoms.center_of_mass()
#        print(gtain_in_lm.atoms)
#        print(gtail_on_integ.atoms)
#        print(mobile)
#        print(ref)
        R, rmsd = align.rotation_matrix(mobile, ref)
        print(rmsd)

        # translate and rotate laminin according to the rotation matrix for gtail
        u_lm.atoms.translate(-gtail_in_lm.center_of_mass())
        u_lm.atoms.rotate(R)
        u_lm.atoms.translate(gtail_on_integ.center_of_mass())
        u_merged = MDAnalysis.core.universe.Merge(u_lm.atoms, u_int_gtail.atoms)
        u_merged.atoms.write("mobile_on_ref.pdb")
#        u_int_gtail.atoms.write("mobile_on_a.pdb")
