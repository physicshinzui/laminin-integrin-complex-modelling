#!/Users/siida/anaconda3/bin/python
import numpy as np
from sys import exit
import os
from glob import glob
import MDAnalysis
from MDAnalysis import Universe
from MDAnalysis.analysis import align
from tqdm import tqdm
import logging

logging.basicConfig(filename='logging.log', level=logging.INFO)
logging.info(f'# snapshot numbers:')
logging.info(f'# Integrin-gtail, laminin')

"""
This program modells laminin-integrin complexes, based on MD trajectories
"""
int_gtail_ref_file  = './templates_traj/integrin_gtail.pdb'
int_gtail_xtc_file  = '../../data_to_be_analysed/traj/xtc_files/n_all_skipped_100.xtc'
lm_ref_file         = './templates_traj/laminin.pdb'
lm_xtc_file         = '../../../laminin/data_to_be_analysed/traj/all_skipped_10.xtc'

u_int_gtail = Universe(int_gtail_ref_file, int_gtail_xtc_file)
u_lm        = Universe(lm_ref_file, lm_xtc_file)

# Create selected objects for gtail in both trajectories.
gtail_on_integ = u_int_gtail.select_atoms('segid I and resid 1605-1609 and name CA')
gtail_in_lm    = u_lm.select_atoms('segid C and resid 1605-1609 and name CA')

snap_no_int = 0
for frame_int_g in tqdm(u_int_gtail.trajectory[0:2]):
    snap_no_lm = 0
    for frame_lm in u_lm.trajectory[0:2]:

        # The positions of the selected atoms above are centred.
        mobile = gtail_in_lm.positions    -    gtail_in_lm.atoms.center_of_mass()
        ref    = gtail_on_integ.positions - gtail_on_integ.atoms.center_of_mass()

        R, rmsd = align.rotation_matrix(mobile, ref)

        # Translate and rotate laminin according to the rotation matrix for gtail
        u_lm.atoms.translate(-gtail_in_lm.center_of_mass())
        u_lm.atoms.rotate(R)
        u_lm.atoms.translate(gtail_on_integ.center_of_mass())

        # Merge two systems laminin system and integrin-gtail one with C-alpha atoms
        u_merged = MDAnalysis.core.universe.Merge(u_lm.atoms.select_atoms('name CA'),
                                                  u_int_gtail.atoms.select_atoms('name CA'))
        u_merged.atoms.write(f"aligned_{snap_no_int}-{snap_no_lm}.pdb")
        logging.info(f'{snap_no_int} - {snap_no_lm}')
        snap_no_lm += 1

    snap_no_int += 1
