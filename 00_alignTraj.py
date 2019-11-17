#!/Users/siida/anaconda3/bin/python
import numpy as np
from sys import exit
import os
from glob import glob
from MDAnalysis import Universe
from MDAnalysis.analysis import align

"""
This program modells laminin-integrin complexes, based on MD trajectories
"""
ref_file  = '../../data_to_be_analysed/npt2_pdbnumberBase_wo_solv.pdb'
xtc_file  = '../../data_to_be_analysed/traj/xtc_files/n_all_skipped_100.xtc'
prefix_out_file = xtc_file.split('/')[-1].split('.')[0]
ref  = Universe(ref_file)
traj = Universe(ref_file, xtc_file)
alignment = align.AlignTraj(traj, ref,
                            select='resid 122-361 and name CA',
                            filename=f'rmsfit_{prefix_out_file}.xtc',
                            verbose=True)
alignment.run()
