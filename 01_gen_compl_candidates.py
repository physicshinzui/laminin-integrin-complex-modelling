#!/Users/siida/anaconda3/bin/python
import numpy as np
import sys
import os
from glob import glob
import MDAnalysis
from MDAnalysis import Universe
from MDAnalysis.analysis import align
from tqdm import tqdm
import logging

def rmsd_fit(mobile_obj, ref_obj, u_mobile):
    """
    this changes state of u_mobile
    """
    # The positions of the selected atoms above are centred.
    mobile = mobile_obj.positions - mobile_obj.atoms.center_of_mass()
    ref    = ref_obj.positions    - ref_obj.atoms.center_of_mass()
    #print(mobile_obj.resids)
    #print(ref_obj.resids)
    R, rmsd = align.rotation_matrix(mobile, ref)
    #print(R, rmsd)
    # Translate and rotate mobile_obj according to the rotation matrix
    u_mobile.atoms.translate(-mobile_obj.center_of_mass())
    u_mobile.atoms.rotate(R)
    u_mobile.atoms.translate(ref_obj.center_of_mass())

#            # The positions of the selected atoms above are centred.
#            mobile = gtail_in_lm.positions    -    gtail_in_lm.atoms.center_of_mass()
#            ref    = gtail_on_integ.positions - gtail_on_integ.atoms.center_of_mass()
#            R, rmsd = align.rotation_matrix(mobile, ref)
#
#            # Translate and rotate laminin according to the rotation matrix for gtail
#            u_lm.atoms.translate(-gtail_in_lm.center_of_mass())
#            u_lm.atoms.rotate(R)
#            u_lm.atoms.translate(gtail_on_integ.center_of_mass())

def main():
    logging.basicConfig(filename='logging.log', level=logging.CRITICAL)
    logging.info(f'# snapshot numbers:')
    logging.info(f'# Integrin-gtail, laminin')

    """
    This program modells laminin-integrin complexes, based on MD trajectories
    """
    # confomers in each cluster: integrin and gtail system
#    int_gtail_ref_file  = './templates_traj/ca_integrin_gtail.pdb'
    int_gtail_ref_file  = './templates_traj/ca_integrin_gtail_renum.pdb'
    int_gtail_conformers_xtc_file = '../1dPMF/conformers/conformers_in_all.xtc'

    # snapshots of lamininE8
    lm_ref_file         = './templates_traj/laminin.pdb'
    lm_xtc_file         = '../../../laminin/data_to_be_analysed/traj/all_skipped_10.xtc'

    # X-ray structures: LM-E8 and integrin a5b1
    lme8     = 'orig_laminin/monomer_lam_remv_anisou.pdb'
    integrin = 'orig_integin/4wjk_remv_anisou.pdb'

    # universe objects generation
    u_int_gtail = Universe(int_gtail_ref_file, int_gtail_conformers_xtc_file)
    u_lm        = Universe(lm_ref_file, lm_xtc_file)
    u_xray_lme8  = Universe(lme8)
    u_xray_integ = Universe(integrin)

    # Create selected objects for gtail in both trajectories.
    gtail_on_integ = u_int_gtail.select_atoms('segid I and resid 1605-1609 and name CA')
    gtail_in_lm    = u_lm.select_atoms('segid C and resid 1605-1609 and name CA')

    mobile_reg_xray_lme8 = u_xray_lme8.select_atoms('resid 2719-2770 or resid 2779-3121 and name CA')
    ref_reg_gtail_in_lm  = u_lm.select_atoms('resid 2719-2770 or resid 2779-3121 and name CA')
    #print(len(mobile_reg_xray_lme8.residues), mobile_reg_xray_lme8.resids)
    #print(len(ref_reg_gtail_in_lm.residues), ref_reg_gtail_in_lm.resids)

    mobile_reg_xray_integ = u_xray_integ.select_atoms('segid B and resid 123-360 and name CA')
    ref_reg_gtail_on_integ= u_int_gtail.select_atoms('segid E and resid 123-360 and name CA')
#    print(len(mobile_reg_xray_integ.residues))
#    print(len(ref_reg_gtail_on_integ.residues))

#    sys.exit()

    snap_no_int = 0
    for frame_int_g in tqdm(u_int_gtail.trajectory):
        snap_no_lm = 0
        for frame_lm in tqdm(u_lm.trajectory):
#        for frame_lm in tqdm(u_lm.trajectory[0:1]):

            #(1) fit gtail_in_lm to gtail_on_integ
            rmsd_fit(gtail_in_lm, gtail_on_integ, u_lm)

            #(2) fit
            rmsd_fit(mobile_reg_xray_lme8, ref_reg_gtail_in_lm, u_xray_lme8)

            # (3) fit integrin Xray to u_int_gtail
            rmsd_fit(mobile_reg_xray_integ, ref_reg_gtail_on_integ, u_xray_integ)

            # Merge two systems laminin system and integrin-gtail one with C-alpha atoms
#            u_merged = MDAnalysis.core.universe.Merge(u_lm.atoms.select_atoms('name CA'),
#                                                      u_int_gtail.atoms.select_atoms('name CA'),
#                                                      u_xray_integ.atoms.select_atoms('name CA'),
#                                                      u_xray_lme8.atoms.select_atoms('name CA'))
#
            u_merged = MDAnalysis.core.universe.Merge(u_xray_integ.atoms.select_atoms('name CA'),
                                                      u_xray_lme8.atoms.select_atoms('name CA'))
            #u_merged.atoms.write(f"complexes/aligned_{snap_no_int}-{snap_no_lm}.pdb")
            u_merged.atoms.write(f"/Volumes/HD-siida/gtail_b1_sys/analysis/aligned_complexes/aligned_{snap_no_int}-{snap_no_lm}.pdb")

            logging.info(f'{snap_no_int} - {snap_no_lm}')
            snap_no_lm += 1

        snap_no_int += 1

if __name__ == '__main__':
    main()
