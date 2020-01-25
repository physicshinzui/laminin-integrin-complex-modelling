from MDAnalysis import Universe
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import pickle

def distances(obj1, obj2):
    """
    24.1.2020
    I aimed to compute all the distances between two protein chains.
    Afterwords, I'm gonna take the minimum distance of an c alpha atom to the other c alpha atoms.
    and then get the minimum distances.

    args:
        obj1: an object of selected_atoms genenrated from Universe
        obj2: another object of that
    return:
        distances: list, which consists of the minimum distances from iatom to other atoms.
            I wanted to get the minimium distances between c alpha atoms in two different interfaces.
    """
    distances = []
    iatom_distances = []
    for iatom in obj1:
        for jatom in obj2:
            iatom_distances.append(np.linalg.norm(iatom.position - jatom.position))
        distances.append(min(iatom_distances))
    return distances

def main():
    #24.1.2020
    #they were downloaded by a certain rule via the adcanced serarch in rcsb pdb,
    #but they contains more than 3 chains, which was out of my scope.
    omited_pdbs = ['4e7u.pdb', '4e7t.pdb', '3exx.pdb', '4fka.pdb', '5ep6.pdb',
                   '4gxv.pdb', '4uqp.pdb', '3uqy.pdb', '4dg4.pdb', '4urh.pdb',
                   '6f4j.pdb', '1xd3.pdb', '3bog.pdb', '6mee.pdb', '4pj2.pdb',
                   '5bpk.pdb', '3cjs.pdb', '4c2v.pdb', '1pid.pdb', '6fu9.pdb',
                   '2oxg.pdb', '1svf.pdb', '6fc1.pdb', '1q7l.pdb', '4kn9.pdb',
                   '4b2b.pdb', '6g6k.pdb', '4m4l.pdb', '4b2c.pdb', '1ben.pdb',
                   '3tt8.pdb', '3fq9.pdb', '5nwg.pdb', '4uql.pdb', '2xkn.pdb',
                   '5nwd.pdb']

    filenames = glob.glob('./interfaces/heterodimer/*.pdb')

    whole_distances = []
    for j, file in enumerate(filenames):
        #print(j, file, file.split('/')[-1] in omited_pdbs)

        if file.split('/')[-1] in omited_pdbs:
            continue

        else:
            print(file.split('/')[-1].split('.')[0].upper()+",",end = '')
            continue
            u = Universe(file)

            nchains = len(set(u.segments.segids))
            if nchains != 2: sys.exit(f'The number of chains = {nchains}. that is out of scope for this program.')

            chain_objs = []
            for i, chain in enumerate(set(u.segments.segids)):
                chain_objs.append(u.select_atoms(f'protein and segid {chain} and name CA'))
                print(f'    *{chain}')
#               print(chain_objs[i].atoms)

            whole_distances.append(distances(chain_objs[0], chain_objs[1]))

    sys.exit('stop! you might have already done this. so i forced you to procced.')
    print(whole_distances)
    whole_distances = np.hstack(whole_distances)
    filtered_dist   = whole_distances[whole_distances<=50.0]

    with open('whole_distances.pkl','wb') as f:
        pickle.dump(whole_distances, f)

    mean = np.mean(whole_distances)
    std  = np.std(whole_distances)

    print(f'mean: {mean}, std:{std}')

if __name__ == '__main__':
    main()
