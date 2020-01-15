#!/Users/siida/anaconda3/bin/python
import numpy as np

def main():
    pcs = np.loadtxt('pca.dat', usecols=(1,2))

    conf_nums = []
    for i, pc in enumerate(pcs):
#        print(pc[0], pc[1])
#        if -45 < pc[0] < -35 and 0 < pc[1] < 5:
        if -5 < pc[0] < 5 and -2 < pc[1] < 2:
            print(f'MODEL NO - 1 = {i}: {pc}')
            conf_nums.append(f'{i:04}')
    print(conf_nums)

    path = '/Volumes/HD-siida/gtail_b1_sys/analysis/models'
    for i in conf_nums:
        cmd.load(f'{path}/{i}.pdb')


main()
