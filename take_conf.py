#!/Users/siida/anaconda3/bin/python
import sys
import numpy as np

def main():
    pcs = np.loadtxt('pca.dat', usecols = (1,2))

    for i, pc in enumerate(pcs,1):
        #if -30 < pc[0] <-10 and -10 < pc[1] < 10: # around 1st
        #if -30 < pc[0] < -10 and 10 < pc[1] < 40: # around 2nd
        #if 30 < pc[0] < 50 and -10 < pc[1] < 10:  # around 3rd
        if 0 < pc[0] < 20 and -25 < pc[1] < -10:  # around 

            print(i, pc)

if __name__ == "__main__":
    main()
