import matplotlib.pyplot as plt
import numpy as np

pc = np.loadtxt('pca.dat', usecols = (1,4))
plt.scatter(pc[:,0], pc[:,1], s=3.0)
plt.show()
