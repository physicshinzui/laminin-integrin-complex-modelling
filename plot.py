import matplotlib.pyplot as plt
import numpy as np

pc = np.loadtxt('pca.dat', usecols = (1,2))
colours = [i for i in range(len(pc[:,0]))]
#plt.scatter(pc[:,0], pc[:,1], s=20.0, c=colours, cmap = 'gnuplot')
plt.scatter(pc[:,0], pc[:,1], s=50.0)

plt.show()
