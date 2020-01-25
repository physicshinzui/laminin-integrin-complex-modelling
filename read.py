import pickle
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gamma

with open('whole_distances.pkl','rb') as fin:
    distances = pickle.load(fin)

filtered_dist = distances[distances < 10.0]
#mean = np.mean(filtered_dist)
#std  = np.std(filtered_dist)
mean = np.mean(distances)
std  = np.std(distances)
print(f'mean: {mean}, std:{std}')

a = 1.99323054838
mean, var, skew, kurt = gamma.stats(a, moments='mvsk')

plt.hist(distances, bins=np.arange(0,max(distances),0.5), color='red')
plt.xlim(0,)
plt.show()
