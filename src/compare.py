import sys

import numpy as np

from qutip import *
import matplotlib.pylab as plt

file1 = sys.argv[1]; file2 = sys.argv[2]

L1 = qload(file1); L2 = qload(file2)
L = L1-L2

print(np.trace(L2.data.todense()))
print(np.sqrt(np.trace((L.dag()*L).data.todense())))
print(np.sqrt(np.trace((L2.dag()*L2).data.todense())))

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.set_aspect('equal')

plt.imshow((L).data.todense().real, interpolation='nearest', cmap=plt.cm.ocean)
plt.colorbar()
plt.show()
