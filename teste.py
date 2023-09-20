import matplotlib.pyplot as plt
import numpy as np

D = np.random.randn(12*72).reshape((12, 72))
D[4, :] = np.nan
D[6, 6] = np.nan

D = np.ma.masked_invalid(D)

cmap = plt.get_cmap('bwr')
cmap.set_bad(color = 'k', alpha = 1.)

xbin = np.linspace(0, 12, 13)
ybin = np.linspace(-90, 90, 73)

fig = plt.figure()
ax = fig.add_subplot(111)
pl = ax.pcolormesh(xbin, ybin, D.T, cmap = cmap, edgecolors = 'None',
                vmin = -5, vmax = 5)
plt.show()
