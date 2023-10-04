import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
# Example data
x = np.array([3, 10, 100, 200, 500, 1000, 1500, 2000, 4000, 6000])
y = np.array([float(0.0006011), float(0.003364), float(0.18128), float(0.29809),  float(0.3962), float(0.654), float(0.724), float(0.786), float(0.89), float(0.97)])

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
x_smooth = np.linspace(x.min(), x.max(), 1000, endpoint=True)
y_smooth = interp1d(x, y, kind='cubic')(x_smooth)
plt.ylim(0, 1)
plt.xlabel(r'\textit{Number of bodies}',fontsize=12)
plt.ylabel(r'\textit{Speedup}',fontsize=12)
# Make room for the ridiculously large title.
#plt.subplots_adjust(top=0.8)
plt.plot(x_smooth, y_smooth)


plt.show()