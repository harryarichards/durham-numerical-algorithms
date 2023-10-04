import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
y1, y2, y3, y4, y5, y6, y7, y8 = [], [], [], [], [], [], [], []
num_lines = 0
with open('5bodye-5.txt') as f:
    for line in f:
        line = line[-26:]
        num = float(line.strip('dx_min='))
        y5.append(num)
with open('5bodye-4.txt') as f:
    for line in f:
        line = line[-26:]
        num = float(line.strip('dx_min='))
        y4.append(num)
with open('5bodye-3.txt') as f:
    for line in f:
        line = line[-26:]
        num = float(line.strip('dx_min='))
        y3.append(num)


x_adaptive = []
adaptive = []
with open('5bodyadaptive.txt') as f:
    for line1 in f:
        line = line1[-26:]
        num = float(line.strip('dx_min='))
        adaptive.append(num)

        line = line1[45:57]
        line = line.strip('t=, ')
        num = float(line.strip('t=,'))
        x_adaptive.append(num)

x_adaptive = np.array(x_adaptive)
x = np.arange(0, 10001, 100)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
x_smooth = np.linspace(x.min(), x.max(), 300, endpoint=True)
y3_smooth = interp1d(x, y3, kind='cubic')(x_smooth)
y4_smooth = interp1d(x, y4, kind='cubic')(x_smooth)
y5_smooth = interp1d(x, y5, kind='cubic')(x_smooth)
x_adaptive_smooth = np.linspace(x_adaptive.min(), x_adaptive.max(), 300, endpoint=True)
adaptive_smooth = interp1d(x_adaptive, adaptive, kind='cubic')(x_adaptive_smooth)

plt.xlabel(r'\textit{time} (t)',fontsize=12)
plt.ylabel(r'\textit{Distance}',fontsize=12)
# Make room for the ridiculously large title.
#plt.subplots_adjust(top=0.8)
plt.plot(x_smooth, y3_smooth,  '-b', label=r'$10^{-3}$')
plt.plot(x_smooth, y4_smooth,  '-r', label=r'$10^{-4}$')
plt.plot(x_smooth, y5_smooth,  '-g', label=r'$10^{-5}$')
plt.plot(x_adaptive_smooth, adaptive_smooth,  'purple', label=r'\textit{adaptive}')

plt.legend(title= 'time step', loc='upper left')
plt.show()