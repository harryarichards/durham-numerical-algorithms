import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)
y1, y2, y3, y4, y5, y6, y7, y8 = [], [], [], [], [], [], [], []
num_lines = 0
with open('1e-1.txt') as f:
    for line in f:
        line = line[-26:]
        num = float(line.strip('dx_min='))
        y1.append(num)
with open('1e-2.txt') as f:
    for line in f:
        line = line[-26:]
        num = float(line.strip('dx_min='))
        y2.append(num)
with open('1e-3.txt') as f:
    for line in f:
        line = line[-26:]
        num = float(line.strip('dx_min='))
        y3.append(num)
with open('1e-4.txt') as f:
    for line in f:
        line = line[-26:]
        num = float(line.strip('dx_min='))
        y4.append(num)
with open('1e-5.txt') as f:
    for line in f:
        line = line[-26:]
        num = float(line.strip('dx_min='))
        y5.append(num)
with open('1e-6.txt') as f:
    for line in f:
        line = line[-26:]
        num = float(line.strip('dx_min='))
        y6.append(num)
with open('1e-7.txt') as f:
    for line in f:
        line = line[-26:]
        num = float(line.strip('dx_min='))
        y7.append(num)
with open('1e-8.txt') as f:
    for line in f:
        line = line[-26:]
        num = float(line.strip('dx_min='))
        y8.append(num)
x_adaptive = []
adaptive = []
with open('adaptive.txt') as f:
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
y1_smooth = interp1d(x, y1, kind='cubic')(x_smooth)
y2_smooth = interp1d(x, y2, kind='cubic')(x_smooth)
y3_smooth = interp1d(x, y3, kind='cubic')(x_smooth)
y4_smooth = interp1d(x, y4, kind='cubic')(x_smooth)
y5_smooth = interp1d(x, y5, kind='cubic')(x_smooth)
y6_smooth = interp1d(x, y6, kind='cubic')(x_smooth)
y7_smooth = interp1d(x, y7, kind='cubic')(x_smooth)
y8_smooth = interp1d(x, y8, kind='cubic')(x_smooth)
x_adaptive_smooth = np.linspace(x_adaptive.min(), x_adaptive.max(), 300, endpoint=True)
adaptive_smooth = interp1d(x_adaptive, adaptive, kind='cubic')(x_adaptive_smooth)

plt.xlabel(r'\textit{time} (t)',fontsize=12)
plt.ylabel(r'\textit{Distance}',fontsize=12)
# Make room for the ridiculously large title.
#plt.subplots_adjust(top=0.8)
plt.plot(x_smooth, y1_smooth,  '-b', label=r'$10^{-1}$')
plt.plot(x_smooth, y2_smooth,  '-r', label=r'$10^{-2}$')
plt.plot(x_smooth, y3_smooth,  '-g', label=r'$10^{-3}$')
plt.plot(x_smooth, y4_smooth,  '-c', label=r'$10^{-4}$')
plt.plot(x_smooth, y5_smooth,  '-k', label=r'$10^{-5}$')
plt.plot(x_smooth, y6_smooth,  'chartreuse', label=r'$10^{-6}$')
plt.plot(x_smooth, y7_smooth,  '-y', label=r'$10^{-7}$')
plt.plot(x_smooth, y8_smooth,  '-m', label=r'$10^{-8}$')
plt.plot(x_adaptive_smooth, adaptive_smooth,  'purple', label=r'\textit{adaptive}')

plt.legend(title= 'time step', loc='upper left')
plt.show()