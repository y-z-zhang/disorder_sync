# <md>
# This script analyzes how the level of synchrony changes with oscillator heterogeneity in our electrochemical experiments.

# In[]
from scipy.signal import hilbert, find_peaks
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from matplotlib import rc
#matplotlib.use('PDF')
rc('text', usetex=True)

# In[]
# These are the "Tableau 20" colors as RGB.
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)

# <md>
# Import time series for coupled oscillators measured from the electrochemical experiments (600 seconds of data for each time series).

# In[]
n = 16  # number of oscillators
m = 8   # number of experiments

# time series for coupled oscillators
data = np.zeros((m, 120480, n))
data[0, :, :] = np.loadtxt('data/oc091319_8.txt')
data[1, :, :] = np.loadtxt('data/oc110819_14.txt')
data[2, :, :] = np.loadtxt('data/oc091319_10.txt')
data[3, :, :] = np.loadtxt('data/oc091319_18.txt')
data[4, :, :] = np.loadtxt('data/oc091319_24.txt')
data[5, :, :] = np.loadtxt('data/oc111119_20.txt')
data[6, :, :] = np.loadtxt('data/oc111119_17.txt')
data[7, :, :] = np.loadtxt('data/oc091319_12.txt')


# <md>
# Calculate the time-averaged synchronization error of coupled oscillators for different heterogeneity

# In[]
# nominal oscillator heterogeneity introduced by setting the resistance of each oscillator to a different value
sigma = [0, .035, .07, .1, .13, .16, .19, .22]

# time-averaged synchronization error (averaged over the last 200 seconds of data for each experiment)
error = np.zeros(m)
for i in range(m):
    error[i] = np.mean(np.std(data[i, 80000:, :], axis=1))


# <md>
# Plot time-averaged synchronization error $⟨e⟩$ as a function of the nominal oscillator heterogeneity $\sigma$

# In[]
fig = plt.figure()

ax = plt.subplot(111)
for axis in ['top', 'right', 'bottom', 'left']:
    ax.spines[axis].set_linewidth(2)

fig.set_size_inches(9, 8)
plt.xticks(fontsize=35)
plt.yticks(fontsize=35)
plt.xlabel(r'$\sigma$', fontsize=40)
plt.ylabel(r'$\langle e \rangle$', fontsize=40)

plt.plot(sigma, error, marker='o', ms=35, ls='none', c=tableau20[8])

plt.gca().tick_params(axis='y', pad=5, size=10, width=2)
plt.gca().tick_params(axis='x', pad=5, size=10, width=2)
plt.xticks([0, .1, .2])
plt.yticks([.01, .04, .07])

fig.set_tight_layout(True)
plt.savefig('error.pdf')
