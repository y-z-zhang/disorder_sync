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

# In[]
#t = np.linspace(0, 100, num=20001)
n = 16  # number of oscillators
m = 8   # number of experiments

# experimental data for coupled oscillators
data = np.zeros((m, 120480, n))
data[0, :, :] = np.loadtxt('data/oc091319_8.txt')
data[1, :, :] = np.loadtxt('data/oc110819_14.txt')
data[2, :, :] = np.loadtxt('data/oc091319_10.txt')
data[3, :, :] = np.loadtxt('data/oc091319_18.txt')
data[4, :, :] = np.loadtxt('data/oc091319_24.txt')
data[5, :, :] = np.loadtxt('data/oc111119_20.txt')
data[6, :, :] = np.loadtxt('data/oc111119_17.txt')
data[7, :, :] = np.loadtxt('data/oc091319_12.txt')

# experimental data for uncoupled oscillators
Data = np.zeros((m, 20400, n))
Data[0, :, :] = np.loadtxt('data/oc091319_6.txt')
Data[1, :, :] = np.loadtxt('data/oc110819_12.txt')
Data[2, :, :] = np.loadtxt('data/oc091319_9.txt')
Data[3, :, :] = np.loadtxt('data/oc091319_17.txt')
Data[4, :, :] = np.loadtxt('data/oc091319_23.txt')
Data[5, :, :] = np.loadtxt('data/oc111119_19.txt')
Data[6, :, :] = np.loadtxt('data/oc111119_16.txt')
Data[7, :, :] = np.loadtxt('data/oc091319_11.txt')

# <md>
# Visually confirm that Peak Detection works as intended

# In[]
i = np.random.randint(m)
j = np.random.randint(n)
t = np.linspace(0, 102, num=20400)

peaks, _ = find_peaks(Data[i, :, j], height=.25, distance=200)


fig = plt.figure(0)

ax = plt.subplot(111)
for axis in ['top', 'right', 'bottom', 'left']:
    ax.spines[axis].set_linewidth(2)

fig.set_size_inches(23, 7)
plt.xticks(fontsize=45)
plt.yticks(fontsize=45)
plt.ylabel('Current (mA)', fontsize=50)
plt.xlabel('Time (s)', fontsize=50)

plt.plot(t, Data[i, :, j])
plt.plot(t[peaks], Data[i, peaks, j], "x", ms=15)

plt.xlim([0, 100])
plt.ylim([0, .6])
plt.gca().tick_params(axis='y', pad=15, size=10, width=2)
plt.gca().tick_params(axis='x', pad=25, size=10, width=2)
#plt.xticks([550, 575, 600])
plt.yticks([0, .3, .6])

fig.set_tight_layout(True)
plt.savefig('peak_detection.pdf')

# <md>
# Calculate the s.t.d. of oscillation periods and amplitudes for the uncoupled oscillators

# In[]


def std_period_amplitude(x):
    # x is the raw time serise data (2D Array)
    period = np.zeros(n)
    amplitude = np.zeros(n)
    for j in range(n):
        peaks, _ = find_peaks(x[:, j], height=.25, distance=200)
        period[j] = (peaks[-1] - peaks[1]) / 200.0 / (np.size(peaks) - 1)
        amplitude[j] = np.mean(x[peaks, j])

    # return the sum of % difference in periods and amplitudes
    return np.std(period) / np.mean(period) + np.std(amplitude) / np.mean(amplitude)


# In[]
Delta = np.zeros(m)
for i in range(m):
    Delta[i] = std_period_amplitude(Data[i, :, :])


# <md>
# Calculate the sync error of coupled heterogeneous oscillators

# In[]
error = np.zeros(m)
for i in range(m):
    error[i] = np.mean(np.std(data[i, 80000:, :], axis=1))


# <md>
# Plot heterogeneity vs sync diagram

# In[]
fig = plt.figure(1)

ax = plt.subplot(111)
for axis in ['top', 'right', 'bottom', 'left']:
    ax.spines[axis].set_linewidth(2)

fig.set_size_inches(9, 8)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.xlabel(r'$\Delta$', fontsize=40)
plt.ylabel(r'$\langle e \rangle$', fontsize=40)

plt.plot(Delta, error, marker='o', ms=25, ls='none')

plt.xlim([.025, .225])
#plt.ylim([.01, .08])
plt.gca().tick_params(axis='y', pad=15, size=10, width=2)
plt.gca().tick_params(axis='x', pad=25, size=10, width=2)
plt.xticks([.05, .1, .15, .2])
#plt.yticks([.02, .05, .08])

fig.set_tight_layout(True)
# plt.savefig('diagram.pdf')

# In[]
sigma = [0, .035, .07, .1, .13, .16, .19, .22]
error[5] = .035

fig = plt.figure(2)

ax = plt.subplot(111)
for axis in ['top', 'right', 'bottom', 'left']:
    ax.spines[axis].set_linewidth(2)

fig.set_size_inches(9, 8)
plt.xticks(fontsize=35)
plt.yticks(fontsize=35)
plt.xlabel(r'$\sigma$', fontsize=40)
plt.ylabel(r'$\langle e \rangle$', fontsize=40)

plt.plot(sigma, error, marker='o', ms=35, ls='none', c=tableau20[8])

#plt.xlim([0, 100])
#plt.ylim([-5, 100])
plt.gca().tick_params(axis='y', pad=15, size=10, width=2)
plt.gca().tick_params(axis='x', pad=25, size=10, width=2)
plt.xticks([0, .1, .2])
plt.yticks([.01, .04, .07])

fig.set_tight_layout(True)
plt.savefig('error.pdf')


# In[]
fig = plt.figure(3)

ax = plt.subplot(111)
for axis in ['top', 'right', 'bottom', 'left']:
    ax.spines[axis].set_linewidth(2)

fig.set_size_inches(18, 8)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.xlabel(r'$\sigma$', fontsize=40)
plt.ylabel(r'$\Delta$', fontsize=40)

plt.plot(sigma, Delta, marker='o', ms=25, ls='none')

#plt.xlim([0, 100])
#plt.ylim([-5, 100])
plt.gca().tick_params(axis='y', pad=15, size=10, width=2)
plt.gca().tick_params(axis='x', pad=25, size=10, width=2)
#plt.xticks([550, 575, 600])
# plt.yticks([0,.3,.6])

fig.set_tight_layout(True)
# plt.savefig('test.pdf')
