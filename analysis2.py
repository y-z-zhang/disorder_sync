# In[]
import os
import seaborn as sns
from scipy.signal import hilbert, find_peaks
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from matplotlib import rc
matplotlib.use('PDF')
rc('text', usetex=False)

#os.chdir('/home/yuanzhao/random_hetero/Istvan')
#print(os.getcwd())

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
m = 12  # number of experiments

# experimental data for coupled oscillators
data = np.zeros((m, 120480, n))
# homogeneous
data[0, :, :] = np.loadtxt('data/oc092019_6.txt')   # trial 1
data[1, :, :] = np.loadtxt('data/oc092019_15.txt')  # trial 2
data[2, :, :] = np.loadtxt('data/oc111919_22.txt')  # trial 6
data[3, :, :] = np.loadtxt('data/oc112519_16.txt')  # trial 11
data[4, :, :] = np.loadtxt('data/oc112519_39.txt')  # trial 14
data[5, :, :] = np.loadtxt('data/oc092019_21.txt')  # trial 4
# heterogeneous
data[6, :, :] = np.loadtxt('data/oc092019_8.txt')   # trial 1
data[7, :, :] = np.loadtxt('data/oc092019_10.txt')  # trial 2
data[8, :, :] = np.loadtxt('data/oc111919_24.txt')  # trial 6
data[9, :, :] = np.loadtxt('data/oc112519_19.txt')  # trial 11
data[10, :, :] = np.loadtxt('data/oc112519_42.txt')  # trial 14
data[11, :, :] = np.loadtxt('data/oc092019_19.txt')  # trial 4

# experimental data for uncoupled oscillators
Data = np.zeros((m, 20400, n))
# homogeneous
Data[0, :, :] = np.loadtxt('data/oc092019_3.txt')   # trial 1
Data[1, :, :] = np.loadtxt('data/oc092019_11.txt')  # trial 2
Data[2, :, :] = np.loadtxt('data/oc111919_19.txt')  # trial 6
Data[3, :, :] = np.loadtxt('data/oc112519_17.txt')  # trial 11
Data[4, :, :] = np.loadtxt('data/oc112519_40.txt')  # trial 14
Data[5, :, :] = np.loadtxt('data/oc092019_20.txt')  # trial 4
# heterogeneous
Data[6, :, :] = np.loadtxt('data/oc092019_7.txt')   # trial 1
Data[7, :, :] = np.loadtxt('data/oc092019_9.txt')   # trial 2
Data[8, :, :] = np.loadtxt('data/oc111919_23.txt')  # trial 6
Data[9, :, :] = np.loadtxt('data/oc112519_18.txt')  # trial 11
Data[10, :, :] = np.loadtxt('data/oc112519_43.txt')  # trial 14
Data[11, :, :] = np.loadtxt('data/oc092019_18.txt')  # trial 4

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

plt.plot(Delta[1:6], error[1:6], marker='o', ms=25, ls='none')
plt.text(.025, .05, 'homogeneous', fontsize=30, c=tableau20[0])
plt.plot(Delta[7:], error[7:], marker='o', ms=25, ls='none')
plt.text(.05, .02, 'heterogeneous', fontsize=30, c=tableau20[2])

#plt.xlim([.025, .075])
#plt.ylim([.01, .08])
plt.gca().tick_params(axis='y', pad=15, size=10, width=2)
plt.gca().tick_params(axis='x', pad=25, size=10, width=2)
#plt.xticks([.03, .05, .07])
#plt.yticks([.02, .05, .08])

fig.set_tight_layout(True)
plt.savefig('diagram.pdf')
