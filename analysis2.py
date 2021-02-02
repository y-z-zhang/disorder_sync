# %% markdown
# ## This script compares homogeneous and heterogeneous oscillators in terms of the measured oscillator heterogeneity (when uncoupled) and the time-averaged synchronization error (when coupled).
# ## Five independent sets of electrochemical experiments were performed, which demonstrate that more heterogeneous oscillators can synchronize better (even when heterogeneity is random).

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

# %% markdown
# Import time series for coupled oscillators measured from the electrochemical experiments.

# In[]
n = 16  # number of oscillators
m = 12  # number of experiments

# experimental time series for coupled oscillators (600 seconds of data for each time series)
data = np.zeros((m, 120480, n))
# homogeneous oscillators
data[0, :, :] = np.loadtxt('data/oc092019_6.txt')
data[1, :, :] = np.loadtxt('data/oc092019_15.txt')
data[2, :, :] = np.loadtxt('data/oc111919_22.txt')
data[3, :, :] = np.loadtxt('data/oc112519_16.txt')
data[4, :, :] = np.loadtxt('data/oc112519_39.txt')
data[5, :, :] = np.loadtxt('data/oc092019_21.txt')
# heterogeneous oscillators
data[6, :, :] = np.loadtxt('data/oc092019_8.txt')
data[7, :, :] = np.loadtxt('data/oc092019_10.txt')
data[8, :, :] = np.loadtxt('data/oc111919_24.txt')
data[9, :, :] = np.loadtxt('data/oc112519_19.txt')
data[10, :, :] = np.loadtxt('data/oc112519_42.txt')
data[11, :, :] = np.loadtxt('data/oc092019_19.txt')

# experimental time series for uncoupled oscillators (200 seconds of data for each time series)
Data = np.zeros((m, 20400, n))
# homogeneous oscillators
Data[0, :, :] = np.loadtxt('data/oc092019_3.txt')
Data[1, :, :] = np.loadtxt('data/oc092019_11.txt')
Data[2, :, :] = np.loadtxt('data/oc111919_19.txt')
Data[3, :, :] = np.loadtxt('data/oc112519_17.txt')
Data[4, :, :] = np.loadtxt('data/oc112519_40.txt')
Data[5, :, :] = np.loadtxt('data/oc092019_20.txt')
# heterogeneous oscillators
Data[6, :, :] = np.loadtxt('data/oc092019_7.txt')
Data[7, :, :] = np.loadtxt('data/oc092019_9.txt')
Data[8, :, :] = np.loadtxt('data/oc111919_23.txt')
Data[9, :, :] = np.loadtxt('data/oc112519_18.txt')
Data[10, :, :] = np.loadtxt('data/oc112519_43.txt')
Data[11, :, :] = np.loadtxt('data/oc092019_18.txt')

# %% markdown
# Visually confirm that Peak Detection works as intended.
# Peak Detection is used to extract the frequency and amplitude of uncoupled oscillators, which in turn are used to calculate the measured oscillator heterogeneity.

# In[]
# randomly pick an oscillator from one of the time series
i = np.random.randint(m)
j = np.random.randint(n)
t = np.linspace(0, 102, num=20400)

# identify all the peaks
peaks, _ = find_peaks(Data[i, :, j], height=.25, distance=200)

# plot the current for that oscillator with all the peaks marked by x
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
plt.gca().tick_params(axis='y', pad=5, size=10, width=2)
plt.gca().tick_params(axis='x', pad=5, size=10, width=2)
plt.yticks([0, .3, .6])

fig.set_tight_layout(True)
plt.savefig('peak_detection.pdf')

# %% markdown
# Calculate the standard deviation of oscillation periods and amplitudes for the uncoupled oscillators.
# These standard deviations are then combined to calculate the measured oscillator heterogeneity.

# In[]
def std_period_amplitude(x):
    # x is the raw time serise data (2D Array)
    period = np.zeros(n)
    amplitude = np.zeros(n)
    for j in range(n):
        peaks, _ = find_peaks(x[:, j], height=.25, distance=200)
        period[j] = (peaks[-1] - peaks[1]) / 200.0 / (np.size(peaks) - 1)
        amplitude[j] = np.mean(x[peaks, j])

    # return the sum of % difference in periods and amplitudes (i.e., measured oscillator heterogeneity)
    return np.std(period) / np.mean(period) + np.std(amplitude) / np.mean(amplitude)


# In[]
# Delta is the measured oscillator heterogeneity
Delta = np.zeros(m)
for i in range(m):
    Delta[i] = std_period_amplitude(Data[i, :, :])


# %% markdown
# Calculate the time-averaged synchronization error of coupled oscillators for different realizations of heterogeneity (one realization per experiment).

# In[]
# time-averaged synchronization error (averaged over the last 200 seconds of data for each experiment)
error = np.zeros(m)
for i in range(m):
    error[i] = np.mean(np.std(data[i, 80000:, :], axis=1))


# %% markdown
# Plot time-averaged synchronization error $⟨e⟩$ vs. measured oscillator heterogeneity $\Delta$, where each dot represents a different realization of heterogeneous (orange) and homogeneous (blue) systems.

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

# plot homogeneous systems
plt.plot(Delta[1:6], error[1:6], marker='o', ms=25, ls='none')
plt.text(.025, .05, 'homogeneous', fontsize=30, c=tableau20[0])
# plot heterogeneous systems
plt.plot(Delta[7:], error[7:], marker='o', ms=25, ls='none')
plt.text(.05, .02, 'heterogeneous', fontsize=30, c=tableau20[2])

#plt.xlim([.025, .075])
#plt.ylim([.01, .08])
plt.gca().tick_params(axis='y', pad=5, size=10, width=2)
plt.gca().tick_params(axis='x', pad=5, size=10, width=2)
#plt.xticks([.03, .05, .07])
#plt.yticks([.02, .05, .08])

fig.set_tight_layout(True)
plt.savefig('diagram.pdf')
