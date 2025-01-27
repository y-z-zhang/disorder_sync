# %% markdown
# This script plots the time series of oscillators measured from the experiments.

# In[]
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
#matplotlib.use('PDF')
sfont = {'fontname': 'serif'}

# %% markdown
# Import time series for coupled oscillators measured from the electrochemical experiments (600 seconds of data for each time series)

# In[]
n = 16 # number of oscillators
t = np.linspace(0, 600, num=120001)

# time series for homogeneous and heterogeneous oscillators, respectively (Fig.7)
current_homo = np.loadtxt('data/oc091319_8.txt')
current_hetero = np.loadtxt('data/oc091319_24.txt')
# time series for homogeneous and heterogeneous oscillators, respectively (Fig.S7)
#current_homo = np.loadtxt('data/oc021121_24.txt')
#current_hetero = np.loadtxt('data/oc021121_26.txt')

# %% markdown
# Plot the evolution of currents

# In[]
fig = plt.figure(1)

ax = plt.subplot(111)
for axis in ['top', 'right', 'bottom', 'left']:
    ax.spines[axis].set_linewidth(2)

fig.set_size_inches(23, 7)
plt.xticks(fontsize=45)
plt.yticks(fontsize=45)
plt.ylabel('Current (mA)', fontsize=50, **sfont)
plt.xlabel('Time (s)', fontsize=50, **sfont)

for i in range(n):
    plt.plot(t, current_homo[:120001, i], ls='-', lw=2)

plt.xlim([550, 600])
plt.ylim([0, .6])
plt.gca().tick_params(axis='y', pad=15, size=10, width=2)
plt.gca().tick_params(axis='x', pad=25, size=10, width=2)
#plt.xticks([550, 575, 600])
plt.yticks([0, .3, .6])

fig.set_tight_layout(True)
plt.savefig('trj_homo.pdf')


# In[]
fig = plt.figure(2)

ax = plt.subplot(111)
for axis in ['top', 'right', 'bottom', 'left']:
    ax.spines[axis].set_linewidth(2)

fig.set_size_inches(23, 7)
plt.xticks(fontsize=45)
plt.yticks(fontsize=45)
plt.ylabel('Current (mA)', fontsize=50, **sfont)
plt.xlabel('Time (s)', fontsize=50, **sfont)

for i in range(n):
    plt.plot(t, current_hetero[:120001, i], ls='-', lw=2)

plt.xlim([550, 600])
plt.ylim([0, .6])
plt.gca().tick_params(axis='y', pad=15, size=10, width=2)
plt.gca().tick_params(axis='x', pad=25, size=10, width=2)
#plt.xticks([550, 575, 600])
plt.yticks([0, .3, .6])

fig.set_tight_layout(True)
plt.savefig('trj_hetero.pdf')

# %% markdown
# Calculate and plot the sync errors

# In[]
error_homo = np.std(current_homo, axis=1)
error_hetero = np.std(current_hetero, axis=1)


# In[]
fig = plt.figure(3)

ax = plt.subplot(111)
for axis in ['top', 'right', 'bottom', 'left']:
    ax.spines[axis].set_linewidth(2)

fig.set_size_inches(23, 7)
plt.xticks(fontsize=45)
plt.yticks(fontsize=45)
plt.ylabel(r'$e$', fontsize=50, **sfont)
plt.xlabel('Time (s)', fontsize=50, **sfont)

plt.plot(t, error_homo[:120001], ls='-', lw=4, label='homogeneous')
plt.plot(t, error_hetero[:120001], ls='-', lw=4, label='heterogeneous')

plt.legend(loc='upper right', ncol=2, prop={'size': 40}, frameon=False)

plt.xlim([550, 600])
plt.ylim([0, .2])
plt.gca().tick_params(axis='y', pad=15, size=10, width=2)
plt.gca().tick_params(axis='x', pad=25, size=10, width=2)
#plt.xticks([0, 300, 600])
plt.yticks([0, .1, .2])

fig.set_tight_layout(True)
plt.savefig('error_t.pdf')

# %% markdown
# Plot the mean currents for the homogeneous and the heterogeneous systems

# In[]
fig = plt.figure(4)

ax = plt.subplot(111)
for axis in ['top', 'right', 'bottom', 'left']:
    ax.spines[axis].set_linewidth(2)

fig.set_size_inches(23, 7)
plt.xticks(fontsize=45)
plt.yticks(fontsize=45)
plt.ylabel('Mean Current (mA)', fontsize=40, **sfont)
plt.xlabel('Time (s)', fontsize=50, **sfont)

plt.plot(t, np.mean(current_homo[:120001, :],
                    axis=1), ls='-', lw=4, label='homogeneous')
plt.plot(t, np.mean(current_hetero[:120001, :],
                    axis=1), ls='-', lw=4, label='heterogeneous')

plt.legend(loc='upper right', ncol=2, prop={'size': 40}, frameon=False)

plt.xlim([550, 600])
plt.ylim([0, .6])
plt.gca().tick_params(axis='y', pad=15, size=10, width=2)
plt.gca().tick_params(axis='x', pad=25, size=10, width=2)
#plt.xticks([0, 300, 600])
plt.yticks([0, .3, .6])

fig.set_tight_layout(True)
plt.savefig('mean_current.pdf')
