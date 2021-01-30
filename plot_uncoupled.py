# In[]
import os
import seaborn as sns
from scipy.signal import hilbert, find_peaks
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
matplotlib.use('PDF')

#os.chdir('/home/yuanzhao/random_hetero/experiments')
#print(os.getcwd())

# In[]
t = np.linspace(0, 100, num=20001)

current_homo = np.loadtxt('oc091319_6.txt')
current_hetero = np.loadtxt('oc091319_23.txt')
current_homo = np.transpose(current_homo)
current_hetero = np.transpose(current_hetero)

# In[]
sfont = {'fontname': 'serif'}

fig = plt.figure(1)

# Remove the plot frame lines. They are unnecessary chartjunk.
ax = plt.subplot(111)
for axis in ['top', 'right', 'bottom', 'left']:
    ax.spines[axis].set_linewidth(2)
# for axis in ['top','right']:
#  ax.spines[axis].set_visible(False)
# ax.get_xaxis().tick_bottom()
# ax.get_yaxis().tick_left()

fig.set_size_inches(23, 7)
plt.xticks(fontsize=45)
plt.yticks(fontsize=45)
plt.ylabel('Current (mA)', fontsize=50, **sfont)
plt.xlabel('Time (s)', fontsize=50, **sfont)

for i in range(16):
    # ,color=(.54,.17,.89))color=tableau20[3])
    plt.plot(t, current_homo[i, :20001], ls='-', lw=2)

#plt.legend(loc='upper left',ncol=1,prop={'size':55},frameon=False)

plt.xlim([50, 100])
plt.ylim([0, .6])
plt.gca().tick_params(axis='y', pad=15, size=10, width=2)
plt.gca().tick_params(axis='x', pad=25, size=10, width=2)
#plt.xticks([550, 575, 600])
plt.yticks([0, .3, .6])

fig.set_tight_layout(True)
plt.savefig('trj_homo.pdf')


# In[]
fig = plt.figure(2)

# Remove the plot frame lines. They are unnecessary chartjunk.
ax = plt.subplot(111)
for axis in ['top', 'right', 'bottom', 'left']:
    ax.spines[axis].set_linewidth(2)
# for axis in ['top','right']:
#  ax.spines[axis].set_visible(False)
# ax.get_xaxis().tick_bottom()
# ax.get_yaxis().tick_left()

fig.set_size_inches(23, 7)
plt.xticks(fontsize=45)
plt.yticks(fontsize=45)
plt.ylabel('Current (mA)', fontsize=50, **sfont)
plt.xlabel('Time (s)', fontsize=50, **sfont)

for i in range(16):
    # ,color=(.54,.17,.89))color=tableau20[3])
    plt.plot(t, current_hetero[i, :20001], ls='-', lw=2)

#plt.legend(loc='upper left',ncol=1,prop={'size':55},frameon=False)

plt.xlim([50, 100])
plt.ylim([0, .6])
plt.gca().tick_params(axis='y', pad=15, size=10, width=2)
plt.gca().tick_params(axis='x', pad=25, size=10, width=2)
#plt.xticks([550, 575, 600])
plt.yticks([0, .3, .6])

fig.set_tight_layout(True)
plt.savefig('trj_hetero.pdf')

# In[]
shifts = np.reshape((np.amax(current_homo, axis=1) +
                     np.amin(current_homo, axis=1)) / 2, (16, 1))
z = hilbert(current_homo - shifts)  # form the analytical signal
inst_amplitude = np.abs(z)  # envelope extraction
inst_phase = np.unwrap(np.angle(z))  # inst phase
# inst_freq = np.diff(inst_phase)/(2*np.pi)*fs 	#inst frequency
r = np.abs(np.mean(np.exp(1j * inst_phase), axis=0))


# In[]
fig = plt.figure(3)

# Remove the plot frame lines. They are unnecessary chartjunk.
ax = plt.subplot(111)
for axis in ['top', 'right', 'bottom', 'left']:
    ax.spines[axis].set_linewidth(2)
# for axis in ['top','right']:
#  ax.spines[axis].set_visible(False)
# ax.get_xaxis().tick_bottom()
# ax.get_yaxis().tick_left()

fig.set_size_inches(23, 7)
plt.xticks(fontsize=45)
plt.yticks(fontsize=45)
plt.ylabel(r'$\phi_i-\phi_1$', fontsize=50, **sfont)
plt.xlabel('Time (s)', fontsize=50, **sfont)

for i in range(16):
    plt.plot(t, inst_phase[i, :20001] - inst_phase[0, :20001], ls='-', lw=2)

#plt.legend(loc='upper left',ncol=1,prop={'size':55},frameon=False)

plt.xlim([0, 100])
#plt.ylim([-5, 100])
plt.gca().tick_params(axis='y', pad=15, size=10, width=2)
plt.gca().tick_params(axis='x', pad=25, size=10, width=2)
#plt.xticks([550, 575, 600])
# plt.yticks([0,.3,.6])

fig.set_tight_layout(True)
plt.savefig('phase_homo.pdf')


# In[]
fig = plt.figure(4)

# Remove the plot frame lines. They are unnecessary chartjunk.
ax = plt.subplot(111)
for axis in ['top', 'right', 'bottom', 'left']:
    ax.spines[axis].set_linewidth(2)
# for axis in ['top','right']:
#  ax.spines[axis].set_visible(False)
# ax.get_xaxis().tick_bottom()
# ax.get_yaxis().tick_left()

fig.set_size_inches(23, 7)
plt.xticks(fontsize=45)
plt.yticks(fontsize=45)
plt.ylabel(r'$A_i$', fontsize=50, **sfont)
plt.xlabel('Time (s)', fontsize=50, **sfont)

for i in range(16):
    plt.plot(t, inst_amplitude[i, :20001], ls='-', lw=2)
    #plt.plot(t, np.real(z[i, :120001]), ls='-', lw=2)
    #plt.plot(t, np.imag(z[i, :120001]), ls='-', lw=2)

#plt.legend(loc='upper left',ncol=1,prop={'size':55},frameon=False)

plt.xlim([50, 100])
# plt.ylim([0,.6])
plt.gca().tick_params(axis='y', pad=15, size=10, width=2)
plt.gca().tick_params(axis='x', pad=25, size=10, width=2)
#plt.xticks([550, 575, 600])
# plt.yticks([0,.3,.6])

fig.set_tight_layout(True)
plt.savefig('amplitude_homo.pdf')


# In[]
fig = plt.figure(5)

# Remove the plot frame lines. They are unnecessary chartjunk.
ax = plt.subplot(111)
for axis in ['top', 'right', 'bottom', 'left']:
    ax.spines[axis].set_linewidth(2)
# for axis in ['top','right']:
#  ax.spines[axis].set_visible(False)
# ax.get_xaxis().tick_bottom()
# ax.get_yaxis().tick_left()

fig.set_size_inches(23, 7)
plt.xticks(fontsize=45)
plt.yticks(fontsize=45)
plt.ylabel(r'$R$', fontsize=50, **sfont)
plt.xlabel('Time (s)', fontsize=50, **sfont)

plt.plot(t, r[:20001], ls='-', lw=2)

#plt.legend(loc='upper left',ncol=1,prop={'size':55},frameon=False)

# plt.xlim([550,600])
# plt.ylim([0,.6])
plt.gca().tick_params(axis='y', pad=15, size=10, width=2)
plt.gca().tick_params(axis='x', pad=25, size=10, width=2)
# plt.xticks([550,575,600])
# plt.yticks([0,.3,.6])

fig.set_tight_layout(True)
plt.savefig('order_homo.pdf')


# In[]
shifts = np.reshape((np.amax(current_hetero, axis=1) +
                     np.amin(current_hetero, axis=1)) / 2, (16, 1))
z = hilbert(current_hetero - shifts)  # form the analytical signal
inst_amplitude = np.abs(z)  # envelope extraction
inst_phase = np.unwrap(np.angle(z))  # inst phase
# inst_freq = np.diff(inst_phase)/(2*np.pi)*fs 	#inst frequency
r = np.abs(np.mean(np.exp(1j * inst_phase), axis=0))


# In[]
fig = plt.figure(6)

# Remove the plot frame lines. They are unnecessary chartjunk.
ax = plt.subplot(111)
for axis in ['top', 'right', 'bottom', 'left']:
    ax.spines[axis].set_linewidth(2)
# for axis in ['top','right']:
#  ax.spines[axis].set_visible(False)
# ax.get_xaxis().tick_bottom()
# ax.get_yaxis().tick_left()

fig.set_size_inches(23, 7)
plt.xticks(fontsize=45)
plt.yticks(fontsize=45)
plt.ylabel(r'$\phi_i-\phi_1$', fontsize=50, **sfont)
plt.xlabel('Time (s)', fontsize=50, **sfont)

for i in range(16):
    plt.plot(t, inst_phase[i, :20001] - inst_phase[0, :20001], ls='-', lw=2)

#plt.legend(loc='upper left',ncol=1,prop={'size':55},frameon=False)

plt.xlim([0, 100])
#plt.ylim([-5, 100])
plt.gca().tick_params(axis='y', pad=15, size=10, width=2)
plt.gca().tick_params(axis='x', pad=25, size=10, width=2)
#plt.xticks([550, 575, 600])
# plt.yticks([0,.3,.6])

fig.set_tight_layout(True)
plt.savefig('phase_hetero.pdf')


# In[]
fig = plt.figure(7)

# Remove the plot frame lines. They are unnecessary chartjunk.
ax = plt.subplot(111)
for axis in ['top', 'right', 'bottom', 'left']:
    ax.spines[axis].set_linewidth(2)
# for axis in ['top','right']:
#  ax.spines[axis].set_visible(False)
# ax.get_xaxis().tick_bottom()
# ax.get_yaxis().tick_left()

fig.set_size_inches(23, 7)
plt.xticks(fontsize=45)
plt.yticks(fontsize=45)
plt.ylabel(r'$A_i$', fontsize=50, **sfont)
plt.xlabel('Time (s)', fontsize=50, **sfont)

for i in range(16):
    plt.plot(t, inst_amplitude[i, :20001], ls='-', lw=2)

#plt.legend(loc='upper left',ncol=1,prop={'size':55},frameon=False)

plt.xlim([50, 100])
# plt.ylim([0,.6])
plt.gca().tick_params(axis='y', pad=15, size=10, width=2)
plt.gca().tick_params(axis='x', pad=25, size=10, width=2)
#plt.xticks([550, 575, 600])
# plt.yticks([0,.3,.6])

fig.set_tight_layout(True)
plt.savefig('amplitude_hetero.pdf')


# In[]
fig = plt.figure(8)

# Remove the plot frame lines. They are unnecessary chartjunk.
ax = plt.subplot(111)
for axis in ['top', 'right', 'bottom', 'left']:
    ax.spines[axis].set_linewidth(2)
# for axis in ['top','right']:
#  ax.spines[axis].set_visible(False)
# ax.get_xaxis().tick_bottom()
# ax.get_yaxis().tick_left()

fig.set_size_inches(23, 7)
plt.xticks(fontsize=45)
plt.yticks(fontsize=45)
plt.ylabel(r'$R$', fontsize=50, **sfont)
plt.xlabel('Time (s)', fontsize=50, **sfont)

plt.plot(t, r[:20001], ls='-', lw=2)

#plt.legend(loc='upper left',ncol=1,prop={'size':55},frameon=False)

# plt.xlim([550,600])
# plt.ylim([0,.6])
plt.gca().tick_params(axis='y', pad=15, size=10, width=2)
plt.gca().tick_params(axis='x', pad=25, size=10, width=2)
# plt.xticks([550,575,600])
# plt.yticks([0,.3,.6])

fig.set_tight_layout(True)
plt.savefig('order_hetero.pdf')
