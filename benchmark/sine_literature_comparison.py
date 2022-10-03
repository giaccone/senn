# load modules
from senn import AxonModel, StimulusModel
from senn import write_ode, evaluate_senn, analyze_solution
import numpy as np
import matplotlib.pylab as plt

# Define axon model
D = 20e-6           # fiber diamater
d = 0.7*D           # axon diamater
L = 100 * D         # internodal length
rhoi = 110*1e-2     # axoplasm resistivity
rhoe = 300 * 1e-2   # external medium resistivity
gm = 30.4*10        # membrane conductance per unit length
l = 2.5e-6          # nodal gap width
cm = 2e-2           # membrane conductance per unit area
Vr = -70e-3         # rest potential
node_num = 51       # total number of node
inl1 = 26 - 1       # first non-linear node
inl2 = 26 + 5       # last non-linear node

axon = AxonModel(D=D, rhoi=rhoi, rhoe=rhoe, gm=gm, l=l, cm=cm, Vr=Vr, node_num=node_num, inl1=inl1, inl2=inl2)

# write & import equation
write_ode(axon.node_num, axon.inl1, axon.inl2)
from eqdiff import eqdiff

# define stimulus
magnitude = -9.105e-3
tp = 20e-6
freq = 1 / tp
tend = 25 * tp
stimulus = StimulusModel(axon, 'sine', magnitude, tp, tend, freq=freq)
ktime = 1e6
umes_time = "($\mu$s)"

# one shot simulation
t, sol = evaluate_senn(axon, stimulus, eqdiff)

# analyze output
analyze_solution(t, sol, axon, stimulus)

# plot stimulus
plt.figure(facecolor='w')
plt.plot(t * ktime, -stimulus.waveform(t) * 1e3, linewidth=2)
plt.xlabel('time ' + umes_time, fontsize=16)
plt.ylabel('stimulus (mA)', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.tight_layout()

# plot stimulus
plt.figure(facecolor='w')
plt.plot(t * ktime, -stimulus.voltage_ext(t, stimulus.inod[0]) * 1e3, linewidth=2)
plt.xlabel('time ' + umes_time,fontsize=16)
plt.ylabel('external potential (mV)', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.tight_layout()

# plot nodes
plt.figure(facecolor='w')
plt.plot(t * ktime, sol[:, stimulus.inod[:3]] * 1e3, linewidth=2)
plt.xlabel('time ' + umes_time, fontsize=16)
plt.ylabel('external potential (mV)', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.tight_layout()

# load reference data
node0 = np.loadtxt('./reference_data/Reilly1985_sine_node0.txt')
plt.plot(node0[:,0]*1e3, node0[:,1], 'C0o', mfc='None')
lg = plt.legend(stimulus.leg[:3] + ['Reilly 1985'],loc='best',bbox_to_anchor=(0.4, 0.55))
lg.set_draggable(True)
plt.xlim(0, 500)
plt.tight_layout()
plt.show()