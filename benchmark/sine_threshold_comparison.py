# load modules
from senn import AxonModel, StimulusModel
from senn import write_ode, evaluate_senn, analyze_solution, find_sub_supra, find_threshold
import numpy as np
import matplotlib.pylab as plt
from rich.table import Table
from rich.console import Console

# Define axon model
D = 20e-6             # fiber diamater
d = 0.7*D             # axon diamater
L = 100 * D           # internodal length
rhoi = 110*1e-2       # axoplasm resistivity
rhoe = 300 * 1e-2     # external medium resistivity
gm = 30.365*10        # membrane conductance per unit length
l = 2.5e-6            # nodal gap width
cm = 2e-2             # membrane conductance per unit area
Vr = -70e-3           # rest potential
node_num = 51         # total number of node
inl1 = 0              # first non-linear node
inl2 = 50             # last non-linear node

axon = AxonModel(D=D, d=d, L=L, rhoi=rhoi, rhoe=rhoe, gm=gm, l=l, cm=cm, Vr=Vr, node_num=node_num, inl1=inl1, inl2=inl2)

# write & import equation
T_senn = 295.16            # value used in Fortran SENN
F_senn = 96487             # value used in Fortran SENN
Nai_senn = 13.74           # value used in Fortran SENN
Vl_senn = 0.0260430075e-3  # value used in Fortran SENN

write_ode(axon.node_num, axon.inl1, axon.inl2, T=T_senn, F=F_senn, Nai=Nai_senn, Vl=Vl_senn)
from eqdiff import eqdiff

# define stimulus
tp_values = np.arange(20, 120, 20) * 1e-6
It_values = []
for tp in tp_values:

    # define stimulus
    magnitude = 0
    freq = 1 / tp
    tend = 25 * tp
    stimulus = StimulusModel(axon, 'sine', magnitude, tp, tend, freq=freq)
    ktime = 1e6
    umes_time = "($\mu$s)"

    # Identification of boundaries for bisection method
    Isub, Isup = find_sub_supra(axon, stimulus, eqdiff, sub_value=0.1e-3, sup_value=0.8e-3)
    
    # Identification of the thresholds
    It, stimulus = find_threshold(axon, stimulus, eqdiff, Isub, Isup, toll=0.5)

    It_values.append(It)


# data obtained with the original Fortran SENN model
refence = [9.125e-3, 3.9062e-3, 2.5938e-3, 1.9531e-3, 1.6094e-3]

# comparison
print("\n\n")
table = Table(title=f"Python SENN vs Fortran (original) SENN", show_lines=True)
table.add_column("Frequency (kHz)", justify="left")
table.add_column("Duration (us)", justify="left")
table.add_column("Python SENN (mA)", justify="left")
table.add_column("Fortran SENN (mA)", justify="left")
table.add_column("deviation (%)", justify="left", style='cyan')
for tp, ele, ref in zip(tp_values, It_values, refence): 
    delta = abs(ele - ref)/ref * 100
    table.add_row("{:.3f}".format(1/tp*1e-3), "{:.0f}".format(tp*1e6), "{:.6f}".format(ele*1e3), "{:.6f}".format(ref*1e3),"{:.2f}".format(delta))

console = Console()
console.print(table)
