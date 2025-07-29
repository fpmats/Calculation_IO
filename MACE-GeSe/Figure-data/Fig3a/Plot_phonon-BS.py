#Name: Tina Mihm
#Date: Nov 21, 2024
#Description: reads in a yaml file from phonopy

import yaml
import json
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#----------------------------------------------------------------------------------
from pylab import *
from matplotlib.colors import *
from matplotlib.font_manager import fontManager, FontProperties
from scipy.optimize import curve_fit
from matplotlib.ticker import MaxNLocator

#------------------------------------------------------------------------------
from ase.phonons import Phonons
from ase.dft.kpoints import ibz_points, bandpath
from ase.phonons import Phonons
from ase.spectrum.band_structure import BandStructure
from ase.io import read
import numpy as np
from ase.io.jsonio import read_json


###
### SET UP FIGURE
###
# rcParams['text.usetex'] = True
###
### Fonts
###
rcParams['axes.labelsize'] = 20
rcParams['xtick.labelsize'] = 20
rcParams['ytick.labelsize'] = 20
rcParams['legend.fontsize'] = 20
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = 'Computer Modern Roman'
rcParams['legend.numpoints'] = 1
rcParams['lines.linewidth'] = 3.0
rcParams['lines.markersize'] = 4.0
# http://stackoverflow.com/questions/7906365/matplotlib-savefig-plots-different-from-show
rcParams['savefig.dpi'] = 400
rcParams['figure.dpi'] = 600
###
### Size of the figure
###
# ratio=(np.sqrt(5)-1)/2.0    # golden ratio
# ratio=1                     # square figure
# plt.rcParams["figure.figsize"] = 3.37, 3.37*ratio
#rcParams[‘figure.figsize’] = 3.37, 3.37
fig = figure()

#---------------------------------------------------------
## High Sym points for VASP data  
#---------------------------------------------------------
High_sym = [[0.0, 0.0, 0.0,], [0.5, 0.0, 0.0], [0.5, 0.5, 0.0], [0.0, 0.5, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.5],  [0.0, 0.5, 0.5],[0.0, 0.0, 0.5]]
High_sym_label = ["$\Gamma$", "X", "S", "Y", "$\Gamma$", "Z", "U", "R", "T", "Z"]
print(High_sym)

#---------------------------------------------------------
## Plot colors and labels 
#---------------------------------------------------------

c_mace_vdw = "#00cc78"
c_mace_novdw = "#00ff96"

c_vasp_vdw = "#6900b2"
c_vasp_novdw = "#9600ff"

mace_color = c_mace_vdw
mace_label = "MACE[PBE]"
vasp_label = "PBE"
vasp_color = "black"

#---------------------------------------------------------
## Read in MACE 
#---------------------------------------------------------

## Reads in band structure information 
path_M = read_json('ASE-phonon-MACE/mybandpath.json')
print(path_M)
bs_M = read_json('ASE-phonon-MACE/mybandstructure.json')

for i in range(0, len(bs_M.path.kpts)): 
    kpnt = bs_M.path.kpts[i]
    Eig = bs_M.energies[0][i]
    for eig in Eig: 
        if eig < 0.0: 
            print("Frequency:", eig, "qpoint:", kpnt)

## Reading in DOS information for MACE
dos_df_M = pd.read_csv('ASE-phonon-MACE/dos_data.csv')
dos_e_M = dos_df_M["energy"]
dos_e_M = dos_e_M/0.0041 ## convert to THz
dos_w_M = dos_df_M["weight"]
dos_w_M_max = max(dos_w_M[:268])
dos_w_M = dos_w_M/dos_w_M_max

## Reading in DOS information for VASP
bs_pdos_vasp = pd.read_excel('projected_dos-VASP-noVDW.xlsx', 'VASP')
V_Freq = bs_pdos_vasp["Freq (THz)"]
V_Freq_ev = V_Freq * 0.0041 ## convert to eV
V_PDOS = bs_pdos_vasp["Pdos av"]

import matplotlib.pyplot as plt  # noqa

## Plot the phonon band structure and the DOS as two subplots: 

plt.figure(1)
fig, ax = plt.subplots()
fig = plt.figure(figsize=(7, 4))
ax = fig.add_axes([.12, .07, .67, .85])


## range for y-axis: 
emax = 7.0 
emin = -1.0

print("plotting MACE data:")
bs_M.plot(ax=ax, color = mace_color, emin=emin, emax=emax, ylabel="Energy [THz]", colorbar = False)

print("plotting MACE DOS")
dosax = fig.add_axes([.8, .07, .17, .85])
dosax.fill_between(dos_w_M, dos_e_M, y2=0, color=mace_color,
                   edgecolor="#00663C", lw=1)


dosax.set_ylim(emin, emax)
dosax.set_yticks([])
dosax.set_xticks([])
dosax.set_xlabel("DOS", fontsize=18)

print("plotting VASP DOS")
dosax.plot(V_PDOS, V_Freq, color = "black", linewidth = 1.0)

fig.savefig("Phonon-bandstructure_MACE-NoVDW_vs_PBE.png", bbox_inches='tight')


#---------------------------------------------------------
## Read in VASP data and plots the band structure
#---------------------------------------------------------

with open('band-VASP-noVDW.yaml', 'r') as file:
    data = yaml.safe_load(file)
    data2 = data["phonon"]

n = 0
frq = []
Xpnt = []
x_label = []
label = []

for i in data2: 
    Qpnt = i["q-position"]
    Dis = i["distance"]
    Dis = Dis*6.4
    x, y, z = Qpnt[0], Qpnt[1], Qpnt[2]
    if Qpnt in High_sym:
        ind = High_sym.index(Qpnt)
        label += [High_sym_label[ind]]
        x_label +=[Dis]
    print(x, y, z)
    n = n+1
    print(n)
    print(Qpnt)
    bands = i["band"]
    for l in bands: 
        t3 = l["frequency"]
        t3 = t3 #* 0.0041 ## convert to eV
        frq +=[t3]
        Xpnt +=[Dis]
        print(t3)
    plt.figure(1)
    
    ax.plot(Xpnt, frq, linestyle = "", marker = "_", color = vasp_color, markersize = 4.0)

    plt.ylabel("Frequencies [eV]")

    fig.savefig("Phonon-bandstructure_MACE-NoVDW_vs_PBE.png", bbox_inches='tight')

    frq = []
    Xpnt = []

handles, labels = ax.get_legend_handles_labels()
font = FontProperties(weight='bold')

plt.axhline(0, color = "#999999", linestyle='--', linewidth = 1.0)
fig.savefig("Phonon-bandstructure_MACE-NoVDW_vs_PBE.png", bbox_inches='tight')