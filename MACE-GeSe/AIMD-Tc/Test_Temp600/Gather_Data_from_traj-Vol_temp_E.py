#Name: Tina Mihm
#Date: 04/17/2025
#Description: Pulls the volume, temp, timesteps,energy and a/b lattice parameters from the traj files and graphs them vs time

import csv

import numpy as np
import subprocess
import re
import pandas as pd
import math as math
from ase.geometry.analysis import Analysis
##import seaborn as sns
import matplotlib.pyplot as plt
from pylab import *
from matplotlib.colors import *
from matplotlib.font_manager import fontManager, FontProperties
from scipy.optimize import curve_fit
from scipy.optimize import least_squares
from matplotlib.ticker import MaxNLocator
import os

#---------------------------
#rcParams['savefig.dpi'] = 100
#rcParams['figure.dpi'] = 200
###
### Size of the figure
###
ratio=(np.sqrt(5)-1)/2.0     # golden ratio
#ratio=1                     # square figure
plt.rcParams["figure.figsize"] = 3.37, (3.37)*ratio
matplotlib.rcParams.update({'errorbar.capsize': 2.5})
#rcParams[‘figure.figsize’] = 3.37, 3.37
fig = figure()
#--------------------------
from ase import Atoms
from ase.io.trajectory import Trajectory
from ase.io.vasp import write_vasp
from ase.io.xyz import write_xyz
from ase.io import read
from ase.io import write
from os import system


#########################################
##  System info  ##
########################################

traj = Trajectory('traj/md.traj')

System = "Original_MD_600K_1000ps"
psmax = 1000
Tmax =600
start_ps = 0

#########################################
##  pull data from trajectory  ##
########################################

print(len(traj))
atoms = traj[0]
angle = atoms.get_cell_lengths_and_angles()
print(angle)

T = []
V = []
A = []
B = []

for i in range(0, len(traj)):
    atoms = traj[i]
    av = atoms.get_volume()
    at = atoms.get_temperature()
    V +=[av]
    Cell = atoms.get_cell()
    ## get lattice parameters to scale coordinates
    A += [Cell[0][0]]
    B += [Cell[1][1]]
    T +=[at]

AB_ratio = [x / y for x, y in zip(A, B)]

cwd = os.getcwd()
dir_path = cwd + r"/log"

Time = []
Etot = []
Epot = []
Ekin = []
Temp = []

for fl in os.listdir(dir_path):
    file = os.path.join(dir_path, fl)
    with open(file, 'r') as fp:
        for l_no, line in enumerate(fp):
            if l_no == 0:
                continue
            #print(l_no)
            #print(line)
            str = line.split()
            #print(str)
            Time += [float(str[0])]
            Etot += [float(str[1])]
            Epot += [float(str[2])]
            Ekin += [float(str[3])]
            Temp += [float(str[4])]
#########################################
##  Save data to a .csv  ##
########################################

start = len(traj) - start_ps

print(len(AB_ratio))
print(len(Time))
print(len(V))
print(len(T))
print(len(Ekin))
print(len(Epot))
print(len(Etot))

df2 = {"Time[ps]": Time, "T (K)": T, "V (A^3)":V,"a/b ratio":AB_ratio, "E_tot (eV)": Etot, "E_pot (eV)": Epot, "E_k (eV)": Ekin }
df2 = pd.DataFrame(df2)
df2.to_csv(System+'-energy-volume-temp-data.csv', index=False)
print(df2)


#########################################
##  Graph data  ##
########################################

### Graph of V vs Time ###

plt.figure(1)
plt.plot(Time, V, linestyle = "-", marker = "o", markersize =2,  color = "#00ff00")
#plt.xlim(0, psmax)
##plt.ylim(-1.5, 2.5)
plt.ylabel(r"Volume [$\AA^3$]")
plt.xlabel(r"Time [ps]")
plt.savefig(System+"-Volumes_vs_Time.png", bbox_inches='tight')
#plt.show()

#---------------------------

plt.figure(2)
plt.plot(Time, Temp, linestyle = "-", marker = "o", markersize =2,  color = "#00ff00")
#plt.xlim(0, psmax)
plt.ylim(Tmax - 50, Tmax + 50)
plt.ylabel(r"Temp [K]")
plt.xlabel(r"Time [ps]")
plt.savefig(System+"-Temp_vs_Time.png", bbox_inches='tight')

#-----------------------------

plt.figure(3)
plt.plot(Time, Etot, linestyle = "-", marker = "o", markersize =2,  color = "#00ff00")
#plt.xlim(0, psmax)
##plt.ylim(-1.5, 2.5)
plt.ylabel(r"E_tot [eV]")
plt.xlabel(r"Time [ps]")
plt.savefig(System+"-TotE_vs_Time.png", bbox_inches='tight')

#-----------------------------

plt.figure(4)
plt.plot("Time[ps]", "a/b ratio",data = df2, linestyle = "-", marker = "o", markersize =2,  color = "#00ff00")
#plt.xlim(0, psmax)
##plt.ylim(-1.5, 2.5)
plt.ylabel(r"a/b lattice ratio")
plt.xlabel(r"Time [ps]")
plt.savefig(System+"-AB-ratio_vs_Time.png", bbox_inches='tight')

