#Name: Tina Mihm
#Date: May 5, 2025
#Description: Pulls in data from the ASE database and graph it as contour plots and scatter plots

import os
import numpy as np
from ase import Atoms
from ase.io import read, write
import matplotlib.pyplot as plt
import pandas as pd
from ase.geometry import wrap_positions
from ase.io.vasp import write_vasp
#-------------------------------------------------------------------------------------------
from pylab import *
from matplotlib.colors import *
from matplotlib.font_manager import fontManager, FontProperties
from scipy.optimize import curve_fit
from matplotlib.ticker import MaxNLocator

from ase.db import connect

###
### SET UP FIGURE
###
# rcParams['text.usetex'] = True
###
### Fonts
###
rcParams['axes.labelsize'] = 12
rcParams['xtick.labelsize'] = 12
rcParams['ytick.labelsize'] = 12
rcParams['legend.fontsize'] = 10
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = 'Computer Modern Roman'
rcParams['legend.numpoints'] = 1
rcParams['lines.linewidth'] = 2.0
rcParams['lines.markersize'] = 6.0
# http://stackoverflow.com/questions/7906365/matplotlib-savefig-plots-different-from-show
rcParams['savefig.dpi'] = 400
rcParams['figure.dpi'] = 600
###
### Size of the figure
###
# ratio=(np.sqrt(5)-1)/2.0    # golden ratio
# ratio=1                     # square figure
# plt.rcParams["figure.figsize"] = 3.37, 3.37*ratio
# fig = figure()

#---------------------------------------------------------
## System Data
#---------------------------------------------------------
Points = "5000-pnts-"
MACE_model = "MACErcut4.0"

dir_path = "./Database/MACE"

from matplotlib.colors import LinearSegmentedColormap


cmap = 'jet'

levels = 30

#--------------------------------------------------------------------------------------
## Pull out data from excel and analize
#--------------------------------------------------------------------------------------
cmax = 0
VMIN = 0

## Pull out the min energy from the delta d_i = 0.0 data base to use as reference energy for all plots ##
for filename in os.listdir(dir_path):
    if filename.startswith("relax-11A"):
        file_path = os.path.join(dir_path, filename)
        print(file_path)
        #------
        db = connect(file_path)
        #-----
        selections = ['model=macecalculator']
        natoms = db[1].natoms
        A_lat = db[1].a
        B_lat = db[1].b
        layers = natoms/4
        print("number of atoms:", natoms)
        for selection in selections:
            energy = np.array([row.energy for row in db.select(selection)])
            print(energy)
            vmin = np.min([row.energy for row in db.select(selection)])
            VMIN = vmin
            print("VMIN", VMIN)

#---------------------------------------------------------------------
## Plot 3D graph from 2D slices
#---------------------------------------------------------------------
Off = []
Dist = []
levels = 100

print("Figure 3")

fig2 = plt.figure(3)
ax2 = fig2.add_subplot(projection = '3d')

## set range for color bar
cmin = 0
cmax = 400 


order = [10.5, 11.5, 11] # defines order of plotting on z axis
ticks = [-0.25, 0.25, 0.0]# defines label for ticks on z axis


for d in order: 
    offset = d
    if offset == 10.5: 
        offset = d
    if offset == 11: 
        offset = d+10
    if offset == 11.5: 
        offset = d+20
    Off +=[offset]
    for filename in os.listdir(dir_path):
        if filename.startswith("relax-"+str(d)):
            Fname = re.findall(r"\d+\.?\d*", filename)
            System = Fname[0]+"A-"
            L_dis = float(Fname[0])
            Dist +=[L_dis]
            file_path = os.path.join(dir_path, filename)
            #------
            db = connect(file_path)
            #-----
            selections = ['model=macecalculator']
                    
            for selection in selections:
                m = selection.split("=")
                model = m[1]+"-"
                print(model)
                dB_array = np.array([row.dB for row in db.select(selection)])
                dA_array = np.array([row.dA for row in db.select(selection)])
                energy = np.array([row.energy for row in db.select(selection)])
                print(energy)
                Ediff = (energy - VMIN)

                print(VMIN)
                print("Min energy @ delta d_i = 0.0", VMIN)

                print(Ediff)
                print("Ediff min and max (meV/layer):", min(Ediff)*1000/layers, max(Ediff)*1000/layers)
                print(cmin, cmax)
                print(cmax)
                print(levels)
                contour2 = ax2.tricontourf(dA_array, dB_array, Ediff*1000/layers, zdir='z', offset=offset,levels=levels, cmap=cmap, vmax = cmax)


# Set axes, color bar, legend and plot name, etc
##---------------------------------------------
print("This is Cmin/max", cmin, cmax)

cbar2 = plt.cm.ScalarMappable(cmap=cmap)

cbar2.set_array([0, cmax])
plt.colorbar(cbar2, shrink=0.5, aspect=10, extend='max', label='$\Delta E$ (meV/layer)', pad=0.0)

            
            #------------------------
ax2.set_xlabel(r'$\Delta a$ ($\AA$)', labelpad=0.5)
ax2.set_ylabel(r'$\Delta b$ ($\AA$)', labelpad=0.5)
ax2.zaxis.set_rotate_label(False)
# ax2.set_zlabel(r'$\Delta c \AA$', rotation = 90, labelpad=0.5)
ax2.set_zlabel(r'$\Delta d_i$ ($\AA$)', rotation = 90, labelpad=1.0)
ax2.set_xticks(np.arange(0, A_lat, 1))
ax2.set_yticks(np.arange(0, B_lat, 1))
ax2.set_zlim(9, 13)
ax2.set_zticks(Off, labels=[l for l in ticks], rotation = 30)

# Change axes tick mark labels
ax2.tick_params(axis='x', pad=-1)
ax2.tick_params(axis='y', pad=-1)
ax2.tick_params(axis='z', pad=1)

# Change axes colors
ax2.xaxis.line.set_color('gray')
ax2.yaxis.line.set_color('gray')
ax2.zaxis.line.set_color('white')

# Set background color
ax2.set_facecolor((0,0,0,0))
ax2.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax2.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax2.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
# make the grid lines transparent
ax2.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
ax2.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
ax2.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)


# sets angle view along Z (elev), and angle along x/y (azim)
ax2.view_init(elev=15, azim=60)

fig2.tight_layout(pad=-0.5)

fig2.savefig("MACE/GeSe-Grid-Search-"+Points+MACE_model+"-3D-Rel-11A-min-energy.png", transparent=True)

##------------------------------------------
