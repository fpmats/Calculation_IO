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
rcParams['axes.labelsize'] = 20
rcParams['xtick.labelsize'] = 18
rcParams['ytick.labelsize'] = 20
rcParams['legend.fontsize'] = 25
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
dir_path2 = "./Database/VASP"

levels = 30

#--------------------------------------------------------------------------------------
## Pull out data from excel and analize
#--------------------------------------------------------------------------------------
cmax = 0
VMIN = 0

## Pull out the Vmin from the delta d_i = 0.0 database for MACE to use as reference energy for all MACE plots ##
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
            vmin = np.min([row.energy for row in db.select(selection)])
            VMIN = vmin

VMIN_V = 0

## Pull out the Vmin from the delta d_i = 0.0 database for PBE to use as reference energy for all PBE plots ##
for filename in os.listdir(dir_path2):
    if filename.startswith("relax-11A"):
        file_path = os.path.join(dir_path2, filename)
        print(file_path)
        #------
        db = connect(file_path)
        #-----
        selections = ['model=vasp']
        natoms = db[1].natoms
        A_lat = db[1].a
        B_lat = db[1].b
        layers = natoms/4
        print("number of atoms:", natoms)
        for selection in selections:
            energy = np.array([row.energy for row in db.select(selection)])
            vmin = np.min([row.energy for row in db.select(selection)])
            VMIN_V = vmin


#---------------------------------------------------------------------
## Plot energy vs dB for one dA for all layers
#---------------------------------------------------------------------

print("Figure 2")

plt.figure(2)
fig, ax = plt.subplots()

DeltaC = []
E_dC = []

t = 0

## Loops over databases, pulls data and graphs as a scatter plot 
##----------------------------------------------------------------

## MACE data: 

for filename in os.listdir(dir_path):
    if filename.startswith("relax"):
        print(filename)
        Fname = re.findall(r"\d+\.?\d*", filename)
        # print(Fname)
        System = Fname[0]+"A-"
        L_dis = float(Fname[0])
        
        file_path = os.path.join(dir_path, filename)
        print("This is data for:", file_path)
        #------
        db = connect(file_path)
        #-----
        selections = ['model=macecalculator']

        for selection in selections:
            m = selection.split("=")
            model = m[1]+"-"

            dB_array = np.array([row.dB for row in db.select(selection)])
            dA_array = np.array([row.dA for row in db.select(selection)])
            dB_sort = np.sort(np.unique(dB_array))
            dA_sort = np.sort(np.unique(dA_array))

            m = selection.split("=")
            model = m[1]+"-"

            ## finds selected dA value for scatter plot and pulls all dB energies and dB values for that point: 
            #------------------------------------------------------------------------------
            print("all dB data pull for dA = ", dA_sort[-7]) ## dA = 3.5226132040826084
            dA_val = dA_sort[-7]
            db_val_list = [dB_sort[1], dB_sort[18], dB_sort[35], dB_sort[51], dB_sort[68]]
            dB_array = np.array([row.dB for row in db.select(f'dA={dA_val},'+selection)])
            E_array = np.array([row.energy for row in db.select(f'dA={dA_val},'+selection)])
            dB_val = min(dB_sort)
            E_dA_dB0 = [row.energy for row in db.select(f'dA={dA_val}, '+f'dB={dB_val}, '+selection)]
            vmin = np.min(E_array)
            vmax = np.max(E_array)
            print("METHOD:", selection, "C LATTIC:",L_dis, "delta a:", dA_val, "MIN E:",  vmin, "MAX E:", vmax )
            print("METHOD:", selection, "C LATTIC:",L_dis, "delta a:", dA_val, "delta b:", dB_val,  "E@dB:",  E_dA_dB0 )
            #--------------------------------------------------------------------------------------------------------------
            

            E_val_list = []
            for i in db_val_list: 
                print(i)
                E_val_list += [db.get(f'dA={dA_val}, dB={i},'+selection).energy]
                if round(i, 8) == 2.29608881: 
                    print(i)
                    DeltaC +=[L_dis]
                    E_dC += [db.get(f'dA={dA_val}, dB={i},'+selection).energy]

            Ediff = (E_array - VMIN)
            print("Difference between min energy for d_i = 11 and min energy for selected database:", min(Ediff))

            Ediff2 = Ediff*1000/layers

            ## pandas dataframe is used to sort data to allow for a line plot to be used 
            #-----------------------------------------------------------------------------
            df = {"x":dB_array, "y": Ediff2}
            df = pd.DataFrame(df)
            df_sort = df.sort_values(by='x')
            #---------------------------------
            print("plotting MACE data")
            if L_dis == 11.0:
                ax.plot("x", "y", data = df_sort, color='#00b269', marker = "o", linestyle='-', label = "MACE[PBE]+D3, $\Delta d_i = 0.0 \AA$",zorder = 20)
            if L_dis == 11.5:
                ax.plot("x", "y", data = df_sort, color='#00663c', marker = "o", linestyle='-', label = "MACE[PBE]+D3, $\Delta d_i = 0.25 \AA$",zorder = 25)
            if L_dis == 10.5:
                ax.plot("x", "y", data = df_sort, color='#00ff96', marker = "o", linestyle='-', label = "MACE[PBE]+D3, $\Delta d_i = -0.25 \AA$", zorder = 15)

            # #------------------------
 

for filename in os.listdir(dir_path2):
    if filename.startswith("relax"):
        Fname = re.findall(r"\d+\.?\d*", filename)
        System = Fname[0]+"A-"
        L_dis = float(Fname[0])
        
        file_path = os.path.join(dir_path2, filename)
        print("This is data for:", file_path)
        #------
        db = connect(file_path)
        #-----
        selections = ['model=vasp']

        for selection in selections:
            m = selection.split("=")
            model = m[1]+"-"
            dB_array = np.array([row.dB for row in db.select(selection)])
            dA_array = np.array([row.dA for row in db.select(selection)])
            dB_sort = np.sort(np.unique(dB_array))
            dA_sort = np.sort(np.unique(dA_array))

            m = selection.split("=")
            model = m[1]+"-"

            ## finds selected dA value for scatter plot and pulls all dB energies and dB values for that point: 
            #------------------------------------------------------------------------------
            print("all dB data pull for dA = ", dA_sort[-7]) ## dA = 3.5226132040826084
            dA_val = dA_sort[-7]
            db_val_list = [dB_sort[1], dB_sort[18], dB_sort[35], dB_sort[51], dB_sort[68]]
            dB_array = np.array([row.dB for row in db.select(f'dA={dA_val},'+selection)])
            E_array = np.array([row.energy for row in db.select(f'dA={dA_val},'+selection)])
            dB_val = min(dB_sort)
            E_dA_dB0 = [row.energy for row in db.select(f'dA={dA_val}, '+f'dB={dB_val}, '+selection)]
            vmin = np.min(E_array)
            vmax = np.max(E_array)
            print("METHOD:", selection, "C LATTIC:",L_dis, "delta a:", dA_val, "MIN E:",  vmin, "MAX E:", vmax )
            print("METHOD:", selection, "C LATTIC:",L_dis, "delta a:", dA_val, "delta b:", dB_val,  "E@dB:",  E_dA_dB0 )
            #--------------------------------------------------------------------------------------------------------------

            Ediff = (E_array - VMIN_V)
            print("Difference between min energy for d_i = 11 and min energy for selected database:", min(Ediff))

            print8=("Plotting PBE data")
            if L_dis == 11.0: 
                ax.plot(dB_array, Ediff*1000/layers, color='#666666', marker = "o", linestyle='-', label = "PBE+D3, $\Delta d_i = 0.0 \AA$",zorder = 5)
            if L_dis == 11.5: 
                ax.plot(dB_array, Ediff*1000/layers, color='black', marker = "o", linestyle='-', label = "PBE+D3, $\Delta d_i = 0.25 \AA$",zorder = 10)
            if L_dis == 10.5: 
                ax.plot(dB_array, Ediff*1000/layers, color='#b2b2b2', marker = "o", linestyle='-', label = "PBE+D3, $\Delta d_i = -0.25 \AA$", zorder = 1)

            # #------------------------


## setting up axies, legend, plot name, etc
#-------------------------------------------
ax.set_xlim(-0.1, max(dB_array)+0.1 )
ax.set_xticks(np.arange(0, B_lat, 0.5))
ax.set_ylim(-0.1, 750)
ax.set_xlabel(r'$\Delta b$')
ax.set_ylabel(r'$\Delta E$ (meV/layer)')

handles,labels = ax.get_legend_handles_labels()

handles = [handles[2], handles[0], handles[1], handles[5], handles[3], handles[4]]
labels = [labels[2], labels[0], labels[1], labels[5], labels[3], labels[4]]

leg = plt.legend(handles,labels, loc='center', bbox_to_anchor=(1.32, 0.5), ncol=1, handlelength = 2.0, handletextpad = 0.5, columnspacing=1.0, fontsize = 12, borderpad=0.1, frameon=False)

fig.savefig("Both/GeSe-Grid-Search-energy-vs-dB_at_dA-"+MACE_model+"-all-Rel-11A-min-energy.png", bbox_inches='tight')
plt.close()
plt.clf() 


