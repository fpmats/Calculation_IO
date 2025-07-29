#Name: Tina Mihm
#Date: May 1, 2025
#Description: takes the Theta gathered across MD trajectories for each temp and graphs them as a density plot 

import numpy as np
import subprocess
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pylab import *
from matplotlib.colors import *
from matplotlib.font_manager import fontManager, FontProperties
from scipy.optimize import curve_fit
from matplotlib.ticker import MaxNLocator

from ase.io import read, write
from ase.io.jsonio import read_json
import matplotlib.pyplot as plt

import os

#-------------------------------------------------------------------------------------------
###
### SET UP FIGURE
###
# rcParams['text.usetex'] = True
###
### Fonts
###
rcParams['axes.labelsize'] = 8
rcParams['xtick.labelsize'] = 8
rcParams['ytick.labelsize'] = 8
rcParams['legend.fontsize'] = 8
rcParams['font.family'] = 'serif'
# rcParams['font.serif'] = 'Computer Modern Roman'
rcParams['legend.numpoints'] = 1
rcParams['lines.linewidth'] = 1.0
rcParams['lines.markersize'] = 2.0
# http://stackoverflow.com/questions/7906365/matplotlib-savefig-plots-different-from-show
rcParams['savefig.dpi'] = 500
rcParams['figure.dpi'] = 600
###
### Size of the figure
# ###
ratio=(np.sqrt(5)-1)/2.0    # golden ratio
# ratio=1                     # square figure
plt.rcParams["figure.figsize"] = 3.37, 3.37*ratio
# #rcParams[‘figure.figsize’] = 3.37, 3.37
fig = figure()

#-------------------------------------------------------------------------------------
from scipy.stats import gaussian_kde
import mpl_scatter_density
from mpl_scatter_density import ScatterDensityArtist
from matplotlib.colors import LinearSegmentedColormap

#--------------------------------------------------------------------------------------
## System information
#--------------------------------------------------------------------------------------

dir_path = "./All_Theta"
filename = "/Original_MD_200K_1000ps-All_Thata.csv"
file_path = dir_path+filename

#--------------------------------------------------------------------------------------
## Read in data from excel 
#--------------------------------------------------------------------------------------

df = pd.read_excel('MACE-PBE_D3-cutoff15-MD.xlsx') ## averaged data

df_all = pd.read_csv(file_path)

#--------------------------------------------------------------------------------------
## Pull out data from excel and analize
#--------------------------------------------------------------------------------------

#Averaged data: 
Temp = df["temp (K)"]
Theta_L1 = df["Av Theta for L1 (rad)"]
Theta_L2 = df["Av Theta for L2 (rad)"]
Theta_L3 = df["Av Theta for L3 (rad)"]
Theta_L4 = df["Av Theta for L4 (rad)"]
Theta_av = df["Av Theta (rad)"]
Theta_abs_av = df["Abs Av Theta (rad)"]

#All Theta data: 

temp_all = []
Theta_all = []


for filename in os.listdir(dir_path):
    if filename.endswith(".csv"):
        Fname = re.findall(r'\d+', filename)
        print(Fname)
        tmp = int(Fname[0])
        file_path = os.path.join(dir_path, filename)
        print(file_path)
        df_all = pd.read_csv(file_path)
        ATheta = df_all["All Theta"]
        ATheta_abs = df_all["All Theta, abs"]
        print(np.mean(ATheta_abs))
        Theta_all += [ATheta]
        temp_all += [tmp]

#--------------------------------------------------------------------------------------
## Plot data as a density map
#--------------------------------------------------------------------------------------
import matplotlib.cm as cm
import matplotlib.colors as colors

white_viridis = LinearSegmentedColormap.from_list('white_viridis', [
    (0, '#ffffff'),
    (0.2, '#404388'),
    (0.4, '#2a788e'),
    (0.6, '#21a784'),
    (0.8, '#78d151'),
    (1, '#fde624'),
], N=len(ATheta))

print(white_viridis)

colors_array = white_viridis(np.linspace(0, 1, 50))
print(f"Array of 5 colors: {colors_array}")

colors_array[0][3] = 0.0

viridis_cmap = cm.get_cmap('viridis')

min_color = viridis_cmap(0.0)

print(min_color)

colors_array_new = viridis_cmap(np.linspace(0, 1, 5))
print(f"Array of 5 colors: {colors_array_new}")

new_color = [0, 0, 0, 0]  # RGBA (white, 0% opacity)
new_colors = np.vstack((new_color, colors_array_new))
print(new_colors)

alpha_viridis = LinearSegmentedColormap.from_list('alpha_viridis', new_colors, N=len(Theta_all[-1]))

# #---------------------------------------------------------
# #Just density map 
# #---------------------------------------------------------

# for i in range(0, len(Theta_all)):
#     print(temp_all[i])
#     fig = plt.figure(1)
#     ax = fig.add_subplot(1, 1, 1, projection='scatter_density')
#     if temp_all[i] < 250: 
#         density = ax.scatter_density([temp_all[i]] * len(Theta_all[i]), Theta_all[i], cmap = white_viridis, vmax = 20000, dpi=50)
#     if temp_all[i] > 250: 
#         density = ax.scatter_density([temp_all[i]] * len(Theta_all[i]), Theta_all[i], cmap = white_viridis, vmax = 20000, dpi=50)
#     ax.set_ylim(-0.75, 0.75)
#     ax.set_xlim(0, 750)
#     ax.set_ylabel(r"Angle, $\theta$ (rad)")
#     ax.set_xlabel(r"Temp (K)")
#     fig.colorbar(density, extend='max', label = r"Number of $\theta$")
#     plt.savefig("Theta_vs_Temp_density-allTemp_Graph-"+str(temp_all[i])+".png", bbox_inches='tight')

#---------------------------------------------------------
#Density map with averages
#---------------------------------------------------------


fig = plt.figure(2)
for i in range(0, len(Theta_all)):
    print(temp_all[i])
    fig = plt.figure(1)
    ax = fig.add_subplot(1, 1, 1, projection='scatter_density')
    if temp_all[i] < 250: 
        density = ax.scatter_density([temp_all[i]] * len(Theta_all[i]), Theta_all[i], cmap = alpha_viridis, vmax = 20000, dpi=50)
    if temp_all[i] > 250: 
        density = ax.scatter_density([temp_all[i]] * len(Theta_all[i]), Theta_all[i], cmap = alpha_viridis, vmax = 20000, dpi=50)
    ax.set_facecolor('none')
    fig.colorbar(density, extend='max', label = r"Number of $\theta$")
    if temp_all[i] == 700: 
        ax.plot(Temp, Theta_L1, linestyle = ":", color = "black")
        ax.plot(Temp, Theta_L2, linestyle = ":", color = "black")
    ax.set_ylim(-0.75, 0.75)
    ax.set_xlim(0, 750)
    ax.set_ylabel(r"Angle, $\theta$ (rad)")
    ax.set_xlabel(r"Temp (K)")

plt.savefig("Theta_vs_Temp_density-allTemp-withAv_Graph.png", bbox_inches='tight', transparent=True)

