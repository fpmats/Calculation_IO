#name: Tina Mihm
#Date: May 8, 2025
#Description: reads in the energy data from each training set and graphs as a histogram

import os
from mpl_toolkits.mplot3d import axes3d
import numpy as np
import matplotlib.pyplot as plt
from ase import Atoms
from ase.io import read, write
import pandas as pd
from ase.geometry import wrap_positions
from ase.io.vasp import write_vasp
#-------------------------------------------------------------------------------------------
from pylab import *
from matplotlib.colors import *
from matplotlib.font_manager import fontManager, FontProperties
from scipy.optimize import curve_fit
from matplotlib.ticker import MaxNLocator
from scipy.stats import boltzmann
import seaborn as sns
from scipy.stats import gaussian_kde
from sklearn.neighbors import KernelDensity

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
rcParams['legend.fontsize'] = 12
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
ratio=1                     # square figure
plt.rcParams["figure.figsize"] = 3.37, 3.37*ratio
# plt.rcParams["figure.figsize"] = 3.37*ratio, 3.37
# fig = figure()

#---------------------------------------------------------
## System Data
#---------------------------------------------------------
System = "Training-Set-"
MACE_model = "MACErcut4.0-"

atom = 72

#---------------------------------------------------------
## Reads in the csv data 
#---------------------------------------------------------
fig = plt.figure()

df = pd.read_excel("AllDistortions-Energy-and-lattice-data.xlsx", "Dist_00")

#-------------------------------------------------------------------------------
## Loop over excel and read in data for each MD run for the distorted structures
#-------------------------------------------------------------------------------
#list for numbered distortion in excel 
Dist = [0, 1, 2, 5, 7, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 21, 22, 23, 24]

E = []
E_tot = []
L_tot = 0
D_v = []
E_max = 0
E_min = 0

for d in Dist: 
    if d < 10: 
        df = pd.read_excel("AllDistortions-Energy-and-lattice-data.xlsx", "Dist_0"+str(d))
    else:
        df = pd.read_excel("AllDistortions-Energy-and-lattice-data.xlsx", "Dist_"+str(d)) 

    E += [df["Energy (eV/atom)"]]
    Max = max(df["Energy (eV/atom)"])
    Min = min(df["Energy (eV/atom)"])
    if d == 0: 
        E_max = Max
        E_min = Min
    if Max > E_max: 
        E_max = Max
    if Min < E_min: 
        E_min = Min
    D_v += [d]
    # Add all the energies into one array to get the total distribution
    L_tot = L_tot + len(df["Energy (eV/atom)"])
    E_tot.extend(df["Energy (eV/atom)"])

print("number of structures in training set:", len(E_tot))

#-------------------------------------------------------------------------------
## Fit data to lorentzian fit vs a gaussian fit
#-------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import math

def lorentzian(x, x0, gamma):
    """Single Lorentzian function."""
    return (gamma / np.pi) / ((x - x0)**2 + gamma**2)
def sum_of_lorentzians(x, values, gamma):
    """Sum of Lorentzian functions centered at each value."""
    spectrum = np.zeros_like(x)
    for v in values:
        spectrum += lorentzian(x, v, gamma)
    return spectrum
#---------------------------------------------------

# # Example usage:
# if __name__ == "__main__":
#     # Discrete values to broaden
#     values = [1.0, 2.0, 3.0, 5.0]
#     # Lorentzian broadening
#     gamma = 0.1
#     # Evaluation grid
#     x = np.linspace(0, 6, 1000)
#     # Compute broadened spectrum
#     spectrum = sum_of_lorentzians(x, values, gamma)
#     # Plotting
#     plt.plot(x, spectrum, label='Sum of Lorentzians')
#     plt.xlabel('x')
#     plt.ylabel('Intensity')
#     plt.title('Lorentzian Broadening of Discrete Values')
#     plt.grid(True)
#     plt.legend()
#     plt.show()

if __name__ == "__main__":
    # Discrete values to broaden
    values = E_tot
    # Lorentzian broadening
    gamma = 0.01
    # Evaluation grid
    x2 = np.linspace(0, 2, 1000)
    print("E_min, E_max", E_min, E_max)

    ## ground state energy for the optimized a=b structure without distortions 
    Grn_E = -35.453553/8
    print("ground state E", Grn_E)
    # shift energies by the ground state 
    x_sub = np.array([i + Grn_E for i in x2])
    print("x max:", max(x_sub))

    # Compute broadened spectrum
    spectrum = sum_of_lorentzians(x_sub, values, gamma)

    #find the max peak value and divide the data by this to norm to one
    max_s = max(spectrum)

    # normalize the data: 

    spectrum_norm = (spectrum[0]/2 + np.sum(spectrum[1:-1]) + spectrum[-1]/2) * (x2[1]-x2[0])
    print("normalization value for spectrum", spectrum_norm)
    print("one over normalization value for spectrum", 1/spectrum_norm)
    spectrum = (1/spectrum_norm)*spectrum

    ## normalized Boltzmann Distributions (BD): 
    x_sub_max300 = [(1/0.025)*np.exp(-i/0.025) for i in x2] ## T = 300
    x_sub_max600 = [(1/0.05)*np.exp(-i/0.05) for i in x2] ## T = 600
    x_sub_max900 = [(1/0.075)*np.exp(-i/0.075) for i in x2] ## T = 900

    a = 0
    b = 0
    c = 0

    # counts up the number of structures that fall between certain ranges on the graph
    diff = [ i - Grn_E for i in E_tot]
    for i in diff:
        ## between T = 300 and T = 900 
        if i < 0.6 and i > 0.2: 
            b = b+1
        # at T = 300
        if i < 0.2: 
            a = a+1
        ## beyond T = 900
        if i > 0.6: 
            c = c+1
    print("Number of structure at T = 300:", a)
    print("Number of structure at between T = 300 and T = 900:", b)
    print("Number of structure above T = 900:", c)

    ## Test normalization results in area under the curve equal to 1: 
    from scipy import integrate
    T300 = integrate.trapezoid(x_sub_max300, x2)
    T600 = integrate.trapezoid(x_sub_max600, x2)
    T900 = integrate.trapezoid(x_sub_max900, x2)
    Training = integrate.trapezoid(spectrum, x_sub)
    print("Result T300 (Trapezoidal):", T300)
    print("Result T600 (Trapezoidal):", T600)
    print("Result T900 (Trapezoidal):", T900)
    print("Result Training Set (Trapezoidal):", Training)
    
    #-----------------------------------------------------
    plt.figure(4)
    fig, ax = plt.subplots()
    # Plotting
    plt.plot(x2, x_sub_max300, label='300K', color = "#ff9999", linestyle = ":")
    plt.plot(x2, x_sub_max600, label='600K', color = "#ff0000", linestyle = ":")
    plt.plot(x2, x_sub_max900, label='900K', color = "#990000", linestyle = ":")

    plt.plot(x2, spectrum, label='Training set', color = "#0000FF")

    plt.xlabel(r"$\Delta E$ (eV/atom)")
    plt.ylabel('Frequency')
    plt.yscale("log")
    plt.ylim(10**(-2.5), 50)
    plt.xlim(0, 0.75)

    leg = plt.legend(loc='center', bbox_to_anchor=(0.7, 0.83), ncol=1, handlelength = 2.0, handletextpad = 0.5, columnspacing=0.5, fontsize = 10, borderpad=0.1, frameon=False)

    fig.savefig("GeSe-TrainingData_Analysis-Total-"+System+MACE_model+"-Loren_density-log.png", bbox_inches='tight')
