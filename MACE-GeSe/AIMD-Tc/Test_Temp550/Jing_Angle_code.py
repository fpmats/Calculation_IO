
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

from ase import Atoms
from ase.io.trajectory import Trajectory
from ase.io.vasp import write_vasp
from ase.io.xyz import write_xyz
from ase.io import read
from ase.io import write
from os import system

import sys
## Import path to the monochalcogenpy code
sys.path.insert(0, '/home/ac.tmihm/Project/GeSe/Orthorombic/ML-ForceField/Original/Run_MD_PVT_TempRange-NoVDW+D3/VASP-OPT-POSCAR')

from monochalcogenpy.tilt_angle import tilt_angle

#########################################
##  System info  ##
########################################

traj = Trajectory('traj/md.traj')

System = "Original_MD_600K_1000ps"
psmax = 1000
Tmax =600
#start_ps = 3
start_ps = 1500

#########################################
##  pull data from trajectory  ##
########################################

Angle = []
Av_angle = []

for i in range((len(traj) - start_ps), len(traj)):
    atoms = traj[i]
    #key = Ge index, values = list of tilt angles in the ac plane pointing to c within 3 Ang radius

    #tilt_angle_dict = tilt_angle(atoms, 'ac')
    #to get all tilt angles

    tilt_angle_array_ac = np.concatenate(list(tilt_angle(atoms, 'ac').values())) 
    print("ac")
    print(tilt_angle_array_ac)
    TA_list_ac = tilt_angle_array_ac.tolist()

    tilt_angle_array = np.concatenate(list(tilt_angle(atoms, 'bc').values())) 
    print("bc")
    #print(tilt_angle_array)
    TA_list = tilt_angle_array.tolist()
    Angle += TA_list
    print("len",len(Angle))

    #to get cell average

    tilt_angle_array_av_ac = np.mean(TA_list_ac)
    #print("Aver title angle ac")
    print(tilt_angle_array_av_ac)

    tilt_angle_array_av = np.mean(TA_list)
    print("Aver title angle bc")
    print(tilt_angle_array_av)
    Av_angle +=[float(tilt_angle_array_av)]

print("average angle", Av_angle)
print("all angles")
print(Angle)

#####################
## Av List        ##
####################

Av_av_angle = np.mean(Av_angle)
Av_Angle = np.mean(Angle)

print("average of the Av angle (deg):", Av_av_angle)
print("average of all angles (deg):", Av_Angle)

