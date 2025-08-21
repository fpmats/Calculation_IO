#Name: Tina Mihm
#Date: 04/15/2025
#Description: reads in the traj file and visualizes it in gui

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
from ase.visualize import view

traj = Trajectory('traj/md.traj')

view(traj)
