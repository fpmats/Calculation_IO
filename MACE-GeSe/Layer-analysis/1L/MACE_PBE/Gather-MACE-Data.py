#Name: Tina Mihm
#Date: Nov 13, 2024
#Description: Pulls the MACE energies from MACE.out and any info from the directory name and puts it into a csv for data analysis

import os
import numpy as np
from ase import Atoms
from ase.io import read, write
import matplotlib.pyplot as plt
import pandas as pd
from ase.geometry import wrap_positions
from ase.io.vasp import write_vasp
import glob
import re

#---------------------------------------------------------
## loops through directories and pulls data from MACE.out 
#---------------------------------------------------------

System = "Non-VDW-Layer-Analysis-Origianl-1L-Streached-"
MACE_model = "MACErcut4.0-"
directory_path = "./"
dir_name = "Test_Theta*"
filename = "MACE.out"
string = "Total Energy from MACE:"

E = []
Diff = []
#Layer = []

Dir = glob.glob(dir_name)
print(Dir)

for file in Dir:
    file_path = os.path.join(directory_path, file)
    # Check if it's a file (not a directory)
    if os.path.isdir(file_path):
        coord = re.findall(r"\d+", file)
        #coord2 = re.findall(r"\d+\.\d+", file)
        print(coord)
        #print(coord2)
        Diff+=[float(coord[0])]
        #Layer+=[float(coord2[0])]
        print(file_path)
        for fl in os.listdir(file_path):
            if fl == filename:
                print("Found the", fl)
                file_path2 = os.path.join(file_path, fl)
                print(file_path2)
                with open(file_path2, 'r') as fp:
                    # read all lines using readline()
                    lines = fp.readlines()
                    for row in lines:
                        # check if string present on a current line
                        if row.find(string) != -1:
                            print(row)
                            s = re.search(r"-?\d+\.\d+", row)
                            en = float(s.group())
                            print(en)
                            E += [en]
print(Diff, E)

df = {"Difference added": Diff, "Energy (eV)": E}
df = pd.DataFrame(df)
df = df.sort_values(by=["Difference added"])
print(df)
df.to_csv(System+MACE_model+"Data.csv", index = False)
