from ase.io import read
from ase.visualize import view
import numpy as np
from ase.build import make_supercell
from ase.calculators.vasp import Vasp
from mace.calculators import MACECalculator
from torch_dftd.torch_dftd3_calculator import TorchDFTD3Calculator
from ase.calculators.mixing import SumCalculator
from ase import Atoms
from ase.db import connect
import sys
import os
import re
from itertools import product
from ase.optimize import LBFGS
from ase.constraints import FixSymmetry

os.system("pwd")

path = './'

db = connect('./relax-11A.db')

def relax(atoms, file_name):
    opt = LBFGS(atoms,
                logfile='./log/'+file_name+'.log',
                trajectory='./traj/'+file_name+'.traj',
                )
    opt.run(fmax=0.01, steps=250)
    atoms.info['opt_converged'] = opt.converged()

## reads in data from POSCAR and generate data for the MACE calculations in the database
pwd = "../../../POSCAR_registry/"
#cwd = os.getcwd()
#print(cwd)
#directory =pwd+"/C-10A"
directory =pwd+"/11A_POSCAR"
print(directory)

for file in os.listdir(directory):
    print(file) 
    if "_POS" in file: 
        continue
    coord = re.findall(r"\d+\.\d+", file)
    print(coord)
    dA = float(coord[0])
    dB = float(coord[1])
    atoms = read(directory+"/"+file)
    pos = atoms.get_positions()
    cell = atoms.get_cell()
    a = cell[0,0]
    b = cell[1,1]
    c = cell[2,2]
    run_name = f'{dA:.02f}_{dB:.02f}_{c:.02f}_opt_mace'    

    mace_calc = MACECalculator(model_paths='./model_swa.model', device='cpu')
    d3 = TorchDFTD3Calculator(atoms=atoms,damping="zero",dispersion_xc = "pbe",dispersion_cutoff = 40.0 ,device = "cpu")
    calc = SumCalculator([mace_calc, d3])

    #bc_rev_atoms.calc = Mixcalc
    #atoms.set_constraint([FixSymmetry(atoms)])
    atoms.calc = calc
    print("***Total Energy from MACE:", dA, dB, c, atoms.get_potential_energy())
    db.write(atoms,a=a, b=b, c=c,dA=dA, dB=dB, model='macecalculator')
    #try:
        #db.get(a=a,b=b,c=c,dA=dA,dB=dB,opt_converged=True, model='macecalculator')
        #db.get(a=a,b=b,c=c,dA=dA,dB=dB, model='macecalculator')
    #except KeyError:
        #relax(atoms, run_name)
        #print("***Total Energy from MACE:", dA, dB, c, atoms.get_potential_energy())
        #db.write(atoms,a=a, b=b, c=c,dA=dA, dB=dB, model='macecalculator', opt_converged=atoms.info['opt_converged'])       
        #db.write(atoms,a=a, b=b, c=c,dA=dA, dB=dB, model='macecalculator')       


