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

db = connect('../relax.db')

def relax(atoms, file_name):
    opt = LBFGS(atoms,
                logfile='../../log/'+file_name+'.log',
                trajectory='../../traj/'+file_name+'.traj',
                )
    opt.run(fmax=0.01, steps=250)
    atoms.info['opt_converged'] = opt.converged()

def vasp_calc(vdw=True, run_name='run'):
    #os.system("mkdir "+run_name)
    #os.system("cd "+run_name)
    #os.system("touch CHGCAR WAVECAR HILLSPOT ICONST")
    #os.system("cd ../../")
    vasp_calc = Vasp(
        system = run_name,
        directory=run_name,
        prec = 'Normal',
        istart = 0,
        icharg = 2,
        isym = 0,
        encut = 360,
        ismear = 0,
        sigma = 0.05,
        npar = 16,
        kpar = 2,
        kpts = (6,6,2),
        #kspacing = 0.25,
        #symprec = 1e-4,
        gamma = True,
        lwave = False,
        lcharg = False,
        lmaxmix = 4,
        xc = 'PBE',
        lasph = True,
        #algo = 'Fast',
        algo = 'Normal', #for fu1.2 R5+
        ediff = 1e-6,
        ispin = 1,
        nelm = 100,
    )
    device = 'cpu'

    calc_d3 = TorchDFTD3Calculator(
            damping="zero",
            dispersion_xc = "pbe",
            dispersion_cutoff = 40.0,
            device = device,
            )
    calc = SumCalculator([vasp_calc, calc_d3])
    #atoms.info['opt_converged'] = opt.converged()
    return calc


## reads in data from POSCAR and generate data for the VASP calculations in the database

cwd = os.getcwd()
print(cwd)
#directory =cwd+"/POSCAR_shrink-F"
 
os.system("touch CHGCAR WAVECAR HILLSPOT ICONST")
atoms = read("POSCAR")
pos = atoms.get_positions()
cell = atoms.get_cell()
a = cell[0,0]
b = cell[1,1]
c = cell[2,2]
run_name = f'{a:.02f}_{b:.02f}_{c:.02f}_opt_mace'    

#bc_rev_atoms.calc = Mixcalc
atoms.set_constraint([FixSymmetry(atoms)])
atoms.calc = vasp_calc(vdw=True, run_name= "./")    
print("***Total Energy from VASP:", b, c, atoms.get_potential_energy())
try:
    db.get(a=a,b=b,c=c, model='vasp',opt_converged=True)
except KeyError:
    relax(atoms, run_name)
    db.write(atoms, a=a, b=b, c=c, model='vasp', opt_converged=atoms.info['opt_converged'])   


