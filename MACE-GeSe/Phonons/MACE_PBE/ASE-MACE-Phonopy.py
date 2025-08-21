from ase.io import read
from ase.visualize import view
import numpy as np
import pandas as pd
from ase.build import make_supercell
from ase.calculators.vasp import Vasp
import ase.io.vasp as vasp
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

from ase.filters import ExpCellFilter
from ase.phonons import Phonons


#----------------------------------------------------------------------------

cwd = os.getcwd()
print(cwd)
directory =cwd

def relax(atoms, file_name):
    opt = LBFGS(atoms,
                logfile='./log/'+file_name+'.log',
                trajectory='./traj/'+file_name+'.traj',
                )
    opt.run(fmax=0.01, steps=250)
    #atoms.info['opt_converged'] = opt.converged()

# Setup structure and calculator
atoms = read(directory+"/POSCAR")

mace_calc = MACECalculator(model_paths='./model_swa.model', device='cpu')
#d3 = TorchDFTD3Calculator(atoms=atoms,damping="zero",dispersion_xc = "pbe",dispersion_cutoff = 40.0 ,device = "cpu")
calc = mace_calc

#  Optomize the structure using chosen calculator

#atoms.set_constraint([FixSymmetry(atoms)])
atoms.set_calculator(calc)
ecf = ExpCellFilter(atoms)

run_name = 'mace_noVDW'
relax(ecf, run_name)

# Phonon calculator
ph = Phonons(atoms, mace_calc, supercell=(6, 6, 2), delta=0.03)
ph.run()

# Read forces and assemble the dynamical matrix
ph.read(acoustic=True)
ph.clean()


#path = [[0, 0, 0], [0.5, 0.0, 0.0], [0.5, 0.5, 0.0],[0.0, 0.5, 0.0], [0, 0, 0], [0.0, 0.0, 0.5], [0.5, 0.0, 0.5],[0.5, 0.5, 0.5], [0.0, 0.5, 0.5], [0.0, 0.0, 0.5]]
#labels = ["$\\Gamma$","X", "S", "Y", "$\\Gamma$", "Z", "U", "R", "T", "Z"]

from ase.io.jsonio import read_json

path = atoms.cell.bandpath('GXSYGZURTZ', npoints=500)

#from ase.dft.kpoints import ibz_points, bandpath, special_paths, sc_special_points
#from ase.dft.kpoints import special_paths

#points = sc_special_points['orc']
#print("points", points)

#energy = ph.band_structure(path)
bs = ph.get_band_structure(path)

path.write('mybandpath.json')
bs.write('mybandstructure.json')


print("bs")
print(bs)

#ph.write_modes([i for i in K])

dos = ph.get_dos(kpts=(41, 41, 41)).sample_grid(npts=500, width=1e-5)
print(dos)


weights = dos.get_weights()
dos_e = dos.get_energies()
print("weights")
print(weights)
print("energy") 
print(dos_e)

import pandas as pd

dos_df = {"weight": weights, "energy": dos_e}
dos_df = pd.DataFrame(dos_df)
dos_df.to_csv("dos_data.csv")

# Plot the band structure and DOS:
import matplotlib.pyplot as plt  # noqa

fig = plt.figure(figsize=(7, 4))
ax = fig.add_axes([.12, .07, .67, .85])

emax = 0.03
bs.plot(ax=ax, emin=-0.004, emax=emax)

dosax = fig.add_axes([.8, .07, .17, .85])
dosax.fill_between(dos.get_weights(), dos.get_energies(), y2=0, color='grey',
                   edgecolor='k', lw=1)

#dosax.plot(dos.get_weights(), dos.get_energies(), color='grey')

dosax.set_ylim(-0.004, emax)
dosax.set_yticks([])
dosax.set_xticks([])
dosax.set_xlabel("DOS", fontsize=18)

fig.savefig('GeSe-MACE_noVDW_phonon.png')

