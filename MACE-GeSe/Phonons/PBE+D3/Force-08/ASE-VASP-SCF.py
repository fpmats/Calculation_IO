from ase.io import read
from ase.visualize import view
import numpy as np
from ase.build import make_supercell
from ase.calculators.vasp import Vasp
from ase import Atoms
import sys

path = './'
charge=0
#atoms=read('01_AOPT/CONTCAR')
atoms=read('POSCAR')

calc_scf = Vasp(system = 'SCF',
             directory='./',
             prec = 'Accurate',
             charge = charge,
             istart = 1,
             icharg = 0,
             ismear = 0,
             sigma = 0.05,
             ispin = 2,
             #npar = 9,
             #kpar = 6,
             kpts = (1,1,1),
             gamma = True,
             lwave = True,
             lcharg = True,
             lreal = False,
             lasph = True,
             algo = 'N',
             ediff = 1e-6,
             nbands = 1728,
             xc = 'PBE',
             encut = 520,
             nelmin = 5, 
             nelm = 100, 
             ivdw=11,
             lorbit = 11)

atoms.calc = calc_scf
atoms.get_potential_energy()


