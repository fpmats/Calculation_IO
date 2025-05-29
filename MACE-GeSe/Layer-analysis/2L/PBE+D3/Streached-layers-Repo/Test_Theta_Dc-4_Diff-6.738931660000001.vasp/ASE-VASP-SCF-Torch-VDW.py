from ase.io import read
from ase.visualize import view
import numpy as np
from ase.build import make_supercell
from ase.calculators.vasp import Vasp
from torch_dftd.torch_dftd3_calculator import TorchDFTD3Calculator
from ase.calculators.mixing import SumCalculator
from ase import Atoms
import sys
import os

path = './'
charge=0
#atoms=read('01_AOPT/CONTCAR')
atoms=read('POSCAR')

calc_scf = Vasp(system = '02_SCF',
             directory='02_SCF',
             prec = 'Accurate',
             charge = charge,
             istart = 1,
             icharg = 0,
             ismear = 0,
             sigma = 0.05,
             ispin = 2,
             #npar = 9,
             kpar = 2,
             kpts = (6,6,2),
             gamma = True,
             lwave = True,
             lcharg = True,
             lreal = False,
             lasph = True,
             algo = 'N',
             ediff = 1e-6,
             nbands = 220,
             xc = 'PBE',
             encut = 520,
             nelmin = 5, 
             nelm = 100, 
             #ivdw=11,
             lorbit = 11)


d3 = TorchDFTD3Calculator(atoms=atoms,damping="zero",dispersion_xc = "pbe",dispersion_cutoff = 40.0 ,device = "cpu")
Mixcalc = SumCalculator([calc_scf, d3])
atoms.calc = Mixcalc
atoms.get_potential_energy()
print("Total Energy from VASP:", atoms.get_potential_energy())


