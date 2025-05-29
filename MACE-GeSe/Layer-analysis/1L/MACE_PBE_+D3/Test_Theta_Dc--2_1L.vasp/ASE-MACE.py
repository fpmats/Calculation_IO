from ase.io import read
from ase.visualize import view
import numpy as np
from ase.build import make_supercell
from ase.calculators.vasp import Vasp
from mace.calculators import MACECalculator
from torch_dftd.torch_dftd3_calculator import TorchDFTD3Calculator
from ase.calculators.mixing import SumCalculator
from ase import Atoms
import sys
import os

os.system("pwd")

path = './'

atoms=read('POSCAR')

mace_calc = MACECalculator(model_paths='./model_swa.model', device='cpu')
d3 = TorchDFTD3Calculator(atoms=atoms,damping="zero",dispersion_xc = "pbe",dispersion_cutoff = 40.0 ,device = "cpu")
Mixcalc = SumCalculator([mace_calc, d3])

atoms.calc = Mixcalc
#atoms.calc = calc
print("Total Energy from MACE:", atoms.get_potential_energy())


