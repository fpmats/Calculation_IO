from ase.io import read
from ase.visualize import view
import numpy as np
from ase.build import make_supercell
from ase.calculators.vasp import Vasp
from mace.calculators import MACECalculator
from ase import Atoms
import sys
import os

os.system("pwd")

path = './'

atoms=read('POSCAR')

calc = MACECalculator(model_paths='./model_swa.model', device='cpu')
atoms.calc = calc

print("Total Energy from MACE:", atoms.get_potential_energy())


