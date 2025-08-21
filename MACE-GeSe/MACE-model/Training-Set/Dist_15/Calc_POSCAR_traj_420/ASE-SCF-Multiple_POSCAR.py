from ase.io import read
from ase.visualize import view
import numpy as np
from ase.build import make_supercell
from ase.calculators.vasp import Vasp
from ase import Atoms
from ase.calculators.emt import EMT
from ase.io.trajectory import Trajectory
import os

directory = "./POSCARs/"
#atoms=read('POSCAR')

calc = Vasp(system = '01_AOPT',
            directory='./',
            prec = 'Accurate',
            istart = 1,
            icharg = 0,
            ismear = 0,
            sigma = 0.05,
            ispin = 2,
            npar = 8,
            kpar = 2,
            kpts = (2,2,2),
            gamma = True,
            lwave = False,
            lcharg = False,
            lreal = False,
            lasph = True,
            algo = 'N',
            ediff = 1e-6,
            xc = 'PBE',
            encut = 520,
            nelmin = 5,
            nbands = 220,
            isym = -1,
            #ivdw = 11,
            nelm = 100,
            lorbit = 11)

#atoms.calc = calc_AOp

traj = Trajectory('GeSe_Dist15.traj', 'w')
for name in os.listdir(directory):
    print(name)
    atoms=read(directory+name)
    atoms.calc = calc
    atoms.get_potential_energy()
    os.system("mkdir Calc_"+name)
    os.system("cp * Calc_"+name+"/.")
    traj.write(atoms)
