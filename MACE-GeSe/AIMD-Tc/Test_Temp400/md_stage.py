from ase.io import read
import numpy as np
import sys
from mace.calculators import MACECalculator
from ase.md.langevin import Langevin
from ase.units import fs
from ase.md.nptberendsen import NPTBerendsen
from ase.md.nptberendsen import Inhomogeneous_NPTBerendsen
from ase.md import velocitydistribution
from ase.units import fs
from ase.calculators.vasp import Vasp
from ase.db import connect
import numpy as np
import sys
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
from time import perf_counter
import shutil
from ase.md import MDLogger
from torch_dftd.torch_dftd3_calculator import TorchDFTD3Calculator
from ase.calculators.mixing import SumCalculator


model_path = './'
model = "model_swa.model"

mace_calc = MACECalculator(model_paths=[model_path+model], device='cuda', default_dtype="float32")
mace_calc.parameters['path'] = f'{model_path}{model}'
mace_calc.parameters['name'] = f'{model}'

path = './'
file_name = 'POSCAR'
atoms = read(path+file_name)
# d3 = DFTD3(method="pbe", damping="d3zero").add_calculator(mace_calc)
calc_d3 = TorchDFTD3Calculator(
        damping="zero",
        xc = "pbe",
        cutoff = 15.0,
        device = "cuda",
        )
Mixcalc = SumCalculator([mace_calc, calc_d3])
# atoms.calc = mace_calc
atoms.calc = Mixcalc

# Print statements
def print_dyn():
    imd = dyn.get_number_of_steps()
    etot  = atoms.get_total_energy()
    temp_K = atoms.get_temperature()
    stress = atoms.get_stress(include_ideal_gas=True) #/GPa
    stress_ave = (stress[0]+stress[1]+stress[2])/3.0
    elapsed_time = perf_counter() - start_time
    print(f"  {imd: >3}   {etot:.3f}    {temp_K:.2f}    {stress_ave:.2f}  {stress[0]:.2f}  {stress[1]:.2f}  {stress[2]:.2f}  {stress[3]:.2f}  {stress[4]:.2f}  {stress[5]:.2f}    {elapsed_time:.3f}")

logfile = f'./log/md.log'

#MaxwellBoltzmannDistribution(atoms, temperature_K=25) ##adding randomize velocity with seed @ t = 0
MaxwellBoltzmannDistribution(atoms, temperature_K=400) ##adding randomize velocity with seed @ t = 0

dyn = Inhomogeneous_NPTBerendsen(atoms,   ## runs MD with Berendsen dynamics in NPT ensemble
               timestep = 10.0 * fs,  ## sets the time steps for the MD
               temperature_K=400,   ## temperature in Kelvin the MD is run at
               taut=100 * fs,
               pressure_au=6.32421e-7, ## atomospheric pressure in eV/A^3
               taup=100 * fs,
               compressibility_au=0.149796219,
               #logfile=f'./log/md.log',
               trajectory=f'./traj/md.traj',
               loginterval=50)

dyn.attach(print_dyn, interval=50)
dyn.attach(MDLogger(dyn, atoms, logfile, header=True, stress=True, peratom=True, mode="a"), interval=50)
start_time = perf_counter()
print(f"    imd     Etot(eV)    T(K)    stress(mean,xx,yy,zz,yz,xz,xy)(GPa)  elapsed_time(sec)")

dyn.run(100000)
