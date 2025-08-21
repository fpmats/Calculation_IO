from ase.io.trajectory import Trajectory
from ase import Atoms
from ase.io.vasp import write_vasp
from ase.io.xyz import write_xyz
from ase.io import read
from ase.io import write
from os import system

System = "Dist07"
traj = Trajectory('GeSe_'+System+'.traj')

print(len(traj))
for i in range(0, len(traj)):
    write("Traj_info_"+str(i)+".xyz", traj[i], format='extxyz')
    system("cat Traj_info_"+str(i)+".xyz >> Traj_info_"+System+".xyz")
