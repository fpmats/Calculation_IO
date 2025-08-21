from ase.io import read, write
from glob import glob
import re

for traj_name in glob('Test_Temp*/traj/md.traj'):
    traj = read(traj_name,'::200')
    T = int(re.findall(r'\d+', traj_name)[0]) #now updated
    write(f'Trunc_Traj/MD_{T:d}K.traj', traj)
