#!/bin/bash
#SBATCH --qos=regular
#SBATCH --nodes=2
###SBATCH --tasks-per-node=1
#SBATCH --account=m2120
#SBATCH -t 01:30:00
#SBATCH --job-name=BGW-QE-02-mergeWF
#SBATCH -C cpu
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tmihm@bu.edu

export OMP_NUM_THREADS=1

module load berkeleygw/3.0.1-cpu

cat kpnt1/vxc.dat kpnt2/vxc.dat >vxc.dat
cp kpnt1/rho.complex ./

srun wfnmerge.x &> wfnmerge.out
