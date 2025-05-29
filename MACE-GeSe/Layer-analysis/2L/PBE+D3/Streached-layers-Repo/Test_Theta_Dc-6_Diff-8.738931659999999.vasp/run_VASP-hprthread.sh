#!/bin/bash -l
#$ -l h_rt=6:00:00
#$ -N GeSe-MP-vdW-k332-VASP
##$ -q qonos
#$ -P fpmats
#$ -m bea
#$ -M tmihm@bu.edu #change it to your email here or delete this line
#$ -pe mpi_16_tasks_per_node 32
##$ -pe mpi_64_tasks_per_node 128
##$ -l mem_per_core=16G

module use /project/fpmats/software/production/modules
module load openmpi/4.1.5_nvidia-2023-23.5
module load vasp/6.4.2
module list 2>&1

##touch CHGCAR WAVECAR HILLSPOT ICONST

# if multi-threading, uncomment and set this
#export OMP_NUM_THREADS=4

export VASP_PP_PATH="/projectnb/fpmats/jing/VASPPAW"
#export ASE_VASP_COMMAND="mpirun --map-by ppr:8:socket:PE=4 --bind-to core  vasp_std > vasp.out"
export ASE_VASP_COMMAND="mpirun vasp_std > vasp.out"


###copy charge density and wavefunction files from previous run
#if you need to
#mkdir 01_AOPT  
mkdir 02_SCF
#mkdir 03_BS
#mkdir 03_BS-ISPIN1

#oldDIR="01_AOPT"
#oldDIR="02_SCF"
#oldDIR="02_SCF-ISPIN1"

#DIR="01_AOPT" 
DIR="02_SCF"
#DIR="03_BS"
#DIR="03_BS-ISPIN1"

#cp $oldDIR/CHGCAR $DIR/.
#cp $oldDIR/WAVECAR $DIR/.

#cp POTCAR $DIR/.

#This is needed for the very first calculation if the 
# CHGCAR and WAVECAR are not being copied from somewhere else
cd $DIR
touch CHGCAR WAVECAR HILLSPOT ICONST
cd ../

#this kernel is required for vdw-df runs
##cp vdw_kernel.bindat $DIR

#run your script
#python ASE-Generate-Strained-Lattice-constants.py 
#python ASE-VASP-Atom_Opt.py
#python ASE-VASP-SCF.py
python ASE-VASP-SCF-Torch-VDW.py > VASP.out
#python ASE-VASP-Bandstructure.py
