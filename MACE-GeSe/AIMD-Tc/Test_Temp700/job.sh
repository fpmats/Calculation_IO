#!/bin/bash
#PBS -l nodes=1:ppn=12:gpus=1
#PBS -l walltime=36:00:00
#PBS -A cnm83533
#PBS -N GeSe-ZG-SCF
#PBS -m bea

cd $PBS_O_WORKDIR

echo "Jobid:	${PBS_JOBID}"
echo "Dir:	${PBS_O_WORKDIR}"
echo "Nodes:"
uniq -c $PBS_NODEFILE

source activate $HOME/miniconda3/envs/ase-mace-cudaTorch
conda  activate $HOME/miniconda3/envs/ase-mace-cudaTorch

mkdir log
mkdir traj

#python md.py
python md_stage.py

##for i in Test*_2band*
##do
##    (	# run in subshell so "cd .." is not needed and to not get lost in the hierarchy
##	cd $i
##	mkdir 02_SCF
##	python ASE-VASP-SCF.py
##    )
##done
