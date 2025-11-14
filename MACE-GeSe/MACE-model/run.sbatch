#!/bin/bash

#PBS -N mace-test
#PBS -A NNFF_mem_sys
#PBS -l select=1:ngpus=1
###PBS -M tmihm@bu.edu # Optional if you require email
###PBS -m bea                 # Optional if you require email
#PBS -l walltime=00:30:00

#export PATH=$PATH:/lcrc/project/2D_dynamics/jing/software/miniconda3/bin/
#which conda
#echo path: $PATH

cd $PBS_O_WORKDIR
echo Working directory is $PBS_O_WORKDIR

#source activate /lcrc/project/2D_dynamics/jing/software/miniconda3/envs/mace-general/
#source activate /lcrc/project/2D_dynamics/jing/software/miniconda3/envs/NNFF-MACE_main/

#export DATA=/lcrc/project/NNFF_mem_sys/workflow_share/01_mace_training/dataset/
#export DATA=/lcrc/project/2D_dynamics/jing/dataset/nickelates/test_trains/combine_AIMD_metastable/old/

conda env list
##conda init
#conda activate ase_mace_env
source activate /home/ac.tmihm/miniconda3/envs/ase_mace_env
echo $CONDA_DEFAULT_ENV

##PATH="/home/ac.tmihm/GeSe/MACE-FF/Original-AvSq-ISIF8/01_mace_training"

mace_run_train \
    --name="model" \
    --train_file="dataset/train.xyz" \
    --test_file="dataset/test.xyz" \
    --energy_key="energy" \
    --forces_key="forces" \
    --stress_key="stress"  \
    --energy_weight=10.0 \
    --forces_weight=1000.0 \
    --stress_weight=10.0 \
    --compute_stress=1 \
    --lr=0.01 \
    --save_cpu \
    --model="MACE" \
    --num_interactions=2 \
    --max_ell=2 \
    --hidden_irreps="64x0e + 64x1o + 64x2e" \
    --num_cutoff_basis=5 \
    --correlation=2 \
    --r_max=4.0 \
    --batch_size=12 \
    --valid_batch_size=10 \
    --eval_interval=1 \
    --max_num_epochs=240 \
    --start_swa=30 \
    --swa_energy_weight=1000.0 \
    --swa_forces_weight=100.0 \
    --swa_stress_weight=1000.0 \
    --swa_lr=0.001 \
    --ema \
    --ema_decay=0.99 \
    --amsgrad \
    --error_table="PerAtomRMSE" \
    --default_dtype="float32" \
    --keep_checkpoints \
    --swa \
    --device=cuda \
    --seed=1234 \
