#!/bin/bash

#PBS -N mace-test
#PBS -A NNFF_mem_sys
#PBS -l select=1:ngpus=1
###PBS -M YOUR-EMAIL# Optional if you require email
###PBS -m bea                 # Optional if you require email
#PBS -l walltime=00:30:00

cd $PBS_O_WORKDIR
echo Working directory is $PBS_O_WORKDIR

conda env list
source activate /home/ac.tmihm/miniconda3/envs/ase_mace_env
echo $CONDA_DEFAULT_ENV

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
