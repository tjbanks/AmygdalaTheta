#!/bin/bash

#SBATCH --partition shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --account=umc113
#SBATCH --job-name=amygdalatheta-tbanks
#SBATCH --output=amygdalatheta_%j.out
#SBATCH --mem=24G
#SBATCH --time 0-48:00

module purge
module load slurm
module load cpu
module load gcc
module load openmpi
module load ncurses

rm -rf outputECP/*

echo "Running model at $(date)"

mpiexec ./components/mechanisms/x86_64/special -mpi run_network.py simulation_configECP_base.json

echo "Done running model at $(date)"
