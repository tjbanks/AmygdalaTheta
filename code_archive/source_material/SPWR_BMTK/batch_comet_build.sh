#!/bin/bash
#SBATCH --partition compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -A TG-DBS180005
#SBATCH --job-name=ca1
#SBATCH --output=ca1%j.out
#SBATCH --time 0-00:30
#SBATCH --qos=normal


module purge
module load python
module load intel
module load openmpi_ib
export PATH=/home/latimerb/neuron/nrn/x86_64/bin:$PATH
export LD_LIBRARY_PATH=/home/latimerb/neuron/nrn/x86_64/lib:$LD_LIBRARY_PATH


rm -rf network

echo "Building model at $(date)"

python3 build_network.py

echo "Done building model at $(date)"


