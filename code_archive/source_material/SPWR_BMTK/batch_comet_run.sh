#!/bin/bash

#SBATCH --partition compute
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=20
#SBATCH -A TG-DBS180005
#SBATCH --job-name=ca1
#SBATCH --output=ca1%j.out
#SBATCH --time 0-08:00
#SBATCH --qos=normal

module purge
module load python
module load intel
module load openmpi_ib
export PATH=/home/latimerb/neuron/nrn/x86_64/bin:$PATH
export LD_LIBRARY_PATH=/home/latimerb/neuron/nrn/x86_64/lib:$LD_LIBRARY_PATH


rm -rf output


echo "Running model at $(date)"

#mpirun nrniv -mpi -quiet -python3 run_network.py simulation_config.json
ibrun nrniv -mpi -python run_network.py simulation_configSAVEWGTS.json

echo "Done running model at $(date)"
