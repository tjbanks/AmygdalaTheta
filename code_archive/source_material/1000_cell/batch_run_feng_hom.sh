#!/bin/bash

#SBATCH -N 1
#SBATCH -n 50
#SBATCH --qos=normal
#SBATCH --job-name=feng_1000
#SBATCH --output=feng1000_%j.out
#SBATCH --time 0-120:00

mpirun nrniv -mpi MC_main_small_forGamma.hoc 
