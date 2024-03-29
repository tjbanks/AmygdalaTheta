#!/bin/bash

#SBATCH -N 1
#SBATCH -n 24
#SBATCH --qos=normal
#SBATCH --job-name=amygdala_theta
#SBATCH --output=amygdala_batch_%j.out
#SBATCH --time 0-24:00

START=$(date)
mpiexec nrniv -mpi -quiet -python run_network.py simulation_configECP_vpsi_homogenous_ach_high.json
#mpiexec ./components_homogenous/mechanisms/x86_64/special -mpi run_network.py simulation_configECP_vpsi_homogenous_ach_high.json
END=$(date)

printf "Start: $START \nEnd:   $END\n"
python analysis_hom.py --save-plots



