#!/bin/bash

#SBATCH -N 1
#SBATCH -n 50
#SBATCH --qos=normal
#SBATCH --job-name=amygdala_theta
#SBATCH --output=amygdala_batch_%j.out
#SBATCH --time 0-24:00

START=$(date)
mpiexec nrniv -mpi -quiet -python run_network.py simulation_configECP_vpsi_homogenous_ach.json
END=$(date)

{ printf "Start: $START \nEnd:   $END\n" & python analysis_hom.py --save-plots & printf "\n\n" & git diff components_homogenous/synaptic_models/; }| mail -r tbg28@mail.missouri.edu -s "AmygdalaTheta VPSI ACH Simulation Results" -a analysis.png tbg28@mail.missouri.edu




