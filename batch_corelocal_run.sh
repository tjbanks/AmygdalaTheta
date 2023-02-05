#!/bin/bash

#SBATCH -N 1
#SBATCH -n 25
#SBATCH --qos=normal
#SBATCH --job-name=amygdala_theta
#SBATCH --output=amygdala_batch_%j.out
#SBATCH --time 0-12:00

START=$(date)
mpiexec ./components/mechanisms/x86_64/special -mpi run_network.py simulation_configECP_base.json
END=$(date)

{ printf "Start: $START \nEnd:   $END\n" & python analysis.py --save-plots & printf "\n\n" & git diff components/synaptic_models/; }| mail -r tbg28@mail.missouri.edu -s "AmygdalaTheta Simulation Results" -a analysis.png tbg28@mail.missouri.edu

echo "Done running model at $(date)"
