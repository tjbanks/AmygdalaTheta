#!/bin/bash

#SBATCH --partition shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --account=umc113
#SBATCH --job-name=amygdalatheta-tbanks
#SBATCH --output=amygdalatheta_%j.out
#SBATCH --mem=24G
#SBATCH --time 0-48:00

START=$(date)
#mpiexec nrniv -mpi -python run_network.py simulation_configECP_base_homogenous.json
mpiexec ./components_homogenous/mechanisms/x86_64/special -mpi run_network.py simulation_configECP_base_homogenous.json
END=$(date)

{ printf "Start: $START \nEnd:   $END\n" & python analysis_hom.py --save-plots & sleep 20 & printf "\n\n" & git diff components_homogenous/synaptic_models/; }| mail -r tbg28@mail.missouri.edu -s "AmygdalaTheta Simulation Results" -a analysis.png tbg28@mail.missouri.edu

echo "Done running model at $(date)"
