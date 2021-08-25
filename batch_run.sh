#!/bin/bash

#SBATCH --partition Lewis
#SBATCH -N 1
#SBATCH -n 50
#SBATCH --qos=normal
#SBATCH --job-name=amygdala_theta
#SBATCH --output=ca1%j.out
#SBATCH --time 0-00:30

echo "Running AmygdalaTheta at $(date)"

mpiexec nrniv -mpi -quiet -python run_network.py simulation_configECP_base_edge_effects.json

{ python analysis.py -no-plots & git diff components/synaptic_models/; }| mail -r tbg28@mail.missouri.edu -s "AmygdalaTheta Simulation Results" tbg28@mail.missouri.edu

echo "Done running model at $(date)"
