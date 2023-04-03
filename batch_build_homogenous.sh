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
python build_network.py homogenous
END=$(date)

echo "Done running build at $(date)"
