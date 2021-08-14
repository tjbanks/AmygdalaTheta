#!/bin/bash
#SBATCH --partition Lewis
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --qos=normal
#SBATCH --job-name=ca1
#SBATCH --output=ca1%j.out
#SBATCH --time 0-00:30



echo "Building model at $(date)"

python3 build_network.py

echo "Done building model at $(date)"

