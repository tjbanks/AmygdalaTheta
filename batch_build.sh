#!/bin/bash

#SBATCH -N 1
#SBATCH -n 1
#SBATCH --qos=normal
#SBATCH --job-name=amygdala_theta_build
#SBATCH --output=amygdala_batch_%j.out
#SBATCH --time 0-2400:00

START=$(date)
python build_network.py
END=$(date)


{ printf "Start: $START \nEnd:   $END\n"; }| mail -r tbg28@mail.missouri.edu -s "AmygdalaTheta Build Complete" tbg28@mail.missouri.edu

echo "Done running build at $(date)"
