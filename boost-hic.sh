#!/bin/bash

## Allocate resources
#SBATCH --time=0-16:00:00
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=4
#SBATCH --partition=all
#SBATCH --job-name="boost-hic"
#SBATCH --mail-user=todor.gitchev@izb.unibe.ch
#SBATCH --mail-type=fail,end

source /home/pmeister/miniconda3/bin/activate boost-hic

echo "CALL: Boost-HiC/main.py $@"
python3 -u /home/pmeister/software/Boost-HiC/main.py "$@"
