#!/bin/bash

## Allocate resources
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=all
#SBATCH --job-name="boost-hic"
##SBATCH --mail-user=your@email
##SBATCH --mail-type=fail,end

source <your_home>/miniconda3/bin/activate boost-hic

echo "CALL: Boost-HiC/main.py $@"
python3 -u <your_instalation_home>/Boost-HiC/boost-hic.py "$@"
