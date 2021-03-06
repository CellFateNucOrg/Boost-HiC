#!/bin/bash

## Allocate resources
#SBATCH --time=5-00:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --cpus-per-task=6
#SBATCH --partition=all
#SBATCH --job-name="boost-hic"
##SBATCH --mail-user=your@email
##SBATCH --mail-type=fail,end

# This line have to be updated to your environment
source ~/miniconda3/bin/activate boost-hic

echo "CALL: Boost-HiC/boost-hic.py $@"
# This line have to be updated to your environment
python3 -u ~/software/github/Boost-HiC/boost-hic.py "$@"
