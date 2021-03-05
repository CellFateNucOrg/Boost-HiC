#!/bin/bash

## Allocate resources
#SBATCH --time=5-00:00:00
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=4
#SBATCH --partition=all
#SBATCH --job-name="boost-hic"
#SBATCH --mail-user=todor.gitchev@izb.unibe.ch
#SBATCH --mail-type=fail,end

source /home/pmeister/miniconda3/bin/activate boost-hic

echo "CALL: Boost-HiC/main.py $@"
python3 -u /home/pmeister/software/github/Boost-HiC/main.py "$@"

## example call
# ./boost-hic.sh -b /mnt/imaging.data/mdas/combine_N2_Arima_hicpro/hic_results/matrix/N2/raw/5000/N2_5000_abs.bed -m /mnt/imaging.data/mdas/combine_N2_Arima_hicpro/hic_results/matrix/N2/raw/5000/N2_5000.matrix -c chrII -o ./results/N2_chrII_ -f matrix boost

## to convert to .cool
# cd results
# singularity exec ../../hicpro_3.0.0_ubuntu.img $HICPRO_PATH/hicpro2higlass.sh -i N2_chrI_boostedmat.matrix -r 10000 -c chromosome_sizes.tsv -n
