#!/usr/bin/python3 -u

import argparse
import logging
import os

import h5py
import numpy as np
import pandas as pd
import sys

# my own toolkit
import HiCutils
import convert
import utils

logging.basicConfig(level=logging.DEBUG)
logging.getLogger("").setLevel(logging.INFO)
logger = logging.getLogger(f'Boos-HiC')

p = argparse.ArgumentParser()
p.add_argument("operation", default="boost", choices=["boost", "sample"],
               help="Operation to be executed")
p.add_argument("-b", "--bedfilename", required=True, help="bed file of genomic coordinate of each bin")
p.add_argument("-m", "--matrixfilename", required=True,
               help="contact map stored in tab separated file as : "
                    "bin_i / bin_j / counts_ij Only no zero values are stored. Contact map are symmetric")
# TODO add support for input .cool file
p.add_argument("-c", "--chromosomes", nargs='+', help="Which chromosomes to boost, otherwise all chromosomes")
p.add_argument("-o", "--output_prefix", default="./boosted/", help="prefix for output files")
p.add_argument("-f", "--format", default="cool", choices=["cool", "hdf5"], help="output file format")
p.add_argument("-g", "--genome_assembly", default="ce11", help="genome assembly as metadata for .cool file")
p.add_argument("-k", "--keep_filtered_bins", action='store_true',
               help="Whether to keep filtered out bins, otherwise they will be removed from the result matrix. "
                    "Not used yet.")
p.add_argument("-a", "--alpha", default=0.2,
               help="AFTER a lot of test : 0.24 is always a good and safe compromise, you must use this value")
args = p.parse_args(sys.argv[1:])

# input file
Operation = args.operation
bedfilename = args.bedfilename
matrixfilename = args.matrixfilename
chromosomes = args.chromosomes
format = args.format
keep_filtered_bins = args.keep_filtered_bins
genome_assembly = args.genome_assembly
alpha = args.alpha

###


def BoostHiC(amat):
    normmat = HiCutils.SCN(np.copy(amat))
    FFmat = np.power(HiCutils.fastFloyd(1 / np.power(normmat.copy(), alpha)),
                     -1 / alpha)  # to dist, FF, to contact in one line
    boostedmat = HiCutils.adjustPdS(normmat, FFmat)
    return boostedmat


def Sample(amat, repositoryout):
    percentofsample = [0.1, 1., 10.]
    for j in percentofsample:
        logger.info(f"Value of sample: {j}")
        chrmat_s = np.copy(amat)
        chrmat = HiCutils.downsample_basic(chrmat_s, j)
        fh5 = h5py.File(repositoryout + "inputmat_sampleat_" + str(j) + "_percent.hdf5", "w")
        fh5['data'] = chrmat
        fh5.close()


# ## CODE EXECUTION ## #
# load the data
logger.info("LOADING MATRIX")
D, total, resolution = convert.loadabsdatafile(bedfilename)
print(*D.items(), sep='\n')
print(f'Total bins:{total} resolution:{resolution}')

bins_boosted = pd.DataFrame(columns=['chrom', 'start', 'end'])
pixels_boosted = pd.DataFrame(columns=['bin1_id', 'bin2_id', 'count'])
chroms = chromosomes if chromosomes else D.keys()
bin_offs = 0
for chrom in chroms:
    repositoryout = f'{args.output_prefix}_{chrom}_'

    beginfend = D[chrom][0]
    endfend = D[chrom][1]
    logger.info(f"Chromosome {chrom} data fend : {beginfend},{endfend}")
    basemat = convert.loadmatrixselected(matrixfilename, beginfend, endfend)

    # matrix filtering
    logger.info("FILTERING")
    bins_num = basemat.shape[0]
    pos_out = HiCutils.get_outliers(basemat)
    utils.savematrixasfilelist3(pos_out, repositoryout + "filteredbin.txt")
    basematfilter = basemat[np.ix_(~pos_out, ~pos_out)]
    basematfilter = np.copy(basematfilter)
    # basematfilter=basematfilter[0:1000,0:1000]
    logger.info(f'len(basemat):{len(basemat)}, len(basematfilter):{len(basematfilter)}')
    if format is None or format == "hdf5":
        fh5 = h5py.File(repositoryout + "inputmat.hdf5", "w")
        fh5['data'] = basemat
        fh5.close()
    if format is None or format == "cool":
        convert.hic_to_cool(basemat, chrom, resolution, repositoryout + "inputmat.cool",
                            genome_assembly=genome_assembly)
    if format is None or format == "hdf5":
        fh5 = h5py.File(repositoryout + "inputmat_filtered.hdf5", "w")
        fh5['data'] = basematfilter
        fh5.close()
    if format is None or format == "cool":
        convert.hic_to_cool(basematfilter, chrom, resolution, repositoryout + "inputmat_filtered.cool",
                            genome_assembly=genome_assembly)

    if Operation == "boost":
        logger.info("Boost Hic")
        boosted = BoostHiC(basematfilter)
        # save
        if format is None or format == "hdf5":
            fh5 = h5py.File(repositoryout + "boostedmat.hdf5", "w")
            fh5['data'] = boosted
            fh5.close()
        if format is None or format == "cool":
            filtered_bins = pos_out if keep_filtered_bins else None
            chrom_bins, chrom_pixels = convert.get_bins_pixels(boosted, chrom, resolution,
                                                               bin_offs=bin_offs, bins_num=bins_num,
                                                               filtered_bins=filtered_bins)
            bins_boosted = pd.concat([bins_boosted, chrom_bins])
            pixels_boosted = pd.concat([pixels_boosted, chrom_pixels])
            bin_offs += bins_num

    elif Operation == "sample":
        logger.info("SAMPLING")
        Sample(basematfilter, repositoryout)

if Operation == "boost" and format is None or format == "cool":  # combined file support only for .cool
    repositoryout = args.output_prefix + (f'_{"_".join(chromosomes)}_' if chromosomes else '_')
    cool_file = f"{repositoryout}boosted{'_kfb' if keep_filtered_bins else ''}.cool"
    convert.create_cool(bins_boosted, pixels_boosted, resolution, cool_file, genome_assembly=genome_assembly)

    cmd = f'cooler balance --cis-only --force {cool_file}'
    logger.info(f'CALL: {cmd}')
    os.system(cmd)

    resolutions = [5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000]
    resolutions_str = ','.join([str(r) for r in resolutions])
    cmd = f'cooler zoomify -r "{resolutions_str}" {cool_file}'
    logger.info(f'CALL: {cmd}')
    os.system(cmd)