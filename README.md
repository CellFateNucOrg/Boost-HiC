# Boost-HiC

Software requiered
=================
Boost-HiC current implementation is in python 3 and need the current list of package :
-h5py
-numpy
-copy
-sklearn
-scipy
-skimage

Installation :
=================

clone repository and install conda environment with the following commands:

    git clone https://github.com/CellFateNucOrg/Boost-HiC.git
    cd Boost-HiC
    ./install.sh

Troubleshooting:

If you encounter some problems running install.sh script, you can execute the following commands manually:

    conda create -n boost-hic python=3.8
    conda activate boost-hic
    pip install -e .

You can activate conda environment like this:

    conda activate boost-hic


You can also add the following line to your .bashrc, or .zshrc, or .bash_profile:

    export PATH=$PATH:<BOOST_INSTALLATION_FOLDER>"

You can check your active shell like this:

    echo $SHELL


You can also change and update your start shell script boost-hic.sh command.


You can complete the installation with the proposed changes at the end.
Help can be recalled with the command:
    
    # if you are in the instalation folder:
    ./boost-hic.py --help
    # or if you properly updated your start shell script meant for running also in a SLURM environment:
    boost-hic.sh --help

Input :
=================
Boost-HiC use HiC-Pro output format (described in : http://nservant.github.io/HiC-Pro/MANUAL.html#browsing-the-results ) for raw contact map.

The contact map is stored in tab separated file as :
bin_i / bin_j / counts_ij
Only no zero values are stored. Contact map are symmetric

The bin are described in a separated bed file which give the genomic coordinate of each bin.

In a first step, the contact map are convert in hdf5 by the pipeline.

How to use it
=================
There is a convenient command line interface (CLI) provided by the python script boost-hic.py 

Important parameter are :

-bedfilename : bed file of genomic coordinate of each bin

-Operation : 'Boost' or 'Sample'

-output_prefix : prefix used to save the output files

-format : output format. Possible values: cool, hdf5. Default: cool

-alpha : the alpha used to Boost-HiC, with preset default value

Some usefull tool are available in HiCtools.py if you need to made your own script.

Here all CLI (command line interface) parameters:


    usage: boost-hic.py [-h] -b BEDFILENAME -m MATRIXFILENAME [-c CHROMOSOMES [CHROMOSOMES ...]] [-o OUTPUT_PREFIX]
    [-f {cool,hdf5}] [-g GENOME_ASSEMBLY] [-k]
    {boost,sample}
    
    positional arguments:
    {boost,sample}        Operation to be executed
    
    optional arguments:
    
    -h, --help            show this help message and exit
    
    -b BEDFILENAME, --bedfilename BEDFILENAME
    bed file of genomic coordinate of each bin
    
    -m MATRIXFILENAME, --matrixfilename MATRIXFILENAME
    contact map stored in tab separated file as : bin_i / bin_j / counts_ij Only no zero
    values are stored. Contact map are symmetric
    
    -c CHROMOSOMES [CHROMOSOMES ...], --chromosomes CHROMOSOMES [CHROMOSOMES ...]
    Which chromosomes to boost, otherwise all chromosomes
    
    -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
    Prefix for output files, including output folder. If not given, it will be the same
    as the input matrixfilename without its file extension butplus '\_boosted\_'.
    
    -f {cool,hdf5}, --format {cool,hdf5}
    output file format
    
    -g GENOME_ASSEMBLY, --genome_assembly GENOME_ASSEMBLY
    genome assembly as metadata for .cool file
    
    -k, --keep_filtered_bins
    Whether to keep filtered out bins, otherwise they will be removed from the result
    matrix. Not used yet.
    
    -a ALPHA, --alpha ALPHA
    AFTER a lot of test : 0.24 is always a good and safe compromise, you must use this
    value

Example call:

    boost-hic.py boost -b ./N2_5000_abs.bed -m N2_5000.matrix -o results/N2_5000 -f cool -k

Or if you are in SLURM environment and adapted your boost-hic.sh script

    sbatch boost-hic.sh boost -b ./N2_5000_abs.bed -m N2_5000.matrix -o results/N2_5000 -f cool -k

Output :
=================
Every output contact map are stored in hdf5 with the key 'data' in the hdf5 dict.
Alternatively, you can have the results in a cooler format by using the '--format cool' option.

-inputmat.hdf5/cool  : The contact map that is load in .matrix at the begin, just convert in hdf5

-inputmat_filtered.hdf5/cool : The original contact map with some filtered bin

-filteredbin.txt : list of raw/col that are filtered or not in the contact map as a boolean list.

-boostedmat.hdf5/cool : The input mat improved by BoostHiC procedure.


