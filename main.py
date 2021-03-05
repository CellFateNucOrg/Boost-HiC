import h5py
import numpy as np
import logging
import argparse
import sys

#my own toolkit
import HiCutils
import utils
import convert

logging.basicConfig(level=logging.DEBUG)
logging.getLogger("").setLevel(logging.INFO)
logger = logging.getLogger(f'Boos-HiC')

p = argparse.ArgumentParser()
p.add_argument("-b", "--bedfilename", required=True, help="bed file of genomic coordinate of each bin")
p.add_argument("-m", "--matrixfilename", required=True,
			   help="contact map stored in tab separated file as : "
					"bin_i / bin_j / counts_ij Only no zero values are stored. Contact map are symmetric")
p.add_argument("-o", "--output", default="./results/", help="output foolder where files are stored")
p.add_argument("--chr", default="all", help="Which chromosome or the whole genome to boost.")
#p.add_argument("-r", "--resolution", default=10000, help="Matrix Resolution")
p.add_argument("operation", default="boost", choices=["boost", "sample"],
			   help="Operation to be executed")
args = p.parse_args(sys.argv[1:])

### YOU ARE SUPPOSED TO ONLY MODIFY VALUE HERE ###
#input file
bedfilename = args.bedfilename  # '/mnt/imaging.data/mdas/combine_N2_Arima_hicpro/hic_results/matrix/N2/raw/5000/N2_5000_abs.bed'
# '/Users/todor/unibe/data/combine_N2_Arima_hicpro/N2_5000_abs.bed'
matrixfilename = args.matrixfilename   # '/mnt/imaging.data/mdas/combine_N2_Arima_hicpro/hic_results/matrix/N2/raw/5000/N2_5000.matrix'
# '/Users/todor/unibe/data/combine_N2_Arima_hicpro/N2_5000.matrix'
repositoryout = args.output   # './results/'
achr = args.chr   # "genome"
Operation = args.operation     # 'Boost'

#default parameter
#resolution=10000 #default : 10kb
alpha=0.2 #AFTER a lot of test : 0.24 is always a good and safe compromise, you must use this value
###


def BoostHiC(amat):
	normmat=HiCutils.SCN(np.copy(amat))
	FFmat=np.power(HiCutils.fastFloyd(1/np.power(normmat.copy(),alpha)),-1/alpha) #to dist, FF, to contact in one line
	boostedmat=HiCutils.adjustPdS(normmat,FFmat)
	return boostedmat

def Sample(amat,repositoryout):
	percentofsample=[0.1,1.,10.]
	for j in percentofsample:
		logger.info(f"Value of sample: {j}")
		chrmat_s=np.copy(amat)
		chrmat=HiCutils.downsample_basic(chrmat_s,j)
		fh5 = h5py.File(repositoryout+"inputmat_sampleat_"+str(j)+"_percent.hdf5", "w")
		fh5['data'] = chrmat
		fh5.close()



### CODE EXECUTION ###

# load the data
logger.info("LOADING MATRIX")
D=convert.loadabsdatafile(bedfilename)
print(*D.items(), sep='\n')
beginfend=D[achr][0]
endfend=D[achr][1]
logger.info(f"Data fend : {beginfend},{endfend}")
basemat=convert.loadmatrixselected(matrixfilename,beginfend,endfend)

#matrix filtering
logger.info("FILTERING")
pos_out=HiCutils.get_outliers(basemat)
basematfilter=basemat[np.ix_(~pos_out, ~pos_out)]
basematfilter=np.copy(basematfilter)
#basematfilter=basematfilter[0:1000,0:1000]
logger.info(f'len(basemat):{len(basemat)}, len(basematfilter):{len(basematfilter)}')
fh5 = h5py.File(repositoryout+"inputmat.hdf5", "w")
fh5['data'] = basemat
fh5.close()
fh5 = h5py.File(repositoryout+"inputmat_filtered.hdf5", "w")
fh5['data']=basematfilter
fh5.close()
utils.savematrixasfilelist3(pos_out,repositoryout+"filteredbin.txt")

if Operation=="boost":
	logger.info("Boost Hic")
	boosted=BoostHiC(basematfilter)
	#save
	fh5 = h5py.File(repositoryout+"boostedmat.hdf5", "w")
	fh5['data']=boosted
	fh5.close()
elif Operation=="sample":
	logger.info("SAMPLING")
	Sample(basematfilter,repositoryout)




