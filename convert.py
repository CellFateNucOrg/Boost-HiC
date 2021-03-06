#Leopold Carron
#december 2016
#CC-By-SA
#now in py 3


import numpy as np
import pandas as pd
import cooler


def loadabsdatafile(filein):
	"""
	in a _abs.bedfile
	out : a dict with chr=[begin,end] in the bed, as finaly in the matrix
	"""
	resolution = 0
	fin=open(filein,"r")
	d={}
	d["Total"]=0
	l=fin.readline()
	while l:
		ls=l.split()
		if d.__contains__(ls[0]):
			d[ls[0]][1]=float(ls[3])
		else: #init
			d[ls[0]]=[float(ls[3]),float(ls[3])]
		if resolution == 0:
			resolution = int(ls[2]) - int(ls[1])
		d["Total"]+=1
		l=fin.readline()
	fin.close()
	d["all"]=[float(1), float(d["Total"])]
	return d, resolution

def loadmatrix(filein,sizemat):
	"""
	in : a matrix file, is size
	out : the matrix generated
	"""
	print(sizemat)
	mat=np.zeros((sizemat,sizemat))
	fin=open(filein,"r")
	l=fin.readline()
	while l:
		ls=l.split()
		i=int(float(ls[0])-1)
		j=int(float(ls[1])-1)
		v=float(ls[2])
		mat[i,j]=v
		mat[j,i]=v
		l=fin.readline()
	fin.close()
	print("Number of contact in map :",np.sum(np.triu(mat)))
	return mat

def loadmatrixselected(filein,B,E):
	"""
	-in : a matrix file, is size
	-out : the matrix generated
	-i-B>=0 j-E>=0
	"""
	#E-=1
	B-=1 #file start at one
	sizemat=int(E-B)
	print("Matrix size :",sizemat)
	mat=np.zeros((sizemat,sizemat))
	fin=open(filein,"r")
	l=fin.readline()
	while l:
		ls=l.split()
		i=int(float(ls[0])-B)
		j=int(float(ls[1])-B)
		if i>=0 and j>=0 and i<sizemat and j<sizemat:
			#print(ls[0],ls[1],i,j,B,E)
			v=float(ls[2])
			mat[i,j]=v
			mat[j,i]=v
		l=fin.readline()
	fin.close()
	print("Number of contact in the map :",np.sum(np.triu(mat)))
	return mat

def hic_to_cool(hic, chrom, resolution, cool_file):

	# build the cooler fields
	N = hic.shape[0]

	bins_index = [[chrom, i * resolution, i * resolution + resolution] for i in range(N)]
	bins = pd.DataFrame(data=bins_index, columns=['chrom', 'start', 'end'])  # , dtype=np.dtype([('','','')]))

	pixels_bin1_id = []
	pixels_bin2_id = []
	pixels_count = []

	tot_iter = (N - 1) * N / 2
	it = 0
	for bin1_id in range(N - 1):
		for bin2_id in range(bin1_id + 1, N):
			it += 1
			progress = (it / tot_iter) * 100
			if (progress % 10) == 0:
				print(f'pixels progress: {progress}%')
			count = hic[bin1_id, bin2_id]
			if count != 0:
				# pixels_pd = pixels_pd.append({'bin1_id': np.int64(bin1_id), 'bin2_id': np.int64(bin2_id),
				# 'count': count}, ignore_index=True)
				pixels_bin1_id.append(np.int64(bin1_id))
				pixels_bin2_id.append(np.int64(bin2_id))
				pixels_count.append(count)

	pixels_dic = {'bin1_id': pixels_bin1_id, 'bin2_id': pixels_bin2_id, 'count': pixels_count}
	metadata = {'format': 'HDF5::Cooler',
				'format-version': '0.8.10',
				'bin-type': 'fixed',
				'bin-size': resolution,
				'storage-mode': 'symmetric-upper',
				'genome-assembly': 'ce11',
				'generated-by': 'boost-hic',
				# 'creation-date': datetime.date.today()
				}

	count_dtypes = {'count': 'float64'}
	cooler.create_cooler(cool_file, bins=bins, pixels=pixels_dic, dtypes=count_dtypes, ordered=True, metadata=metadata)
	return cool_file
	# problem with showing .cool file in higlass but with .mcool it works
