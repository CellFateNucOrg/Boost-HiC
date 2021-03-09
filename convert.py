# Leopold Carron
# december 2016
# CC-By-SA
# now in py 3


import numpy as np
import pandas as pd
import cooler


def loadabsdatafile(filein):
    """
	in a _abs.bedfile
	out : a dict with chr=[begin,end] in the bed, as finaly in the matrix
	"""
    resolution = 0
    fin = open(filein, "r")
    d = {}
    d_Total = 0
    l = fin.readline()
    while l:
        ls = l.split()
        if d.__contains__(ls[0]):
            d[ls[0]][1] = float(ls[3])
        else:  # init
            d[ls[0]] = [float(ls[3]), float(ls[3])]
        if resolution == 0:
            resolution = int(ls[2]) - int(ls[1])
        d_Total += 1
        l = fin.readline()
    fin.close()
    # d["all"]=[float(1), float(d_Total)]
    return d, d_Total, resolution


def loadmatrix(filein, sizemat):
    """
	in : a matrix file, is size
	out : the matrix generated
	"""
    print(sizemat)
    mat = np.zeros((sizemat, sizemat))
    fin = open(filein, "r")
    l = fin.readline()
    while l:
        ls = l.split()
        i = int(float(ls[0]) - 1)
        j = int(float(ls[1]) - 1)
        v = float(ls[2])
        mat[i, j] = v
        mat[j, i] = v
        l = fin.readline()
    fin.close()
    print("Number of contact in map :", np.sum(np.triu(mat)))
    return mat


def loadmatrixselected(filein, B, E):
    """
	-in : a matrix file, is size
	-out : the matrix generated
	-i-B>=0 j-E>=0
	"""
    # E-=1
    B -= 1  # file start at one
    sizemat = int(E - B)
    print("Matrix size :", sizemat)
    mat = np.zeros((sizemat, sizemat))
    fin = open(filein, "r")
    l = fin.readline()
    while l:
        ls = l.split()
        i = int(float(ls[0]) - B)
        j = int(float(ls[1]) - B)
        if i >= 0 and j >= 0 and i < sizemat and j < sizemat:
            # print(ls[0],ls[1],i,j,B,E)
            v = float(ls[2])
            mat[i, j] = v
            mat[j, i] = v
        l = fin.readline()
    fin.close()
    print("Number of contact in the map :", np.sum(np.triu(mat)))
    return mat


def correct_bin_id(bin_id: int, filtered_bins: np.array) -> int:
    corrected_id = bin_id
    if filtered_bins is not None:
        filtered_offset = np.sum(filtered_bins[:corrected_id])
        while (corrected_id - bin_id) != filtered_offset:
            corrected_id += filtered_offset - (corrected_id - bin_id)
            filtered_offset = np.sum(filtered_bins[:corrected_id])
    return corrected_id


def get_bins_pixels(hic, chrom, resolution, bin_offs=0, bins_num=None, filtered_bins=None):
    """
    Get bins and pixels as needed for cooler format
    :param hic: matrix
    :param chrom: chromosome name
    :param resolution: resolution
    :param bin_offs: accumulated bin offest from the beginning t
    :param bins_num: chromosome bins numbers, otherwise take the size of the matrix
    :return: bins, pixels as pd.DataFrames
    """
    N_bins = bins_num if bins_num else hic.shape[0]
    bins_index = [[chrom, i * resolution, i * resolution + resolution] for i in range(N_bins)]
    bins = pd.DataFrame(data=bins_index, columns=['chrom', 'start', 'end'])

    N = hic.shape[0]
    pixels_bin1_id = []
    pixels_bin2_id = []
    pixels_count = []

    tot_iter = (N - 1) * N / 2
    it = 0
    print(f'get bins and pixels for {chrom} with resolution {resolution} with matrix size N={N}, '
          f'number of bins={N_bins}, max pixels: {int(tot_iter)}')
    for bin1_id in range(N - 1):
        for bin2_id in range(bin1_id + 1, N):
            it += 1
            if (it % 1000000) == 0:
                progress = (it / tot_iter) * 100
                print(f'pixels progress: {int(progress)}%')
            count = hic[bin1_id, bin2_id]
            if count != 0:
                # sum filtered up to bin_id
                bin1_id_corrected = correct_bin_id(np.int64(bin1_id), filtered_bins)
                pixels_bin1_id.append(bin_offs + bin1_id_corrected)
                bin2_id_corrected = correct_bin_id(np.int64(bin2_id), filtered_bins)
                pixels_bin2_id.append(bin_offs + bin2_id_corrected)
                pixels_count.append(count)
    print(f'pixels progress: 100%')
    pixels = pd.DataFrame({'bin1_id': pixels_bin1_id, 'bin2_id': pixels_bin2_id, 'count': pixels_count},
                          columns=['bin1_id', 'bin2_id', 'count'])
    return bins, pixels


def create_cool(bins, pixels, resolution, cool_file, genome_assembly):
    metadata = {'format': 'HDF5::Cooler',
                'format-version': '0.8.10',
                'bin-type': 'fixed',
                'bin-size': resolution,
                'storage-mode': 'symmetric-upper',
                'genome-assembly': genome_assembly,
                'generated-by': 'boost-hic',
                # 'creation-date': datetime.date.today()
                }

    count_dtypes = {'count': 'float64'}
    bins = bins.astype({'chrom': str, 'start': int, 'end': int})
    pixels = pixels.astype({'bin1_id': int, 'bin2_id': int, 'count': float})
    cooler.create_cooler(cool_file, bins=bins, pixels=pixels, dtypes=count_dtypes, ordered=True, ensure_sorted=True,
                         metadata=metadata)
    # cooler.create_cooler(cool_file, bins=bins, pixels=pixels, dtypes=count_dtypes, ordered=True, metadata=metadata)
    return cool_file


def hic_to_cool(hic, chrom, resolution, cool_file, bins_num=None, genome_assembly=''):
    bins, pixels = get_bins_pixels(hic, chrom, resolution, bins_num=bins_num)
    create_cool(bins, pixels, resolution, cool_file, genome_assembly=genome_assembly)
