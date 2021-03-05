#! /usr/bin/env python

import numpy as np
import pandas as pd
import math
import cooler
import h5py
from matplotlib import pyplot as plt

import sys
import os
import re
import logging
import argparse

PUBLISHED_FOLDER = './data'

# Initialization
logger = logging.getLogger('hic_converters')

# used to rename pseudo chromosome so that they are alphabetic ordered
# allowing the produced .bedgraph files to be converted into .bw (bigwig)
CHR_TO_INT_DICT = {'chrI': '1', 'chrII': '2', 'chrIII': '3', 'chrIV': '4', 'chrV': '5', 'chrX': 'X'}


def diff_matrix_to_cool(matrix1_file, matrix2_file, chrom, resolution):
    # bin_w = resolution
    if matrix1_file.endswith('gz'):
        matrix1_pd = pd.read_table(matrix1_file, delimiter='\t', index_col=0, comment='#', compression='gzip')
    else:
        matrix1_pd = pd.read_table(matrix1_file, delimiter='\t', index_col=0, comment='#')
    hic1 = matrix1_pd.to_numpy()
    np.fill_diagonal(hic1, 0)

    if matrix2_file.endswith('gz'):
        matrix2_pd = pd.read_table(matrix2_file, delimiter='\t', index_col=0, comment='#', compression='gzip')
    else:
        matrix2_pd = pd.read_table(matrix2_file, delimiter='\t', index_col=0, comment='#')
    hic2 = matrix2_pd.to_numpy()
    np.fill_diagonal(hic2, 0)

    hic_diff = hic2 - hic1
    cool_file = matrix1_file + f'_diff.{resolution}.cool'

    return hic_to_cool(hic_diff, chrom, resolution, cool_file)


def matrix_to_cool(matrix_file, chrom, resolution):
    if matrix_file.endswith('gz'):
        matrix_pd = pd.read_table(matrix_file, delimiter='\t', index_col=0, comment='#', compression='gzip')
    else:
        matrix_pd = pd.read_table(matrix_file, delimiter='\t', index_col=0, comment='#')

    logger.info(f"data.shape: ${matrix_pd.shape}")

    hic = matrix_pd.to_numpy()
    np.fill_diagonal(hic, 0)

    cool_file = matrix_file + f'.{resolution}.cool'
    return hic_to_cool(hic, chrom, resolution, cool_file)


def load_hic_hdf5(hic_file):
    with h5py.File(hic_file, 'r') as f:
        # List all groups
        logger.info(f"Load hic file {hic_file}.")
        a_group_key = list(f.keys())[0]

        logger.info("a_group_key0: %s" % a_group_key)

        # Get the data
        data = list(f[a_group_key])
        hic = np.array(data)
        logger.info(f"hic.shape: ${hic.shape}")

        return hic


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
                logger.info(f'pixels progress: {progress}%')
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


def cool_chrom_to_matrix(cool_file, chrom=None, balanced=False, res=10000, matrix_file=None):
    """
    Extract chromosome matrix in TSV file format from a cooler file.
    :param cool_file: the cooler file to extract from. It could be .cool or .mcool file
    :param chrom: the chromosome to extract. Default: the all chromosomes in the cooler file
    :param balanced: to extract balanced date, otherwise not
    :param res: resolution if given a .mcool file
    :param matrix_file: the output matrix file. Default name: f'{cool_file}_{chrom}_res{res}.matrix.tsv'
    :return: the name of the
    """
    mat_cooler = get_cooler(cool_file, res=res)
    if chrom:
        logger.info(f'fetching chromosome {chrom} from {cool_file} ...')
    else:
        # default
        logger.info(f'fetching all chromosomes from {cool_file} ...')
    if not matrix_file:
        # default
        matrix_file = f'{cool_file}{f"_{chrom}" if chrom else ""}_res{res}{"b" if balanced else ""}.matrix.tsv'
    # name = os.path.splitext(os.path.basename(cool_file))[0]
    # species = "ce11"
    if chrom:
        mat = mat_cooler.matrix(balance=balanced).fetch(chrom)
    else:
        mat = mat_cooler.matrix(balance=balanced, sparse=False)[:, :]  # all chromosomes
    # res = mat_cooler.binsize
    # bin_names = [f'{name}|{species}|{chrom}:{bi*res}-{(bi+1)*res}' for bi in range(mat.shape[0])]
    mat_df = pd.DataFrame(data=mat)   # TODO add it optional: , index=bin_names, columns=bin_names)
    mat_df.to_csv(matrix_file, sep='\t', header=False, index=False)
    logger.info(f'Saved cool to matrix file {matrix_file}')
    return matrix_file


def matrix_to_mcool(matrix_file, chrom, resolution, factors):
    cool_file = matrix_file + f'.{resolution}.cool'
    if not os.path.isfile(cool_file):
        matrix_to_cool(matrix_file, chrom, resolution)  # == cool_file
    resolutions = [int(i * resolution) for i in factors]
    mcool_file = matrix_file + f'.{resolutions[0]}.mcool'
    cooler.zoomify_cooler(cool_file, mcool_file, resolutions=resolutions, chunksize=int(10e6))
    return mcool_file


def get_cooler(cool_file, res=10000):
    if cool_file.endswith('.mcool'):
        hic_cooler = cooler.Cooler(f'{cool_file}::/resolutions/{res}')
    else:
        hic_cooler = cooler.Cooler(f'{cool_file}::/')
    return hic_cooler


def hdf5_to_chrom_cool(hic_file, chrom, resolution=2000, cool_file=None):
    if not cool_file:
        # default
        cool_file = f'{hic_file}.cool'

    hic = load_hic_hdf5(hic_file)
    if not os.path.isfile(cool_file):
        cool_file = hic_to_cool(hic, chrom=chrom, resolution=resolution, cool_file=cool_file)
    return cool_file


def cool_to_mcool(cool_file, mcool_file=None, resolutions=[10000, 20000, 50000, 100000]):
    if not mcool_file:
        mcool_file = re.sub(r'.cool', f'.mcool', cool_file)
    if not os.path.isfile(mcool_file):
        hic_cooler = get_cooler(cool_file)
        resolutions.insert(0, hic_cooler.binsize)
        cooler.zoomify_cooler(cool_file, mcool_file, resolutions=resolutions,
                              chunksize=int(10e6))
    return mcool_file


def hic_prob_model(cool_file, chrom, balanced=False, res=10000):
    hic_cooler = get_cooler(cool_file, res=res)
    logger.info(f'hic_prob_model({cool_file})')
    # mat_pixels = hic_cooler.matrix(balance=balanced, as_pixels=True).fetch(chrom)
    mat_pixels = hic_cooler.matrix(balance=balanced, sparse=True).fetch(chrom)
    logger.info(f'{mat_pixels}')
    chrom_size = hic_cooler.chromsizes[chrom]
    alpha = 0
    for i, j, c in zip(mat_pixels.row, mat_pixels.col, mat_pixels.data):
        # logger.info(f'#{mat_pixels[i]}')
        if i != j:
            alpha += c*math.log(abs(i-j))
        # logger.info(f'{i}-{j}:#{c} \t alpha={alpha}')
    logger.info(f'alpha={alpha}')
    dp = []
    di = []
    for i in range(chrom_size):
        di.append(i)
        dp.append(math.pow(i, alpha))

    plt.plot(di, dp)

    # for (i, j, counts) in mat_pixels:
    return


def extract_AB_compartment(cool_file, compAB_bedgraph, res=5000):

    hic_cooler = get_cooler(cool_file, res=res)
    chrom_sizes = hic_cooler.chromsizes
    ab_pd = pd.read_csv(compAB_bedgraph, delim_whitespace=True, names=['chrom', 'start', 'end', 'value'],
                        encoding='utf-8',
                        dtype={"start": "int64", "end": "int64", "value": "float64"})

    resolution = int(hic_cooler.binsize)
    bins = pd.DataFrame(columns=['chrom', 'start', 'end', 'weight'])
    pixels = pd.DataFrame(columns=['bin1_id', 'bin2_id', 'count'])

    # find Bleft/right, A
    idx_offs = 0
    for chrom in hic_cooler.chromnames:
        ab_chrom_pd = ab_pd[ab_pd['chrom'] == chrom]
        chrom_size = chrom_sizes[chrom]
        chrom_mid = chrom_size//2

        chrom_b_left = 0
        # from B to A
        # for index, row in ab_chrom_pd.iterrows():
        #    if row['value'] < 0:
        #        chrom_b_left = row['start']
        #        break
        # from A back to B
        ab_chrom_left_pd = ab_chrom_pd[ab_chrom_pd['start'] <= chrom_mid]
        for index, row in ab_chrom_left_pd.iloc[::-1].iterrows():  # reversed
            if row['value'] > 0:
                chrom_b_left = row['end']
                break

        # find the B-right compartment strating from the middle
        ab_chrom_right_pd = ab_chrom_pd[ab_chrom_pd['start'] > chrom_mid]
        chrom_b_right = 0
        for index, row in ab_chrom_right_pd.iterrows():
            if row['value'] > 0:
                chrom_b_right = row['start']
                break

        # logger.info(f'{chrom}: compBleft={chrom_b_left}, chrom_b_right={chrom_b_right}')

        # crop compBleft
        chrom_b_left_pixels = hic_cooler.matrix(balance=False, as_pixels=True).fetch((chrom, 0, chrom_b_left))
        chrom_b_left_bins = hic_cooler.bins().fetch((chrom, 0, chrom_b_left))  # (f'{chrom}:0-{chrom_b_left}')
        # rename
        chrom_b_left_bins['chrom'] = chrom_b_left_bins['chrom'].replace([chrom], f'{CHR_TO_INT_DICT[chrom]}_1B')
        # shift
        chrom_b_left_bins, chrom_b_left_pixels, idx_offs = shift(chrom_b_left_bins, chrom_b_left_pixels, idx_offs)

        # crop compA between Bleft and Bright
        chrom_a_pixels = hic_cooler.matrix(balance=False, as_pixels=True).fetch((chrom, chrom_b_left, chrom_b_right))
        chrom_a_bins = hic_cooler.bins().fetch((chrom, chrom_b_left, chrom_b_right))
        # rename
        chrom_a_bins['chrom'] = chrom_a_bins['chrom'].replace([chrom], f"{CHR_TO_INT_DICT[chrom]}_2A")
        # shift
        chrom_a_bins, chrom_a_pixels, idx_offs = shift(chrom_a_bins, chrom_a_pixels, idx_offs)

        # crop compBright
        chrom_b_righ_pixels = hic_cooler.matrix(balance=False, as_pixels=True).fetch((chrom, chrom_b_right, chrom_size))
        chrom_b_righ_bins = hic_cooler.bins().fetch((chrom, chrom_b_right, chrom_size))  # (f'{chrom}:{chrom_b_right}-')
        # rename
        chrom_b_righ_bins['chrom'] = chrom_b_righ_bins['chrom'].replace([chrom], f'{CHR_TO_INT_DICT[chrom]}_3B')
        # shift
        chrom_b_righ_bins, chrom_b_righ_pixels, idx_offs = shift(chrom_b_righ_bins, chrom_b_righ_pixels, idx_offs)

        bins = pd.concat([bins, chrom_b_left_bins, chrom_a_bins, chrom_b_righ_bins])
        pixels = pd.concat([pixels, chrom_b_left_pixels, chrom_a_pixels, chrom_b_righ_pixels])

        logger.info(f'{chrom}: compBleft={chrom_b_left}, chrom_b_right={chrom_b_right}')

    compAB_file_name = f'{cool_file.replace(".cool", "")}_compAB'
    compAB_cool_file = f'{compAB_file_name}.cool'
    metadata = {'format': 'HDF5::Cooler',
                'format-version': '0.8.10',
                'bin-type': 'fixed',
                'bin-size': resolution,
                'storage-mode': 'symmetric-upper',
                'genome-assembly': 'ce11',
                'generated-by': 'boost-hic',
                # 'creation-date': datetime.date.today()
                }

    # count_dtypes = {'count': 'float64'}
    count_dtypes = {'count': 'int'}
    # hic_bins = hic_cooler.bins()[:]
    # hic_bins = hic_bins.astype({'chrom': str, 'start': int, 'end': int, 'weight': float})
    # remove wight as we'll calculate it later
    bins = bins[['chrom', 'start', 'end']]
    bins = bins.astype({'chrom': str, 'start': int, 'end': int})   # , 'weight': float})
    pixels = pixels.astype({'bin1_id': int, 'bin2_id': int, 'count': int})
    cooler.create_cooler(compAB_cool_file, bins=bins, pixels=pixels, dtypes=count_dtypes,
                         ordered=True, ensure_sorted=True, metadata=metadata)

    # save chromosome_sizes
    hicBlr_cooler = cooler.Cooler(f'{compAB_cool_file}::/')
    chrom_names = hicBlr_cooler.chromnames
    chrom_sizes = hicBlr_cooler.chromsizes
    chrom_sizes_df = pd.DataFrame({"chromnames": chrom_names, "chromsizes": chrom_sizes})
    chrom_sizes_df.to_csv(f"{compAB_file_name}_chrom_sizes.tsv", sep="\t", index=False, header=0)

    resolutions = [5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000]

    cmd = f'cooler balance --cis-only --force {compAB_cool_file}'
    logger.info(f'CALL: {cmd}')
    os.system(cmd)

    # create .mcool file to ingest in HiGlass, could possible resolutions abbreviation: 5000N
    resolutions_str = ','.join([str(r) for r in resolutions] )
    cmd = f'cooler zoomify -r "{resolutions_str}" {compAB_cool_file}'
    logger.info(f'CALL: {cmd}')
    os.system(cmd)

    # compBlr_mcool_file = f'{compAB_file_name}.mcool'
    # cooler.zoomify_cooler(compAB_cool_file, compBlr_mcool_file, resolutions=resolutions, chunksize=int(10e6))

    return compAB_cool_file


def shift(chrom_bins, chrom_pixels, idx_offs):
    if idx_offs > 0:
        idx_0 = chrom_bins.index[0]
        idx_skip = idx_0 - idx_offs
        if idx_skip > 0:
            chrom_bins.index -= idx_skip
            chrom_pixels = chrom_pixels.apply(lambda v: v - idx_skip if v.name == 'bin1_id' else v)
            chrom_pixels = chrom_pixels.apply(lambda v: v - idx_skip if v.name == 'bin2_id' else v)

    # shift bins starting from 0
    start_0 = chrom_bins.iloc[0]['start']
    if start_0 > 0:
        chrom_bins = chrom_bins.apply(lambda v: v - start_0 if (v.name == 'start') or (v.name == 'end') else v)

    idx_offs += len(chrom_bins.index)
    return chrom_bins, chrom_pixels, idx_offs


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    logging.getLogger("").setLevel(logging.INFO)

    p = argparse.ArgumentParser()
    p.add_argument("-cf", "--cool_file", default=f"{PUBLISHED_FOLDER}/wt_N2_Brejc2017_5000.cool",
                   help="Experimental cool file to plot")
    p.add_argument("-o", "--output_folder", default=".", help="output folder")
    p.add_argument("-b", "--balanced", action='store_true', help="Balanced or not")
    p.add_argument("-z", "--z_score", action='store_true', help="Z-score normalized otherwise original Hi-C")
    p.add_argument("-chrs", "--chr_synonyms", nargs="+", default=['chrX', 'X', '6'],
                   help="List of chromosome synonyms to look for")
    p.add_argument("-ab", "--ab_compartment",
                   help="A/B compartment profile in bredgraph file format.")
    p.add_argument("-c", "--cmap", default="hot_r",
                   help="Color map: cool, hot_r, gist_heat_r, afmhot_r, YlOrRd, Greys, gist_yarg, seismic")
    p.add_argument("-f", "--file_extension", default="png", help="File format extension: png, tif")
    args = p.parse_args(sys.argv[1:])
