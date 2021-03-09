import argparse

import numpy as np
import pandas as pd
import sys
from matplotlib import pyplot as plt

p = argparse.ArgumentParser()
p.add_argument("filtered_bins", help="Filepath to a *_filteredbin.txt file")  # 'results/N2_5000_chrI_filteredbin.txt'
args = p.parse_args(sys.argv[1:])

filtered_bins = pd.read_csv(args.filtered_bins, names=['bin'])

filtered_bins = filtered_bins[['bin']]['bin'].astype(np.int64)
fig = plt.figure()
plt.title(f'Filtered bins #{np.sum(filtered_bins)} for {args.filtered_bins}')
plt.plot(filtered_bins, color='black', linewidth=1)
plt.xlabel('bin position')
plt.ylabel('filtered bins')

plt.show()
plt.close()
