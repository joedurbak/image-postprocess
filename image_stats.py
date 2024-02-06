import os
import re
from glob import glob

from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt

file_suffixes = ['Leach_001.fits', '_Leach2_001.fits']
hdu_numbers = (4, 2)
x_ranges = ((1698, 1762), (2046, 2124))
y_ranges = ((2029, 2093), (1867, 1948))
intensity_ranges = ((0, 40000), (0, 40000))
focus_files = []
for suffix in file_suffixes:
    focus_files.append(glob('F:/Leach data/processed/20220614--YJ-Y--HK-unknown+2100*'+suffix))

[lst.sort() for lst in focus_files]

focus_positions = []
for focus_file_list in focus_files:
    focus_position_strings = [re.findall(r'\d+-\d+', f)[0] for f in focus_file_list]
    focus_positions.append(focus_position_strings)
focus_positions = [[int(s.split('-')[0]) for s in focus_position_list] for focus_position_list in focus_positions]

focus_positions_arr = np.asarray(focus_positions)
argsort = np.argsort(focus_positions_arr[0])
focus_positions_arr = focus_positions_arr[:, argsort]
focus_files = np.asarray(focus_files)[:, argsort]
print(focus_files)
print(focus_positions_arr)
"""
table_data = []
for positions, files in zip(focus_positions, focus_files):
    table1_data = []
    for p, f in zip(positions, files):
        hdus = fits.open(f)
        print(p, f)
        fmt = p + "\t{index}\t{mean}\t{median}\t{std}"
        for i, hdu in enumerate(hdus[3:]):
            data = hdu.data
            flat_data = data.reshape((data.shape[0]*data.shape[1],))
            std = np.std(data)
            mean = np.mean(data)
            median = np.median(data)
            header = fmt.format(index=i, std=std, mean=mean, median=median)
            hist = plt.hist(flat_data, bins=1000)
            plt.title(header)
            plt.show()
"""

for files, positions, hdu, x_range, y_range, intensity_range in zip(focus_files, focus_positions_arr, hdu_numbers, x_ranges, y_ranges, intensity_ranges):
    _f, axarr = plt.subplots(3, 5)
    n_plots = len(focus_files[0])
    plot_indices = [divmod(i + 1, 5) for i in range(n_plots)]
    _f.suptitle(files[0])
    for f, p, subplot in zip(files, positions, plot_indices):
        hdu_data = fits.getdata(f, hdu)
        cropped_data = hdu_data[y_range[0]:y_range[1], x_range[0]:x_range[1]]
        axarr[subplot[0], subplot[1]].imshow(cropped_data, vmin=intensity_range[0], vmax=intensity_range[1])
        axarr[subplot[0], subplot[1]].set_title(p)
    plt.show()
