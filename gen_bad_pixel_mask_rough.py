import os

import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt

cutoff_percentile = 30
directory = r'G:\RIMAS_data\reduced\20250906'
prefix = '20250906.rimas.0009'
suffixes = ['.YJ.fits', '.HK.fits']
out_fnames = ['bad_pixels.imaging.YJ.fits', 'bad_pixels.imaging.HK.fits']
files = [os.path.join(directory, prefix+suffix) for suffix in suffixes]

for f, out_f in zip(files, out_fnames):
    print(f)
    data = fits.getdata(f)
    cutoff = np.nanpercentile(data, cutoff_percentile)
    print(cutoff)
    bad_pixels = data < cutoff
    bad_pixel_map = np.zeros(data.shape, dtype=np.uint16)
    bad_pixel_map[bad_pixels] = 1
    for axis in (0,1):
        for roll in (-1, 1):
            bad_pixel_map += np.roll(bad_pixel_map, roll, axis=axis)
    bad_pixel_map_with_neighbors = np.zeros(data.shape, dtype=np.uint16)
    bad_pixels = bad_pixel_map != 0
    bad_pixel_map_with_neighbors[bad_pixels] = 1
    plt.imshow(bad_pixel_map_with_neighbors)
    plt.show()
    fits.HDUList([fits.PrimaryHDU(data=bad_pixel_map_with_neighbors)]).writeto(out_f, overwrite=True)
