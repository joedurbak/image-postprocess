import os

import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt

cutoff_percentile = 99.3
directory = r'F:\Leach data\fowler\20230306'
prefix = '20230306.rimas.0002.cds_100x60_Leach'
suffixes = ['_001.fits', '2_001.fits']
out_fnames = ['bad_pixels.YJ.fits', 'bad_pixels.HK.fits']
files = [os.path.join(directory, prefix+suffix) for suffix in suffixes]

for f, out_f in zip(files, out_fnames):
    print(f)
    data = fits.getdata(f)
    cutoff = np.percentile(data, cutoff_percentile)
    print(cutoff)
    bad_pixels = data > cutoff
    bad_pixel_map = np.zeros(data.shape, dtype=np.uint16)
    bad_pixel_map[bad_pixels] = 1
    plt.imshow(bad_pixel_map)
    plt.show()
    fits.HDUList([fits.PrimaryHDU(data=bad_pixel_map)]).writeto(out_f, overwrite=True)
