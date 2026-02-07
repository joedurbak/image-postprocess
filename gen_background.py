import os

from astropy.io import fits
import numpy as np

date = '20251009'
file_numbers = range(68, 77)
cameras = ('YJ', 'HK')
file_format = '{}.rimas.{:04d}.{}.fits'
# basedir = r'G:\RIMAS_data\reduced\{}'.format(date)
basedir = r'D:\reduced\{}'.format(date)
file_format = os.path.join(basedir, file_format)
bad_maps = {cam: fits.getdata('bad_pixels.imaging.{}.fits'.format(cam)) for cam in cameras}

for cam in cameras:
    data = np.asarray([fits.getdata(file_format.format(date, n, cam)) for n in file_numbers])
    background = np.nanmedian(data, axis=0)
    background[bad_maps[cam].astype(bool)] = np.nan
    med = np.nanmedian(background)
    std = np.nanstd(background)
    print(cam, med, std)

