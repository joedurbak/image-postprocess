import os

import numpy as np
from astropy.io import fits

date = '20250619'
data_dir = r'E:\asdetector-data\output\res\{}'.format(date)
frames = range(0, 106)
exposures = range(67, 74)
bands = ('YJ', 'HK')
temperatures = (95, 97)
file_format = '{}.rimas.{:04d}.{}.{:04d}.fits'
output_dir = r'E:\asdetector-data\output\dark-archive'
output_format = os.path.join(output_dir, '{}.{}.{}K.fits')

if not os.path.exists(output_dir):
    os.mkdir(output_dir)

for band, temperature in zip(bands, temperatures):
    print(band)
    average_frames = []
    for frame in frames:
        print(frame)
        arrays = [fits.getdata(os.path.join(data_dir, file_format.format(date, exposure, band, frame))) for exposure in exposures]
        average_frames.append(np.median(arrays, axis=0))
    output_name = output_format.format(date, band, temperature)
    fits.HDUList([fits.PrimaryHDU(np.asarray(average_frames))]).writeto(output_name)
    print(output_name)

