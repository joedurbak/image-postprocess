import os

import numpy as np
from astropy.io import fits

date = '20251009'
data_dir = r'E:\asdetector-data\output\raw\{}'.format(date)
#           "E:\asdetector-data\output\raw\20251009\20251009.rimas.0038.YJ.0000.fits"
#           'E:\asdetector-data\output\raw\202501009\202501009.rimas.0038.YJ.0000.fits'
frames = range(0, 43)
exposures = range(38, 43)
exposures = [38, 40, 41, 42]
bands = ('YJ', 'HK')
temperatures = (100, 100)
file_format = '{}.rimas.{:04d}.{}.{:04d}.fits'
output_dir = r'E:\asdetector-data\output\dark-archive'
output_format = os.path.join(output_dir, '{}.{}.{}K_v2.fits')

if not os.path.exists(output_dir):
    os.mkdir(output_dir)

for band, temperature in zip(bands, temperatures):
    print(band)
    average_frames = []
    for frame in frames:
        print(frame)
        arrays = np.asarray([fits.getdata(os.path.join(data_dir, file_format.format(date, exposure, band, frame))).astype(float) for exposure in exposures])
        arrays = arrays - arrays[0]
        average_frames.append(np.median(arrays, axis=0))
    output_name = output_format.format(date, band, temperature)
    fits.HDUList([fits.PrimaryHDU(np.asarray(average_frames))]).writeto(output_name)
    print(output_name)

