import os

import numpy as np
from astropy.io import fits
from pandas import read_csv

overwrite=True
obsdate = '20230802'
obslog_file_format = '{}-obslog.tsv'
# obsdate_file_abs = os.path.join('F:', 'leach-obslogs', obslog_file_format.format(obsdate))
obsdate_file_abs = os.path.join('E:', 'obslogs', obslog_file_format.format(obsdate))

# output_directory = os.path.join('F:', 'Leach data', 'obslog_reduce', obsdate)
output_directory = os.path.join('E:', 'asdetector-data', 'output', 'obslog_reduce', obsdate)

if not os.path.isdir(output_directory):
    os.makedirs(output_directory)

# fits_file_format = '{}.rimas.{:04d}.{}.fits'
fits_file_format = '{}.rimas.{:04d}.{}.ramp.fits'
# fits_file_format_abs = os.path.join('F:', 'Leach data', 'fowler', obsdate, fits_file_format)
fits_file_format_abs = os.path.join('E:', 'asdetector-data', 'output', 'ramps', obsdate, fits_file_format)

bands = ['YJ', 'HK']

obslog = read_csv(obsdate_file_abs, sep='\t', header=0)
# obslog['output_files'] = obslog['Source'] + '.' + obslog['Exposure Frames'].astype(str) + '.' + obslog['YJ Wheel'] +\
#                          '.' + obslog['HK Wheel'] + '.' + obslog['Slit Wheel'] + '.' + obslog['Aux wheel'] + '.{}.fits'

obslog['output_files'] = obslog['Source'] + '.' + obslog['commanded exposure time (s)'].astype(str) + '.' + obslog['YJ Wheel'] +\
                         '.' + obslog['HK Wheel'] + '.' + obslog['Slit Wheel'] + '.' + obslog['Aux wheel'] + '.{}.fits'

for output_file in obslog.output_files.unique():
    # print(output_file)
    observation_df = obslog[obslog['output_files'] == output_file]
    for band in bands:
        output_file_band = os.path.join(output_directory, output_file.format(band))
        print(output_file_band)
        if not os.path.exists(output_file_band) or overwrite:
            obsfiles = \
                [fits_file_format_abs.format(obsdate, obsnum, band) for obsnum in observation_df['Observation number']]
            observation_df[band] = obsfiles
            source_on_df = observation_df[observation_df['ON?']]
            source_off_df = observation_df[observation_df['ON?'] == False]
            print(source_on_df[band])
            source_on_array = np.mean(np.asarray([fits.getdata(f) for f in source_on_df[band]]), axis=0)
            if not source_off_df.empty:
                source_off_array = np.mean(np.asarray([fits.getdata(f) for f in source_off_df[band]]), axis=0)
            else:
                source_off_array = np.zeros(source_on_array.shape)
                print("zeros shape {}".format(source_off_array.shape))
            output_array = source_on_array - source_off_array
            fits.HDUList(fits.PrimaryHDU(data=output_array)).writeto(output_file_band, overwrite=overwrite)
            print(output_file_band)
