import os
import re

import numpy as np
from astropy.io import fits


def fowler(input_filename, output_filename, fowler_number=2, base_hdu=1, end_hdu=-1):
    input_data = fits.open(input_filename)
    bias_frames = np.asarray([hdu.data for hdu in input_data[base_hdu:base_hdu+fowler_number-1]])
    end_frames = np.asarray([hdu.data for hdu in input_data[end_hdu-fowler_number+1:end_hdu]])
    bias_mean = np.mean(bias_frames, axis=0)
    end_mean = np.mean(end_frames, axis=0)
    dir_path = os.path.dirname(output_filename)
    if not os.path.isdir(dir_path):
        os.makedirs(dir_path)
    fits.HDUList(fits.PrimaryHDU(data=(bias_mean-end_mean))).writeto(output_filename, overwrite=True)


if __name__ == '__main__':
    replacement_dict = {
        '_Leach_001.fits': '.YJ.fits',
        '_Leach2_001.fits': '.HK.fits'
    }

    image_directory = r'F:\Leach data'
    obsdate = '20230316'
    image_prefix = f'{obsdate}.rimas.'
    output_dir = os.path.join(image_directory, 'fowler', obsdate)
    files = [f for f in os.listdir(image_directory) if f.startswith(image_prefix)]
    # files = [f for f in os.listdir(image_directory) if f.endswith('___Leach_001.fits')]
    input_files = [os.path.join(image_directory, f) for f in files]
    print(files)
    output_files = []
    # shortened_exposures = np.arange(0, 8)
    shortened_exposures = []
    for f, input_file in zip(files, input_files):
        # suffix = '___L' + f.split('___L')[-1]
        suffix = '_L' + f.split('_L')[-1]
        output_file = os.path.join(output_dir, f.replace(suffix, replacement_dict[suffix]))
        # output_file = os.path.join(output_dir, f)
        if os.path.isfile(output_file):
            continue
        obs_number = int(re.findall('\d+', f)[1])
        if obs_number in shortened_exposures:
            final_hdu = 8
        else:
            final_hdu = -1
        fowler(input_file, output_file, 2, 1, final_hdu)
        print(output_file)
