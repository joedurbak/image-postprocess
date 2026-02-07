import os

import numpy as np
from astropy.io import fits


_raw_dir_format = r'E:\asdetector-data\output\raw\{date}'
_bands = ('YJ', 'HK')
_raw_file_format = '{date}.rimas.{exposure:04d}.{band}.{index:04d}.fits'


def cubify(filelist):
    cube_arr = np.asarray([fits.getdata(f) for f in filelist], dtype=np.uint16)
    if cube_arr.shape[2] == 4230:
        cube_arr = cube_arr[:,:,6:]  # cropping telemetry headers
    hdr = fits.getheader(filelist[-1])
    return cube_arr, hdr


def save_cube(cube_arr, hdr, outfile, overwrite=False):
    fits.HDUList([fits.PrimaryHDU(header=hdr), fits.ImageHDU(cube_arr)]).writeto(outfile, overwrite=overwrite)


def cubify_filenames(filenames, save_dir, suffix='.cube.fits', overwrite=False):
    assert isinstance(filenames, dict)
    for band, exposures in filenames.items():
        prefix = band + 15*'?'
        assert isinstance(exposures, dict)
        for prefix, exposure in exposures.items():
            save_name = os.path.join(save_dir, prefix + suffix)
            if not os.path.isfile(save_name) or overwrite:
                print('creating cube for {}'.format(prefix))
                cube_arr, hdr = cubify(exposure)
                save_cube(cube_arr, hdr, save_name, overwrite=overwrite)


def get_filenames(dates, exposures, indices, dir_format=_raw_dir_format, file_format=_raw_file_format, bands=_bands):
    assert len(dates) == len(exposures)
    filenames = {}
    combined_format = os.path.join(dir_format, file_format)
    for band in bands:
        filenames[band] = {}
        for date, exposure_list in zip(dates, exposures):
            for exposure in exposure_list:
                exposure_filenames = []
                for index in indices:
                    filename = combined_format.format(date=date, band=band, index=index, exposure=exposure)
                    if not os.path.exists(filename):
                        print('{} does not exist. Skipping...'.format(filename))
                    else:
                        exposure_filenames.append(filename)
                    filenames[band]['{date}.rimas.{exposure:04d}.{band}'.format(date=date, exposure=exposure, band=band)] = exposure_filenames
    return filenames


def main():
    output_dir = r'G:\RIMAS_data\irrc_files\darks'
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    total_noise_dates = [20240724]
    total_noise_exposures = [np.arange(5, 105)]
    total_noise_indices = np.arange(0, 55)
    total_noise_filenames = get_filenames(total_noise_dates, total_noise_exposures, total_noise_indices)
    cubify_filenames(total_noise_filenames, output_dir, overwrite=True)


if __name__ == '__main__':
    main()
