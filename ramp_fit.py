import os


def calculate_ramp(filenames, exposure_times, save_dir='.', suffix='.dark.fits', overwrite=False):
    assert isinstance(filenames, dict)
    for band, exposures in filenames.items():
        total_noise_arrays = []
        prefix = band + 15*'?'
        assert isinstance(exposures, dict)
        for prefix, exposure in exposures.items():
            save_name = os.path.join(save_dir, prefix + suffix)
            if not os.path.isfile(save_name) and not overwrite:
                print('creating array for {}'.format(prefix))
                slope_array, slope_array_var = calculate_total_noise_exposure(exposure, exposure_times)
                save_array_as_fits(slope_array, save_name)
                save_array_as_fits(slope_array_var, save_name.replace('.fits', '.var.fits'))
                total_noise_arrays.append(slope_array)
            else:
                total_noise_arrays.append(fits.getdata(save_name))
        total_noise = np.median(np.asarray(total_noise_arrays), axis=0)
        save_name = os.path.join(save_dir, prefix[:14] + '.' + band + suffix)
        save_array_as_fits(total_noise, save_name)
