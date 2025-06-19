import os

from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt

_raw_dir_format = r'E:\asdetector-data\output\raw\{date}'
_bands = ('YJ', 'HK')
_raw_file_format = '{date}.rimas.{exposure:04d}.{band}.{index:04d}.fits'


def ecdf_f(x):
    return lambda t: (1.0/x.shape[0])*x[np.logical_or((x < t), (x == t))].shape[0]


def plot_img_histogram(data, sig=3, nbins=50, save=False, size=(8, 6), xlabel='Data', unit='counts', title=None,
                       yscale=None, xscale=None, ax_range=None):
    """Plots histogram and saves PDF of given image data."""

    # Read in image
    data = data[1000:2900, 1000:2600]
    data = data.squeeze()

    # Calculate statistics
    mean_data = np.nanmean(data)
    median_data = np.nanmedian(data)
    percentile = np.nanpercentile(data, [16, 84])
    per_range = percentile[1] - percentile[0]
    print("Percentile is ", percentile)
    print(" Range is %f" % (percentile[1] - percentile[0]))
    if ax_range is None:
        edges = np.linspace(median_data - sig * per_range, median_data + sig * per_range, nbins)
    else:
        edges = np.linspace(ax_range[0], ax_range[1], nbins)
    print("Mean data is %.3f\nMedian data is %.3f" % (mean_data, median_data))

    # Plot
    fig, axh = plt.subplots(figsize=size)
    n, bins, patches = plt.hist(data.flatten(), bins=edges, histtype='step')
    if title is not None:
        plt.title(title)
    axh.set_xlabel(xlabel)
    axh.set_ylabel("Frequency")
    axh.text(0.03, 0.8,
             "$\mathrm{Mean} =\/ %.3f$ %s\n$\mathrm{Median} = \/%.3f$ %s" % (mean_data, unit, median_data, unit),
             transform=axh.transAxes)

    # Plot cumulative function
    ecdf = ecdf_f(data.flatten())
    if ax_range is None:
        t = np.linspace(median_data - sig * per_range, median_data + sig * per_range, 2 * nbins)
    else:
        t = np.linspace(ax_range[0], ax_range[1], 2 * nbins)
    cumul = [ecdf(i) for i in t]
    ax2 = axh.twinx()
    ax2.step(t, cumul, 'r-', label='Cumulative Distribution')
    ax2.set_ylabel(r'Cumulative')
    if ax_range is None:
        plt.xlim(median_data - sig * per_range, median_data + sig * per_range)
    else:
        plt.xlim(ax_range[0], ax_range[1])
    axh.ticklabel_format(axis='y', style='sci', scilimits=(0, 6), useOffset=False)
    axh.ticklabel_format(axis='x', useOffset=False)

    if yscale is not None:
        axh.set_yscale(yscale)
        ax2.set_yscale(yscale)
    if xscale is not None:
        axh.set_xscale(xscale)
        axh.set_xscale(xscale)

    plt.legend(loc=0)
    plt.tight_layout()

    if save:
        name = title + '_hist.pdf'
        plt.savefig(name)
    else:
        plt.show()


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


def get_filenames_dir(data_dir):
    filenames = {}
    ls = [f for f in os.listdir(data_dir) if f.endswith('.fits') and not f.startswith('reset')]
    for f in ls:
        date, instrument, exposure, band, index, suffix = f.split('.')
        prefix = '.'.join((date, instrument, exposure, band))
        if band not in filenames.keys():
            filenames[band] = {}
        if prefix not in filenames[band].keys():
            filenames[band][prefix] = {'files': [], 'indices': []}
        filenames[band][prefix]['files'].append(os.path.join(data_dir, f))
        filenames[band][prefix]['files'].sort()
        filenames[band][prefix]['indices'].append(int(index))
        filenames[band][prefix]['files'].sort()
    return filenames


def save_array_as_fits(array, filename):
    print('Writing file {}'.format(filename))
    if not os.path.isdir(os.path.dirname(filename)):
        os.makedirs(os.path.dirname(filename))
    fits.HDUList([fits.PrimaryHDU(data=array)]).writeto(filename, overwrite=True)


def ref_pixel_correct_frame(frame, readout_channels=33, tel_columns=True):
    if tel_columns:
        frame = frame[:, 6:]
    ref_pixels_0 = frame[2:4, :]  # acquiring the desired reference pixels
    ref_pixels_1 = frame[-4:-2, :]  # acquiring the desired reference pixels
    ref_pixels = np.vstack((ref_pixels_0, ref_pixels_1))  # stacking reference pixels
    channel_width = int(frame.shape[1] / readout_channels)  # getting the width of each readout channel
    channel_starts = range(0, frame.shape[1], channel_width)  # getting the starting pixel value of each channel
    subtraction_row = np.zeros(frame.shape[1])  # initializing row array
    for channel_start in channel_starts:
        # finding the median reference pixel value for each channel and setting the subtraction reference
        # subtraction for each channel based on that
        subtraction_row[channel_start:channel_start+channel_width] = \
            np.median(ref_pixels[:, channel_start:channel_start+channel_width])
    subtraction_array = np.zeros(frame.shape)
    # subtraction_array[:] = subtraction_row  # creating subtraction array
    frame = frame.astype('float') - subtraction_array
    return frame[:, :-channel_width]


def get_ref_pix_corrected_array(filename):
    frame = fits.getdata(filename)
    # print(frame.shape)
    saturated_pixels = frame[4:4092, 4:4092] > 40000
    saturation_percentage = 100 * np.sum(saturated_pixels) / np.prod(saturated_pixels.shape)
    print(filename, saturation_percentage)
    return ref_pixel_correct_frame(frame)


def get_cds_pair_array(exposure_list):
    start_indices = np.arange(0, len(exposure_list), 2)
    cds_pair_array = []
    for start_index in start_indices:
        start_name = exposure_list[start_index]
        end_name = exposure_list[start_index+1]
        print('creating cds pair {} - {}'.format(end_name, start_name))
        start_array = get_ref_pix_corrected_array(start_name)
        end_array = get_ref_pix_corrected_array(end_name)
        cds_pair_array.append(end_array-start_array)
    return np.asarray(cds_pair_array)


def calculate_cds_noise_frame(filenames, save_dir='.', suffix='.cds.fits', overwrite=False, gain=1.0):
    assert isinstance(filenames, dict)
    for band, exposures in filenames.items():
        cds_noise_arrays = []
        assert isinstance(exposures, dict)
        for prefix, exposure in exposures.items():
            save_name = os.path.join(save_dir, prefix+suffix)
            if not os.path.isfile(save_name) or overwrite:
                print('creating array for {}'.format(prefix))
                cds_pair_array = get_cds_pair_array(exposure)
                cds_noise_array = np.std(cds_pair_array, axis=0) * gain
                cds_noise_arrays.append(cds_noise_array)

                save_array_as_fits(cds_noise_array, save_name)
            else:
                cds_noise_arrays.append(fits.getdata(save_name))
        cds_noise = np.median(np.asarray(cds_noise_arrays), axis=0)
        save_name = os.path.join(save_dir, prefix[:14]+'.'+band+suffix)
        save_array_as_fits(cds_noise, save_name)


def calculate_total_noise_exposure(exposure, exposure_times):
    # exposure_arrays = np.asarray([fits.getdata(f)[:, 6: -128] for f in exposure])
    exposure_arrays = np.asarray([get_ref_pix_corrected_array(f) for f in exposure])
    reshaped_exposures = exposure_arrays.reshape((exposure_times.shape[0], -1))
    # print(reshaped_exposures.shape)
    # print(exposure_times.shape)
    try:
        popt, pcov = np.polyfit(exposure_times, reshaped_exposures, deg=1, cov=True)
        # print(popt.shape, pcov.shape)
    except ValueError:
        popt = np.polyfit(exposure_times, reshaped_exposures, deg=1, cov=False)
        # print(popt.shape)
        pcov = np.zeros((2, *popt.shape))
        # print(pcov.shape)
    # print("popt shape is", popt[0].shape)
    slope_arr = popt[0].reshape((exposure_arrays.shape[1], -1))
    # print("Slope arr shape is: ", slope_arr.shape)
    slope_var_arr = pcov[0, 0, :].reshape((exposure_arrays.shape[1], -1))
    return slope_arr, slope_var_arr


def calculate_total_noise_frame(
        filenames, exposure_times, save_dir='.', suffix='.total.fits', overwrite=False, gain=1.0
):
    assert isinstance(filenames, dict)
    for band, exposures in filenames.items():
        total_noise_arrays = []
        total_noise_var_arrays = []
        prefix = band + 15*'?'
        assert isinstance(exposures, dict)
        for prefix, exposure in exposures.items():
            save_name = os.path.join(save_dir, prefix + suffix)
            if not os.path.isfile(save_name) or overwrite:
                print('creating array for {}'.format(prefix))
                slope_array, slope_array_var = calculate_total_noise_exposure(exposure, exposure_times)
                save_array_as_fits(slope_array*gain, save_name)
                save_array_as_fits(slope_array_var, save_name.replace('.fits', '.var.fits'))
                total_noise_arrays.append(slope_array)
            else:
                total_noise_arrays.append(fits.getdata(save_name))
        dt = exposure_times[-1] - exposure_times[0]
        total_noise = np.std(np.asarray(total_noise_arrays), axis=0, ddof=1) * dt
        save_name = os.path.join(save_dir, prefix[:14] + '.' + band + suffix)
        save_array_as_fits(total_noise, save_name)


def calculate_exposure_ramps(
        data_dir, output_dir, frame_time_s=2.86, suffix='.ramp.fits', overwrite=False, start_frame=0
):
    filenames = get_filenames_dir(data_dir)
    for band, exposures in filenames.items():
        for prefix, exposure in exposures.items():
            output_name = os.path.join(output_dir, prefix+suffix)
            print(output_name)
            if not os.path.isfile(output_name) or overwrite:
                exposure_times = np.asarray(exposure['indices']) * frame_time_s
                exposure_times = exposure_times[start_frame:]
                files = exposure['files'][start_frame:]
                try:
                    slope, var = calculate_total_noise_exposure(files, exposure_times)
                except np.linalg.LinAlgError:
                    print('LinAlgError, skipping...')
                    continue
                save_array_as_fits(slope, output_name)
                save_array_as_fits(var, output_name.replace('.fits', '.var.fits'))
            else:
                print('File exists, skipping...')


def calculate_ramp(filenames, exposure_times, save_dir='.', suffix='.dark.fits', overwrite=False, gain=1.0):
    assert isinstance(filenames, dict)
    for band, exposures in filenames.items():
        total_noise_arrays = []
        total_noise_var_arrays = []
        prefix = band + 15*'?'
        assert isinstance(exposures, dict)
        for prefix, exposure in exposures.items():
            save_name = os.path.join(save_dir, prefix + suffix)
            if not os.path.isfile(save_name) or overwrite:
                print('creating array for {}'.format(prefix))
                slope_array, slope_array_var = calculate_total_noise_exposure(exposure, exposure_times)
                save_array_as_fits(slope_array*gain, save_name)
                save_array_as_fits(slope_array_var*gain, save_name.replace('.fits', '.var.fits'))
                total_noise_arrays.append(slope_array)
                total_noise_var_arrays.append(slope_array_var)
            else:
                total_noise_arrays.append(fits.getdata(save_name))
        total_noise = np.median(np.asarray(total_noise_arrays), axis=0)
        total_noise_var = np.median(np.asarray(total_noise_var_arrays), axis=0)
        save_name = os.path.join(save_dir, prefix[:14] + '.' + band + suffix)
        save_array_as_fits(total_noise, save_name)
        save_array_as_fits(total_noise_var, save_name.replace('.fits', '.var.fits'))


def calculate_mean_std(
        filenames, save_dir='.', std_suffix='.std.fits', mean_suffix='.mean.fits', mid_point=None, box_size=None
):
    for band, exposures in filenames.items():
        n_indices = len(list(exposures.values())[0])
        prefix = list(exposures.keys())[0]
        for index in range(1, n_indices):
            same_time_exposures = [(exposure[index], exposure[0]) for exposure in exposures.values()]
            same_time_exposure_cube = []
            for f in same_time_exposures:
                cds_array = get_ref_pix_corrected_array(f[0]).astype(np.double) - \
                            get_ref_pix_corrected_array(f[1]).astype(np.double)
                if mid_point is not None and box_size is not None:
                    cds_array = cds_array[
                        mid_point[1]-box_size//2:mid_point[1]+box_size//2,
                        mid_point[0]-box_size//2:mid_point[0]+box_size//2
                    ]
                same_time_exposure_cube.append(cds_array)
            same_time_exposure_cube = np.asarray(same_time_exposure_cube)
            _std = same_time_exposure_cube.std(axis=0, ddof=1)
            std_name = os.path.join(save_dir, prefix+'.{:04d}'.format(index)+std_suffix)
            save_array_as_fits(_std, std_name)
            _mean = same_time_exposure_cube.mean(axis=0)
            mean_name = os.path.join(save_dir, prefix + '.{:04d}'.format(index) + mean_suffix)
            save_array_as_fits(_mean, mean_name)


def calculate_gain(mean_std_dir, slice_middles, box_widths, std_suffix='.std.fits', mean_suffix='.mean.fits'):
    ls = [f for f in os.listdir(mean_std_dir) if f.endswith('.fits')]
    prefixes = set([f[:19] for f in ls])
    print(prefixes)
    gains = []
    for prefix, middle, width in zip(prefixes, slice_middles, box_widths):
        mean_files = [os.path.join(mean_std_dir, f) for f in ls if f.endswith(mean_suffix) and f.startswith(prefix)]
        mean_files.sort()
        means = np.asarray(
            [np.median(fits.getdata(f)[middle[1]-width//2:middle[1]+width//2, middle[0]-width//2:middle[0]+width//2]) for f in mean_files]
        )
        std_files = [os.path.join(mean_std_dir, f) for f in ls if f.endswith(std_suffix) and f.startswith(prefix)]
        std_files.sort()
        stds = np.asarray(
            [np.median(fits.getdata(f)[middle[1] - width // 2:middle[1] + width // 2,
             middle[0] - width // 2:middle[0] + width // 2]**2) for f in std_files]
        )
        fit = np.polyfit(means, stds, 1)
        print(fit)
        gains.append(1/fit[0])
        frames = np.arange(0, means.shape[0])
        mean_fit = np.polyfit(frames, means, 1)
        plt.plot(frames, means, 'o')
        plt.plot(frames, np.polyval(mean_fit, frames))
        plt.xlabel("Frame Number")
        plt.ylabel("Intensity (ADU)")
        plt.title('Mean intensity of {} superpixel over time'.format(width))
        # plt.plot(means-np.polyval(mean_fit, frames))
        plt.show()
        plt.plot(means, stds, '.', label='measurement')
        plt.plot(means, np.polyval(fit, means), '-', label='fit')
        plt.xlabel("Mean Intensity (ADU)")
        plt.ylabel("Mean Std^2 (ADU^2)")
        plt.title("Intensity vs Std^2 of {} superpixel".format(width))
        plt.show()
    return gains


def plot_noise_distribution(filenames, save_dir, suffix, xlabel, unit='counts'):
    for band, exposures in filenames.items():
        prefix = band + 15 * '?'
        for prefix, exposure in exposures.items():
            pass
        save_name = os.path.join(save_dir, prefix[:14] + '.' + band + suffix)
        data_array = fits.getdata(save_name)
        title = os.path.basename(save_name)
        print(title)
        plot_img_histogram(data_array, unit=unit, title=title, xlabel=xlabel)


def plot_ramp_linearity(filenames, middles, box_widths, plot_individual_ramps=False):
    for band, exposures in filenames.items():
        slopes = []
        ends = []
        for prefix, exposure in exposures.items():
            middle = middles[band]
            width = box_widths[band]
            means = np.asarray(
                [fits.getdata(f)[middle[1] - width // 2:middle[1] + width // 2,
                 middle[0] - width // 2:middle[0] + width // 2].mean() for f in exposure]
            )
            ends.append(means[-1])
            frames = np.arange(0, means.shape[0])
            mean_fit = np.polyfit(frames, means, 1)
            slopes.append(mean_fit[0])
            if plot_individual_ramps:
                plt.plot(means - np.polyval(mean_fit, frames))
                plt.title(str(mean_fit))
                plt.xlabel("Frame")
                plt.ylabel("deviation from linear fit")
                plt.show()
        plt.plot(ends)
        plt.title('intensity for each exposure')
        plt.show()
        plt.plot(slopes)
        plt.title('slopes for each exposure')
        plt.show()


def cubify(
        filenames, outfile
):
    """

    :param filenames: list of files to combine into a cube
    :return:
    """
    arrays = np.asarray([fits.getdata(f) for f in filenames])
    # save_array_as_fits(arrays, outfile)
    if not os.path.isdir(os.path.dirname(outfile)):
        os.makedirs(os.path.dirname(outfile))
    fits.writeto(outfile, arrays, overwrite=True)


if __name__ == '__main__':
    frame_time = 2.86
    theoretical_gain = 1.8
    # cds_dates = [20230801]
    # cds_exposures = [np.arange(2, 7)]
    # cds_indices = np.arange(1, 101)
    cds_dates = [20240724]
    cds_exposures = [np.arange(0, 5)]
    cds_indices = np.arange(0, 100)
    cds_save_dir = os.path.join(_raw_dir_format.format(date=cds_dates[-1]), 'cds')
    cds_filenames = get_filenames(cds_dates, cds_exposures, cds_indices)

    # calculate_cds_noise_frame(cds_filenames, save_dir=cds_save_dir, gain=theoretical_gain)
    # plot_noise_distribution(cds_filenames, cds_save_dir, '.cds.fits', 'CDS noise (e-)', unit='e-')

    # dark_dates = [20230731, 20230801]
    # dark_dates = [20240223]
    # dark_dates = [20240229]
    dark_dates = [20240724]
    dark_exposures = [
        # [4, 5],
        # [0, 1]
        # [0, 1, 2, 3]
        [105, 106, 107, 108]
    ]
    # dark_indices = np.arange(0, 2626+1, 26)
    # dark_indices = np.arange(0, 2525+1, 25)
    dark_indices = np.arange(0, 2600+1, 25)
    # dark_indices = np.arange(0, 106, 1)
    dark_filenames = get_filenames(dark_dates, dark_exposures, dark_indices)
    dark_save_dir = os.path.join(_raw_dir_format.format(date=dark_dates[-1]), 'dark')
    dark_exposure_times = dark_indices * frame_time
    # calculate_ramp(dark_filenames, dark_exposure_times, save_dir=dark_save_dir, gain=theoretical_gain)
    # plot_noise_distribution(dark_filenames, dark_save_dir, '.dark.fits', unit='e-/s', xlabel='Dark Current (e-/s)')
    # plot_noise_distribution(dark_filenames, dark_save_dir, '.dark.var.fits', unit='e-/s')

    # total_noise_dates = [20230801]
    # total_noise_exposures = [np.arange(7, 107)]
    # total_noise_indices = np.arange(1, 56)
    total_noise_dates = [20240724]
    total_noise_exposures = [np.arange(5, 105)]
    total_noise_indices = np.arange(1, 55)
    noise_save_dir = os.path.join(_raw_dir_format.format(date=total_noise_dates[-1]), 'noise')
    total_noise_exposure_times = np.arange(
        total_noise_indices[0]*frame_time, (total_noise_indices[-1]+0)*frame_time+1, frame_time
    )
    total_noise_filenames = get_filenames(total_noise_dates, total_noise_exposures, total_noise_indices)
    # calculate_total_noise_frame(
    #     total_noise_filenames, total_noise_exposure_times, save_dir=noise_save_dir, gain=theoretical_gain)
    # plot_noise_distribution(total_noise_filenames, noise_save_dir, '.total.fits', xlabel='Total Noise (e-)', unit='e-')

    # gain_exposures = [np.arange(3, 103)]
    # gain_exposures = [np.arange(10, 110)]  # HK run numbers
    gain_exposures = [np.arange(120, 219)]  # YJ run numbers
    gain_indices = np.arange(2, 22)
    gain_dates = [20240820]
    gain_directory = r'G:\asdetector-data\output\raw\{date}'
    gain_bands = ['YJ']
    gain_filenames = get_filenames(gain_dates, gain_exposures, gain_indices, gain_directory, bands=gain_bands)
    gain_save_dir = os.path.join(
        gain_directory.format(date=gain_dates[-1]), 'gain',
        'exposures_{}-{}'.format(gain_exposures[0][0], gain_exposures[0][-1]))
    # calculate_mean_std(gain_filenames, gain_save_dir, mid_point=(2062,1869), box_size=1000)  # HK coords
    # calculate_mean_std(gain_filenames, gain_save_dir, mid_point=(1828, 2058), box_size=1000)  # YJ coords
    gains = calculate_gain(gain_save_dir, [(310, 700)], [128])
    print(gains)
    # alphas = np.asarray([0.0161])  # my calculation based on all 8 neighbors: 0.0161. DCL calculations based on 4: 0.023 HK
    alphas = np.asarray([0.0139])  # my calculation based on 4 neighbors: 0.0139 YJ
    print(np.asarray(gains)*(1-4*alphas))
    # plot_ramp_linearity(gain_filenames, middles={'HK': (2032, 1847)}, box_widths={'HK': 256})

    # ramp_files = get_filenames_dir(r'E:\asdetector-data\output\raw\20230803')
    # print(ramp_files)
    # raw_dirs = [
        # r'E:\asdetector-data\output\raw\20230802',
        # r'E:\asdetector-data\output\raw\20230803'
        # r'E:\asdetector-data\output\raw\20240223'
    # ]
    # output_dirs = [
        # r'E:\asdetector-data\output\ramps\20230802',
        # r'E:\asdetector-data\output\ramps\20230803'
        # r'.\darks'
    # ]
    # for raw_dir, output_dir in zip(raw_dirs, output_dirs):
    #     pass
        # calculate_exposure_ramps(raw_dir, output_dir, overwrite=True, start_frame=1)
    cube_data_dir = r'G:\asdetector-data\output\raw\20240820\gain\exposures_10-109'
    # cube_files = [os.path.join(cube_data_dir, f) for f in os.listdir(cube_data_dir) if f.endswith('.mean.fits')]
    # cubify(cube_files, os.path.join(cube_data_dir, 'cube', '20240820.rimas.0010-0109.HK.mean.fits'))
