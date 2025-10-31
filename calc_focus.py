import os
import re

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from scipy.optimize import curve_fit, minimize

show_images = True
show_line_plots = False
show_optimization_fit_plot = True
crop = False
# obsdate = '20250616'
obsdate = '20251007'
# obsdate = '20240805'
# obsdate = '20240524'
focus_data_dir = f'E:\\asdetector-data\\output\\reduced\\{obsdate}\\focus2'
# focus_data_dir = f'E:\\asdetector-data\\output\\reduced\\{obsdate}'
focus_data_files = [os.path.join(focus_data_dir, f) for f in os.listdir(focus_data_dir) if f.endswith('.fits')]
yj_focus_data_files = [f for f in focus_data_files if f.endswith('YJ.fits')]
hk_focus_data_files = [f for f in focus_data_files if f.endswith('HK.fits')]
yj_focus_data_files.sort()
hk_focus_data_files.sort()
focus_data_info = {
    'YJ': {
        'y_range': (1900-30, 1900+30), 'x_range': (2040-60, 2040+60), 'data_files': yj_focus_data_files, 'vmin': 0, 'vmax': 5000
    },
    # 'HK': {'y_range': (2265, 2358), 'x_range': (2176, 2319), 'data_files': hk_focus_data_files}
    'HK': {
        # 'y_range': (1855-30, 1855+30), 'x_range': (2088-80, 2088+80), 'data_files': hk_focus_data_files, 'vmin': 0, 'vmax': 30000
        'y_range': (1805 - 30, 1805 + 30), 'x_range': (2017 - 30, 2017 + 30), 'data_files': hk_focus_data_files,
        'vmin': 0, 'vmax': 15000

    }
}


def gauss(x, amplitude, mean, std_dev, y0):
    return amplitude * np.exp(-((x - mean) / np.sqrt(2) / std_dev)**2) + y0


def tophat(x, hat_level, hat_mid, hat_width, base_level):
    return np.where((hat_mid-hat_width/2. < x) & (x < hat_mid+hat_width/2.), hat_level, base_level)


def parabola(x, a, b, c):
    return a*x*x + b*x + c


def minimize_tophat_fit(x, y):
    def objective(params, _x, _y):
        return np.sum(np.abs(tophat(_x, *params) - _y))
    guess = np.asarray(np.median(x), (np.max(y), 1.5 * np.median(x), np.percentile(y, 10)))
    # print(guess)
    res = minimize(objective, guess, args=(x, y), method='Nelder-Mead')
    # print(res.x)
    return res.x


def minimize_gauss_fit(x, y, guess=None):
    if guess is None:
        y_max = y.max()
        _edges = x[np.isclose(y, 0.5 * y_max, rtol=0.1)]
        try:
            _width = np.max((_edges[-1] - _edges[0], 8))
        except IndexError:
            _width = 8
        guess = (
            y_max,
            x[y == y_max][0],
            _width,
            np.percentile(y, 10)
        )
        print(guess)
    _popt, _ = curve_fit(gauss, x, y, guess)
    return _popt


def crop_image(img, plot_title='', show=True):
    _image = img.copy()
    hori_array = np.median(img, axis=0)
    hori_pix = np.arange(hori_array.shape[0])
    _fit = minimize_gauss_fit(hori_pix, hori_array)
    hori_amplitude, hori_mean, hori_std_dev, hori_y0 = _fit
    if show_line_plots:
        plt.plot(hori_pix, hori_array, label='data')
        plt.plot(hori_pix, gauss(hori_pix, *_fit), label='fit')
        # plt.plot(fit_pix, tophat(fit_pix, *popt), label='fit')
        plt.legend()
        plt.title(str(plot_title))
        plt.show()
    hori_cutoff = np.abs(int(3 * hori_std_dev))
    hori_mean_int = np.abs(int(hori_mean))
    vert_image = img[:, hori_mean_int - hori_cutoff: hori_mean_int + hori_cutoff]
    if show:
        plt.imshow(vert_image)
        plt.title(str(plot_title))
        plt.show()
    vert_array = np.median(vert_image, axis=1)
    vert_pix = np.arange(vert_array.shape[0])
    _fit = minimize_tophat_fit(vert_pix, vert_array)
    if show:
        plt.plot(vert_pix, vert_array, label='data')
        plt.plot(vert_pix, tophat(vert_pix, *_fit), label='fit')
        plt.legend()
        plt.title(str(plot_title))
        plt.show()
    hat_level, hat_mid, hat_width, base_level = _fit
    low_cutoff = int(hat_mid - hat_width/2)
    high_cutoff = int(hat_mid + hat_width/2)
    cutoff_img = _image[low_cutoff:high_cutoff, :]
    if show:
        plt.imshow(cutoff_img)
        plt.title(str(plot_title))
        plt.show()
    return cutoff_img




def main():
    for band, info in focus_data_info.items():
        _guess = None
        fit_info = []
        fit_number = []
        for data_file in info['data_files']:
            # _fit_num = float(re.findall('\d+', data_file)[-1])
            bandpass = data_file.split('.')[-2]
            bandpass_focus_dict = {'YJ': 'FOCUS0', 'HK': 'FOCUS1'}
            _fit_num = fits.getval(data_file, bandpass_focus_dict[bandpass])
            fit_number.append(_fit_num)
            image = fits.getdata(data_file)[info['y_range'][0]:info['y_range'][1], info['x_range'][0]:info['x_range'][1]]
            if show_images:
                plt.imshow(image, vmin=info['vmin'], vmax=info['vmax'])
                plt.title(str(_fit_num))
                plt.show()
            if crop:
                image = crop_image(image, str(_fit_num), show_images)
            fit_array = np.nanmedian(image, axis=0)
            fit_pix = np.arange(fit_array.shape[0])
            popt = minimize_gauss_fit(fit_pix, fit_array, _guess)
            # _guess = popt
            # popt = minimize_tophat_fit(fit_pix, fit_array)
            fit_info.append(np.abs(popt[2]))  # adding standard deviation to list TODO: uncomment
            # fit_info.append(np.std(image))  # TODO: comment out
            # fit_info.append(np.abs(popt[1]))  # adding mean to list
            # fit_info.append(np.abs(popt[0]))  # adding amplitude to list
            if show_line_plots:
                plt.plot(fit_pix, fit_array, label='data')
                # plt.plot(fit_pix, gauss(fit_pix, *popt), label='fit')
                # plt.plot(fit_pix, tophat(fit_pix, *popt), label='fit')
                plt.legend()
                plt.title(str(_fit_num))
                plt.show()
        fit_number = np.asarray(fit_number)
        fit_info = np.asarray(fit_info)
        fit_degree = 4
        p = np.polyfit(fit_number, fit_info, fit_degree)
        if show_optimization_fit_plot:
            plt.plot(fit_number, fit_info)
            plt.plot(fit_number, np.polyval(p, fit_number))
            plt.show()
        _min = minimize(np.poly1d(p), fit_number[fit_info == fit_info.min()][0])
        print(band + " minimum", _min.x)
        _max = minimize(-np.poly1d(p), fit_number[fit_info == fit_info.max()][0])
        print(band + " maximum", _max.x)


main()
