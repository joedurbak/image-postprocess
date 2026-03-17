from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt


def crop_data(data):
    """
    Removes science headers and refout from h4rg data
    """
    if data.shape[1] == 4030:
        data = data[:, 6:]
    if data.shape[1] == 4224:
        data = data[:, :-128]
    return data


def get_refpix_arr(data, refpix_border_pix=4, skip_border_pix=0):
    refpix_arr = np.zeros(data.shape, dtype=bool)
    refpix_arr[skip_border_pix:refpix_border_pix] = True
    refpix_arr[:, skip_border_pix:refpix_border_pix] = True
    if skip_border_pix > 0:
        refpix_arr[:, -refpix_border_pix:-skip_border_pix] = True
        refpix_arr[-refpix_border_pix:-skip_border_pix] = True
    else:
        refpix_arr[:, -refpix_border_pix:] = True
        refpix_arr[-refpix_border_pix:] = True
    return refpix_arr


def gen_refpix_histogram(data, bins=1000, refpix_border_pix=4, skip_border_pix=0):
    refpix_arr = get_refpix_arr(data, refpix_border_pix, skip_border_pix)
    hist = np.histogram(data[refpix_arr], bins=bins)
    return hist


def plot_refpix_histogram(hist, channel):
    plt.plot(hist[1][1:], hist[0])
    plt.title('{} refpix var histogram'.format(channel))
    plt.xlabel('variance (ADU)')
    plt.ylabel('pixel count')
    plt.show()


def gen_bad_refpix_map(data, cutoff, outfile, refpix_border_pix=4, skip_border_pix=0):
    refpix_arr = get_refpix_arr(data, refpix_border_pix, skip_border_pix)
    data_refpix = data.copy()
    data_refpix[np.invert(refpix_arr)] = -1000
    bad_refpix_map = np.zeros(data.shape, dtype=np.uint8)
    bad_refpix_map[data_refpix > cutoff] = 1
    fits.HDUList([fits.PrimaryHDU(bad_refpix_map)]).writeto(outfile, overwrite=True)
    return bad_refpix_map


def main():
    var_filename = '~/Downloads/20251209.rimas.0003.YJ.var.fits'
    var_data = fits.getdata(var_filename)
    var_data = crop_data(var_data)
    hist = gen_refpix_histogram(var_data)
    plot_refpix_histogram(hist, 'YJ')
    refpix_map = gen_bad_refpix_map(var_data, -1, '~/Documents/refpix.border2.YJ.fits', skip_border_pix=2)


if __name__ == '__main__':
    main()
