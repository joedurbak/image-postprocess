import os

from argparse import ArgumentParser
from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import minimize
from photutils.background import Background2D, MedianBackground
from photutils.segmentation import detect_sources, SourceFinder, SourceCatalog
from astropy.stats import sigma_clipped_stats


def minimize_and_plot(x, y, poly_order=3, plot_fit=True):
    p = np.polyfit(x, y, poly_order)
    _min = minimize(np.poly1d(p), y[x == x.min()][0], bounds=((x.min(), x.max()), ))
    if plot_fit:
        x_order = np.argsort(x)
        x = x[x_order]
        y = y[x_order]
        fig, ax = plt.subplots()
        _min_x_label = ','.join(['{:.1f}'.format(m) for m in _min.x])
        ax.plot(x, y, label='input data', linestyle=':')
        ax.plot(x, np.polyval(p, x), label='fit data', linestyle='-')
        ax.axvline(_min.x, color='grey', linestyle='--', label='min focus {}'.format(_min_x_label))
        # ax.axhline(_min.y, color='grey', linestyle='--', label='min radius {} pix'.format(_min_y))
        ax.grid(True)
        ax.legend()
        plt.show()
    return _min.x


def get_radius(data, flux_frac=0.5, boxsize=400):
    center_coords = [int(d/2) for d in data.shape]
    half_width = int(boxsize/2)
    crop_data = data[center_coords[0]-half_width:center_coords[0]+half_width, center_coords[1]-half_width:center_coords[1]+half_width]
    bkg_estimator = MedianBackground()
    bkg = Background2D(crop_data, (64, 64), filter_size=(3, 3), bkg_estimator=bkg_estimator)
    threshold = bkg.background + (20.0 * bkg.background_rms)
    data_sub = crop_data - bkg.background
    try:
        segm = detect_sources(data_sub, threshold, npixels=5)
        cat = SourceCatalog(data_sub, segm)
    except TypeError:
        return -1
    radii = cat.fluxfrac_radius(flux_frac).value
    diff = cat.fluxfrac_radius(0.8).value - radii
    print(radii)
    print(diff)
    tab=cat.to_table()
    radius = np.median(radii[radii > 1.5])
    return radius


def calc_focus(images, focus_vals, poly_order=3, plot_fit=True):
    radii = np.asarray([get_radius(image) for image in images])
    focus_vals = np.array(focus_vals)
    good_radii = radii > 0
    radii = radii[good_radii]
    focus_vals = focus_vals[good_radii]
    _min = minimize_and_plot(focus_vals, radii, poly_order=poly_order, plot_fit=plot_fit)
    return _min


def get_focus_curve(
        path, camera_names=('YJ','HK'), camera_name_key='CAMNAME',
        focus_objtype_key='OBJTYPE', focus_objyype_val='tel_focus',
        focus_val_key='TELFOCUS', header_ext=0, image_ext=0,
    ):
    print('Focus path: {}'.format(path))
    if os.path.isdir(path):
        print('Found directory: {}'.format(path))
        file_list = [os.path.join(path, f) for f in os.listdir(path) if f.endswith('.fits')]
    else:
        file_list = path.split(',')
        for f in file_list:
            if not os.path.isfile(f):
                raise FileNotFoundError('File not found: {}'.format(f))
    print('Found {} files'.format(file_list))
    cam_image_lists = {cam: ([], []) for cam in camera_names}
    for f in file_list:
        hdus = fits.open(f, mode='readonly')
        hdr = hdus[header_ext].header
        if hdr.get(focus_objtype_key, '').lower() == focus_objyype_val.lower():
            print('{}'.format(f))
            cam_image_lists[hdr[camera_name_key]][0].append(hdus[image_ext].data)
            cam_image_lists[hdr[camera_name_key]][1].append(hdr[focus_val_key])
    focus_vals = []
    for k in camera_names:
        focus_args = cam_image_lists[k]
        if focus_args[0]:
            focus_vals.append(calc_focus(*focus_args))
        else:
            focus_vals.append(-1)
    return focus_vals


def main():
    parser = ArgumentParser()
    parser.add_argument("path", help="path to directory containing focus curve images, or comma separated list of files")
    args = parser.parse_args()
    focus_vals = get_focus_curve(args.path)
    print('Focus values:\n{}'.format(focus_vals))
    print('Average focus:\n{}'.format(np.mean(focus_vals)))


if __name__ == '__main__':
    main()
