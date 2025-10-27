import os
import math

import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
# from skimage.measure import regionprops

date = '20250628'
data_dir = r'G:\RIMAS-commissioning-data\{}'.format(date)
crop_center = [1988, 1866]
initial_crop_size = 30
final_crop_size = 20
file_numbers = range(29, 44)
filenames = [os.path.join(data_dir, '{}.rimas.{:04d}.HK.fits'.format(date, n)) for n in file_numbers]
plt_columns = 5
plt_row = math.ceil(len(filenames)/plt_columns)


def find_nearest(array, value):
    array = np.asarray(array)
    close_array = np.abs(array - value)
    return np.argwhere(close_array == close_array.min())[0]


def plot3d(array, subplot, title=''):
    ny, nx = array.shape
    x = range(nx)
    y = range(ny)
    X, Y = np.meshgrid(x, y)
    subplot.plot_surface(X, Y, array, cmap=cm.viridis)
    subplot.set_zlim(0, 10000)
    subplot.set_title(title)



fig, axes = plt.subplots(plt_row, plt_columns, sharex=True, sharey=True, subplot_kw={'projection': '3d'})
for f, ax in zip(filenames, axes.flatten()):
    hdu = fits.open(f)[0]
    start_image = hdu.data[
        crop_center[1]-initial_crop_size:crop_center[1]+initial_crop_size,
        crop_center[0]-initial_crop_size:crop_center[0]+initial_crop_size
    ]
    start_image = start_image - np.median(start_image)
    # centroid = regionprops(start_image)[0].centroid
    centroid = find_nearest(start_image, np.percentile(start_image,97))
    final_image = start_image[
        int(centroid[1])-final_crop_size:int(centroid[1])+final_crop_size,
        int(centroid[0])-final_crop_size:int(centroid[0])+final_crop_size
    ]
    plot3d(final_image, ax, hdu.header['COMMENT1'].replace('BD+19 2769 secondary ', '').replace(' focus 725', ''))
    # plt.imshow(final_image)
plt.show()

plt.show()
