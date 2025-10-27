import os
from glob import glob

from astropy.io import fits
from matplotlib import pyplot as plt

data_dir = r'G:\asdetector-data\output\obslog_reduce\combined_test_set'
os.chdir(data_dir)
crop_settings = {
    'YJ': {
        'Vph30': {
            'long': {
                'x': [1725, 2125],
                'y': [1750, 2350]
            },
            '80 um': {
                'x': [1725, 2125],
                'y': [1975, 2125]
            },
            '130 um': {
                'x': [1725, 2125],
                'y': [1975, 2125]
            },
            '250 um': {
                'x': [1725, 2125],
                'y': [1975, 2125]
            },
        },
        'Vph300': {
            'long': {
                'x': [900, 3200],
                'y': [1700, 2500]
            },
            '80 um': {
                'x': [900, 3200],
                'y': [1900, 2300]
            },
            '130 um': {
                'x': [900, 3200],
                'y': [1900, 2300]
            },
            '250 um': {
                'x': [900, 3200],
                'y': [1900, 2300]
            },
        },
        'grism': {
            '80 um': {
                'x': [600, 3600],
                'y': [950, 3350]
            },
            '130 um': {
                'x': [600, 3600],
                'y': [950, 3350]
            },
            '250 um': {
                'x': [600, 3600],
                'y': [950, 3350]
            },
        }
    },
    'HK': {
        'Vph30': {
            'long': {
                'x': [1900, 2300],
                'y': [1500, 2100]
            },
            '80 um': {
                'x': [1900, 2300],
                'y': [1780, 1930]
            },
            '130 um': {
                'x': [1900, 2300],
                'y': [1780, 1930]
            },
            '250 um': {
                'x': [1900, 2300],
                'y': [1780, 1930]
            },
        },
        'Vph300': {
            'long': {
                'x': [900, 3200],
                'y': [2000, 2800]
            },
            '80 um': {
                'x': [900, 3200],
                'y': [2150, 2650]
            },
            '130 um': {
                'x': [900, 3200],
                'y': [2150, 2650]
            },
            '250 um': {
                'x': [900, 3200],
                'y': [2150, 2650]
            },
        },
        'grism': {
            '80 um': {
                'x': [800, 3800],
                'y': [200, 2600]
            },
            '130 um': {
                'x': [800, 3800],
                'y': [200, 2600]
            },
            '250 um': {
                'x': [800, 3800],
                'y': [200, 2600]
            },
        }
    },
}

cameras = ['HK', 'YJ']
gratings = ['Vph30', 'Vph300', 'grism']
slits = ['long', '80 um', '130 um', '250 um']
sources = [
    ['SLS201L', 'HgArKrXeNe', 'Kr'],
    # ['Kr', 'Ne', 'Xe']
]

# filename = r'G:\asdetector-data\output\obslog_reduce\20230726\argon.9.medresVPH.medresVPH.long slit.open.YJ.fits'
# filename = r'G:\asdetector-data\output\obslog_reduce\20230726\argon.9.medresVPH.medresVPH.long slit.open.HK.fits'
# filename = r'G:\asdetector-data\output\obslog_reduce\20230726\argon.3.lowresVPH.lowresVPH.80 micron slit.open.YJ.fits'
filename = r'G:\asdetector-data\output\obslog_reduce\20230726\argon.3.lowresVPH.lowresVPH.80 micron slit.open.HK.fits'
# title = 'Argon YJ Medium Res'
# title = 'Argon HK Medium Res'
# title = 'Argon YJ Low Res'
title = 'Argon HK Low Res'
regex_format = '{source}.*.{grating}.{grating}.{slit}.open.{camera}.fits'  # SLS201L.20.Vph30.Vph30.250 um.open.HK.fits

subplot_xs = range(2)
subplot_ys = range(3)
subplot_coords = [[(y, x) for y in subplot_ys] for x in subplot_xs]

for grating in gratings:
    for slit in slits:
        if grating == 'grism' and slit == 'long':
            continue
        for source_list in sources:
            fig, axs = plt.subplots(3, 2)
            fig.set_figwidth(8)
            fig.set_figheight(10)
            suptitle = '{} grating with {} slit'.format(grating, slit)
            fig.suptitle(suptitle)
            for camera, subplots in zip(cameras, subplot_coords):
                for source, subplot in zip(source_list, subplots):
                    regex_dict = {'source': source, 'grating': grating, 'slit': slit, 'camera': camera}
                    regex = regex_format.format(**regex_dict)
                    files = glob(regex)
                    image = fits.getdata(files[0])
                    crop = crop_settings[camera][grating][slit]
                    # image = image[crop['y'][0]:crop['y'][1], crop['x'][0]:crop['x'][1]]
                    print(files)
                    ax = axs[subplot[0], subplot[1]]
                    if source == 'SL201L':
                        vmax = 50000
                    else:
                        vmax = 1500
                    ax.imshow(image, vmin=0, vmax=vmax, interpolation='nearest')
                    ax.set_xlim(*crop['x'])
                    ax.set_ylim(*crop['y'])
                    ax.set_xticks([])
                    ax.set_yticks([])
                    # ax.set_vlim(0, vmax)
                    # ax.set_title(source)
                    ax.set_xlabel(camera)
                    ax.set_ylabel(source)
                    # ax.colorbar()
            for ax in axs.flat:
                ax.label_outer()
            plt.tight_layout()
            plt.savefig('jpg\\{}.{}.jpg'.format(grating, slit))
            plt.show()


# vmin = 0
# vmax = 2000
# xlim = (2000, 2400)
# ylim = (2100, 1700)
# data = fits.getdata(filename)
# plt.imshow(data, vmin=vmin, vmax=vmax)
# plt.xlim(*xlim)
# plt.ylim(*ylim)
# plt.colorbar()
# plt.title(title)
# plt.savefig(filename.replace('.fits', '.png'))
# plt.show()


