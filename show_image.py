from astropy.io import fits
from matplotlib import pyplot as plt


# filename = r'G:\asdetector-data\output\obslog_reduce\20230726\argon.9.medresVPH.medresVPH.long slit.open.YJ.fits'
# filename = r'G:\asdetector-data\output\obslog_reduce\20230726\argon.9.medresVPH.medresVPH.long slit.open.HK.fits'
# filename = r'G:\asdetector-data\output\obslog_reduce\20230726\argon.3.lowresVPH.lowresVPH.80 micron slit.open.YJ.fits'
filename = r'G:\asdetector-data\output\obslog_reduce\20230726\argon.3.lowresVPH.lowresVPH.80 micron slit.open.HK.fits'
# title = 'Argon YJ Medium Res'
# title = 'Argon HK Medium Res'
# title = 'Argon YJ Low Res'
title = 'Argon HK Low Res'
vmin = 0
vmax = 2000
xlim = (2000, 2400)
ylim = (2100, 1700)
data = fits.getdata(filename)
plt.imshow(data, vmin=vmin, vmax=vmax)
plt.xlim(*xlim)
plt.ylim(*ylim)
plt.colorbar()
plt.title(title)
plt.savefig(filename.replace('.fits', '.png'))
plt.show()


