from astropy.io import fits

date = 20240909
on = 21
off = 20
# exposure = 29
off_scale = 1.0
yj = 'YJ'
hk = 'HK'
suffix = '.sub'
# f_format = r'G:\asdetector-data\output\raw\{}\{}.rimas.{:04d}.{}.{:04d}.fits'
f_format = r'G:\asdetector-data\output\reduced\{}\{}.rimas.{:04d}.{}.fits'

yj_sub = fits.getdata(f_format.format(date, date, on, yj)).astype(float) - fits.getdata(f_format.format(date, date, off, yj)).astype(float) * off_scale
fits.writeto(f_format.format(date, date, on, yj+suffix), yj_sub, overwrite=True)
hk_sub = fits.getdata(f_format.format(date, date, on, hk)).astype(float) - fits.getdata(f_format.format(date, date, off, hk)).astype(float) * off_scale
fits.writeto(f_format.format(date, date, on, hk+suffix), hk_sub, overwrite=True)
