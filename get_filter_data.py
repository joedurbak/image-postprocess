import os

from astropy.io import fits

parent_dir = r'G:\asdetector-data\output\raw'
headers = ('OBJNAME', 'OBJTYPE', 'FILTER1', 'FILTER2', 'FILTER3', 'FILTER4',)

sub_dirs = [os.path.join(parent_dir, d) for d in os.listdir(parent_dir) if not d.startswith('_')]
print(f'FILENAME\t' + '\t'.join(headers))

for sub_dir in sub_dirs:
    ls = [os.path.join(sub_dir, f) for f in os.listdir(sub_dir) if f.endswith('.fits')]
    for f in ls:
        try:
            hdr = fits.getheader(f)
            f_data = [hdr[k].strip() for k in headers]
            frame = hdr['FRAME']
            adr = hdr['ASICADDR']
        except (KeyError, OSError):
            continue
        filter_3 = f_data[-2].upper()
        filter_4 = f_data[-1].upper()
        if filter_3 == 'OPEN' and filter_4 == 'PINHOLE' and frame==0 and adr == 2:
        # if frame == 0 and adr == 2:
            print(f'{f}\t' + '\t'.join(f_data))
