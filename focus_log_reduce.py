import os

import numpy as np
from astropy.io import fits
from pandas import read_csv

obsdate = '20230810'
logfile = f'C:\\PythonProjects\\image-postprocess\\focus_logs\\focus_log_{obsdate}.csv'
bands = ['YJ', 'HK']
# data_dir = f'F:\\Leach data\\fowler\\{obsdate}'
data_dir = f'E:\\asdetector-data\\output\\reduced\\{obsdate}'
# output_dir = data_dir.replace('fowler', 'focus')
output_dir = os.path.join(data_dir, 'focus')
if not os.path.isdir(output_dir):
    os.makedirs(output_dir)
file_format = os.path.join(data_dir, '{}.rimas.{:04d}.{}.fits')

log_df = read_csv(logfile)
if 'abs_pos' in log_df.keys():
    log_df['abs_pos_yj'] = log_df['abs_pos']
    log_df['abs_pos_hk'] = log_df['abs_pos']
obs_nums = log_df['obs_num']  # FIXME: remove "+ 130"
for band in bands:
    log_df[band] = [file_format.format(obsdate, obsnum, band) for obsnum in obs_nums]
unique_positions_yj = log_df['abs_pos_yj'].unique()
unique_positions_hk = log_df['abs_pos_hk'].unique()
unique_positions_all = (unique_positions_yj, unique_positions_hk)
save_name_format = os.path.join(output_dir, 'f.{:04d}.{}.fits')

for band, unique_positions in zip(bands, unique_positions_all):
    for pos in unique_positions:
        save_name = save_name_format.format(pos, band)
        if os.path.isfile(save_name):
            continue
        pos_df = log_df[log_df['abs_pos_{}'.format(band.lower())] == pos]
        # save_name = os.path.join(output_dir, 'f.{:04d}.{}.fits')
        on_frames = []
        off_frames = []
        for i, row in pos_df.iterrows():
            image = fits.getdata(row[band])
            if row['on_off'].lower() == 'on':
                on_frames.append(image)
            else:
                off_frames.append(image)
        print(save_name)
        output_image = np.median(np.asarray(on_frames), axis=0) - np.median(np.asarray(off_frames), axis=0)
        fits.HDUList(fits.PrimaryHDU(data=output_image)).writeto(save_name, overwrite=True)
