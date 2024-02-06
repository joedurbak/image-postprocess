import argparse
import os

from astropy.io import fits
import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--obs', help='observation file prefix being processed')
    parser.add_argument('-b', '--baseframe', help='frame index (base 1) to subtract, default=2', default=2)
    parser.add_argument('-d', '--dark', default=None, help='file prefix for frames to be subtracted from the obs frame')
    parser.add_argument('--datadir', default=r'F:\Leach data')
    parser.add_argument('--outdir', default=r'F:\Leach data\processed')
    parser.add_argument('--outputfile', default=None)
    args = parser.parse_args()

    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)
    # observation = fits.open(os.path.join(args.datadir, args.obs))
    observation_file = os.path.join(args.datadir, args.obs)
    file_suffixes = ('_Leach_001.fits', '_Leach2_001.fits')
    observations = [fits.open(observation_file+suffix) for suffix in file_suffixes]
    if args.outputfile is None:
        output_files = [os.path.join(args.outdir, args.obs+suffix) for suffix in file_suffixes]
    else:
        output_files = [os.path.join(args.outputfile, args.obs + suffix) for suffix in file_suffixes]
    if args.dark is not None:
        dark_file = os.path.join(args.datadir, args.dark)
        darks = [fits.open(dark_file+suffix) for suffix in file_suffixes]
    else:
        darks = []
        for observation in observations:
            _dark = [fits.PrimaryHDU()]
            for _frame in observation[1:]:
                arr = np.zeros(_frame.data.shape).astype(np.int32)
                _dark.append(fits.ImageHDU(arr))
            darks.append(fits.HDUList(_dark))
    for observation, dark, outfile in zip(observations, darks, output_files):
        n_frames = len(observation)
        hdus = [fits.PrimaryHDU()]
        base_data = observation[int(args.baseframe)].data.astype(np.int32) - \
            dark[int(args.baseframe)].data.astype(np.int32)
        for _frame in range(1, n_frames):
            dark_sub = observation[_frame].data.astype(np.int32) - dark[_frame].data.astype(np.int32)
            hdus.append(fits.ImageHDU(base_data - dark_sub))
        fits.HDUList(hdus).writeto(outfile, overwrite=True)


if __name__ == '__main__':
    main()
