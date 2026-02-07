import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from photutils.aperture import CircularAperture, aperture_photometry

# 1. Load the PSF image
# Replace 'proto.fits' with your PSF model image or extracted array
# snap_file = '/Users/jdurbak/Documents/RIMAS_data/test/coadd/NGC7027/psf/snap_coaddNGC7027_20251009T031656570_H.fits.fits'
snap_file = '/Users/jdurbak/Documents/RIMAS_data/test/coadd/NGC7027/psf/snap_coaddNGC7027_20251009T031709307_J.fits.fits'
hdu = fits.open(snap_file)
data = hdu[0].data
header = hdu[0].header

# Ensure data is 2D (sometimes PSFEx check-images have a dummy 3rd dim)
if data.ndim == 3:
    data = data[0]

# 2. Define the center and the radii to measure
# PSFEx models are usually centered exactly in the middle of the vignette
cy, cx = np.array(data.shape) / 2.0 - 0.5
radii = np.arange(0.5, 12.5, 0.5)  # From 0.5 to 15 pixels

# 3. Calculate flux in concentric apertures
ee_flux = []
for r in radii:
    aperture = CircularAperture((cx, cy), r=r)
    phot_table = aperture_photometry(data, aperture)
    ee_flux.append(phot_table['aperture_sum'][0])

# 4. Normalize the Encircled Energy (0.0 to 1.0)
ee_normalized = np.array(ee_flux) / np.max(ee_flux)

# 5. Plot the results
plt.figure(figsize=(8, 5))
plt.plot(radii, ee_normalized, marker='o', linestyle='-', color='darkblue')
plt.axhline(0.5, color='red', linestyle='--', label='50% EE')
plt.axhline(0.8, color='green', linestyle='--', label='80% EE')
plt.xlabel('Radius (pixels)')
plt.ylabel('Fraction of Enclosed Energy')
plt.title('PSF Encircled Energy Profile - J band')
plt.grid(alpha=0.3)
plt.legend()
plt.show()

# Print specific metrics
r80 = radii[np.where(ee_normalized >= 0.8)[0][0]]
print(f"80% Encircled Energy Radius: {r80:.2f} pixels")
