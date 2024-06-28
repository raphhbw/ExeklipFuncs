from spaceKLIP import plotting
from spaceklip_functions import plot_subimages, display_coron_image, display_image
import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as fits
import sys

# plotting data before psf subtraction
# test_file = '/Users/rbendahan/Astronomy/data/JWST/MIRI-1150/HD92945/spaceklip/to_plot/psf_sub/ADI+RDI_NANNU6_NSUBS6_JWST_MIRI_MIRIMAGE_F1140C_NONE_4QPM_1140_MASK1140-KLmodes-all.fits'
# ref = '/Users/rbendahan/Astronomy/data/JWST/MIRI-1150/HD92945/spaceklip/bgsub/jw01668009001_04101_00007_mirimage_calints.fits'
# unsub = '/Users/rbendahan/Astronomy/data/JWST/MIRI-1150/HD92945/spaceklip/bgsub/jw01668007001_04101_00001_mirimage_calints.fits'

# display_coron_image(unsub, ffov=24, unsubt=True)
# # plotting.display_coron_image(test_file)
# plt.show()
# sys.exit()

# plotting grid with unsub, and different substraction methods
star = 'HD206893'

unsub = '/Users/rbendahan/Astronomy/data/JWST/MIRI-1140/{}/spaceklip/to_plot/unsub/'.format(star)
# unsub = '/Users/rbendahan/Astronomy/data/JWST/MIRI-1140/HD92945/spaceklip/bgsub/jw01668009001_04101_00007_mirimage_calints.fits'
sub = '/Users/rbendahan/Astronomy/data/JWST/MIRI-1140/{}/spaceklip/to_plot/psf_sub6/'.format(star)
# sub = '/Users/rbendahan/Astronomy/data/JWST/MIRI-1140/{}/spaceklip/to_plot/psf_sub1/'.format(star)

fig, ax = plot_subimages(imgdirs=[unsub], subdirs=[sub], filts=['F1140C'], submodes=['ADI', 'RDI', 'ADI+RDI'], numKL=[10], 
                   window_size=24, cmaps_list=['gist_heat'],
                   imgVmin=[-1], imgVmax=[2], subVmin=[-1], subVmax=[2],
                   labelpos=[0.02, 0.02],
                    imtext_col='w', showKL=True, useticklabels=True, cbar_textoff=0.96,
                       hspace=0.1, wspace=0.1, target=star)
plt.show()