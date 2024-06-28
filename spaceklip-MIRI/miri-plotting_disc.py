from spaceKLIP import plotting
from spaceklip_functions import plot_subimages, display_coron_image, display_image
import matplotlib.pyplot as plt
import numpy as np
# import astropy.io.fits as fits
import sys

# star = 'HD206893'

sub = '/Users/rbendahan/Astronomy/data/JWST/MIRI-1140/Cycle1/sub6/ADI+RDI/'
# sub = '/Users/rbendahan/Astronomy/data/JWST/MIRI-1140/Cycle1/sub1/RDI/'

discdata = {'HD92945':{'rin':54.,
                        'rout':133.,
                        'gapin':62.6,
                        'gapout':81.5,
                        'PA':100.,
                        'inc':65.4,
                        'dpc': 21.5},
            'HD107146':{'rin':44.,
                        'rout':144.3,
                        'gapin':52.1,
                        'gapout':99.5,
                        'PA':153.3,
                        'inc':19.9,
                        'dpc': 27.4},
            'HD206893':{'rin':34.8,
                        'rout':120.,
                        'gapin':49.,
                        'gapout':89.,
                        'PA':61.7,
                        'inc':40.,
                        'dpc': 40.7}}

fig, ax = display_image(subdirs=[sub], filts=['F1140C'], submodes=['ADI+RDI'], numKL=[10], 
                   ffov=24, cmaps_list=['gist_heat'], subVmin=[-1], subVmax=[2],
                   labelpos=[0.04, 0.02],
                    imtext_col='w', showKL=True, useticklabels=True, cbar_textoff=0.97,
                       hspace=0., wspace=0.2, stars=3, discdata=discdata)
plt.show()