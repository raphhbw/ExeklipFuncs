# =============================================================================
# IMPORTS
# =============================================================================

import os
import pdb
import sys

import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np

from spaceKLIP import database, pyklippipeline


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    
    # Set the input and output directories and grab the input FITS files.
    # input_dir = '/Users/rbendahan/Astronomy/data/JWST/NIRCam-444/HD195627/spaceklip/stage1_v3/'
    input_dir = '/Users/rbendahan/Astronomy/data/JWST/NIRCam-200/HD195627/spaceklip/stage1/'
    output_dir = '/Users/rbendahan/Astronomy/data/JWST/NIRCam-200/HD195627/spaceklip/'

    idir = input_dir.replace('stage1', 'padded')
    # idir = input_dir.replace('stage1_v3', 'padded_v3')
    Database = database.Database(output_dir=output_dir)
    fitsfiles = sorted([idir + f for f in os.listdir(idir) if f.endswith('_calints.fits')])
    Database.read_jwst_s012_data(datapaths=fitsfiles, bgpaths=None)
    
    pyklippipeline.run_obs(database=Database,
                        kwargs={'mode': ['RDI'],
                                 'annuli': [1,5],
                              # annuli & subsections defines the different sections where 
                              # klip does subtraction
                                 'subsections': [1],
                                 'numbasis': [1,2,3,4,5,10,20, 50, 100, 200],
                              # numbasis = klmode
                                 'algo': 'klip',
                                 'save_rolls': True},
                        subdir='klipsub')
                        # subdir='klipsub_v3')