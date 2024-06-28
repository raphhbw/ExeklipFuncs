import os
import pdb
import sys

import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np

from spaceKLIP import database, imagetools, coron1pipeline, coron2pipeline,\
    coron3pipeline, pyklippipeline, classpsfsubpipeline, analysistools

if __name__ == '__main__':
    star = 'HD195627'
    odir = '/Users/rbendahan/Astronomy/data/JWST/NIRCam-200/{}/spaceklip/'.format(star)

    Database = database.Database(output_dir=odir)
    # datapaths = [odir + 'klipsub/ADI+RDI_NANNU6_NSUBS6_JWST_MIRI_MIRIMAGE_F1140C_NONE_4QPM_1140_MASK1140-KLmodes-all.fits']
    datapaths = [odir + 'klipsub_v1/RDI_NANNU1_NSUBS1_JWST_NIRCAM_NRCA2_F200W_MASKRND_MASKA335R_SUB320A335R-KLmodes-all.fits']
    Database.read_jwst_s3_data(datapaths)

    # Initialize spaceKLIP analysis tools class.
    AnalysisTools = analysistools.AnalysisTools(Database)

    # Compute raw contrast.
    AnalysisTools.raw_contrast(starfile='/Users/rbendahan/Astronomy/data/photometry/hd195627_phoenix_sol+modbb_disk_r_.txt',
                               spectral_type='F0V',
                            #    companions=[[0.431, -0.717, 2.]],
                               subdir='rawcon_v1')

    AnalysisTools.calibrate_contrast(rawcon_subdir='rawcon_v1',
                                    #  companions=[[0.431, -0.717, 2.]],
                                     use_saved=False )