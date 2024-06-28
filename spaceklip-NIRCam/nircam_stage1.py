import os
import pdb
import sys

import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np

from spaceKLIP import database, imagetools, coron1pipeline, coron2pipeline,\
    coron3pipeline, pyklippipeline, classpsfsubpipeline, analysistools

# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":

    # Set input and output directories and get uncal files.
    idir = '/Users/rbendahan/Astronomy/data/JWST/NIRCam-200/HD195627/uncal/'
    odir = '/Users/rbendahan/Astronomy/data/JWST/NIRCam-200/HD195627/spaceklip/'
    
    fitsfiles = sorted([idir + f for f in os.listdir(idir) if f.endswith('_uncal.fits')])

    # Initialize spaceKLIP database and read uncal files.
    Database = database.Database(output_dir=odir)
    Database.read_jwst_s012_data(datapaths=fitsfiles,
                                 bgpaths=None)
    
    Database.summarize()
    # sys.exit()

    # v1
    coron1pipeline.run_obs(database=Database,
                           steps={'saturation': {'n_pix_grow_sat': 1,
                                                 'grow_diagonal': False},
                                  'refpix': {'odd_even_columns': True,
                                             'odd_even_rows': True,
                                             'nlower': 4,
                                             'nupper': 4,
                                             'nleft': 4,
                                             'nright': 4,
                                             'nrow_off': 0,
                                             'ncol_off': 0},
                                  'dark_current': {'skip': True},
                                  'jump': {'rejection_threshold': 4.,
                                           'three_group_rejection_threshold': 4.,
                                           'four_group_rejection_threshold': 4.},
                                  'ramp_fit': {'save_calibrated_ramp': False}},
                           subdir='stage1_v1')


    # v2
    # coron1pipeline.run_obs(database=Database,
    #                        steps={'saturation': {'n_pix_grow_sat': 1,
    #                                              'grow_diagonal': True},
    #                               'refpix': {'odd_even_columns': True,
    #                                          'odd_even_rows': True,
    #                                          'nlower': 4,
    #                                          'nupper': 4,
    #                                          'nleft': 4,
    #                                          'nright': 4,
    #                                          'nrow_off': 0,
    #                                          'ncol_off': 0},
    #                               'dark_current': {'skip': False},
    #                               'jump': {'rejection_threshold': 4.,
    #                                        'three_group_rejection_threshold': 4.,
    #                                        'four_group_rejection_threshold': 4.},
    #                               'ramp_fit': {'save_calibrated_ramp': False}},
    #                        subdir='stage1_v2')

    # v3
    # coron1pipeline.run_obs(database=Database,
    #                        steps={'saturation': {'n_pix_grow_sat': 1,
    #                                              'grow_diagonal': False},
    #                               'refpix': {'odd_even_columns': True,
    #                                          'odd_even_rows': True,
    #                                          'nlower': 4,
    #                                          'nupper': 4,
    #                                          'nleft': 4,
    #                                          'nright': 4,
    #                                          'nrow_off': 0,
    #                                          'ncol_off': 0},
    #                               'dark_current': {'skip': False},
    #                               'jump': {'rejection_threshold': 4.,
    #                                        'three_group_rejection_threshold': 4.,
    #                                        'four_group_rejection_threshold': 4.},
    #                               'ramp_fit': {'save_calibrated_ramp': False}},
    #                        subdir='stage1_v3')

    # normal
    # coron1pipeline.run_obs(database=Database,
    #                        steps={'saturation': {'n_pix_grow_sat': 1,
    #                                              'grow_diagonal': True},
    #                               'refpix': {'odd_even_rows': True},
    #                               'dark_current': {'skip': False},
    #                               'jump': {'rejection_threshold': 8., 
    #                                        'find_showers':False,
    #                                   'expand_large_events':False,
    #                                   'after_jump_flag_time1':10,
    #                                   'after_jump_flag_time2':10,
    #                                   'maximum_cores':'half'},
    #                               'ramp_fit': {'save_calibrated_ramp': False, 
    #                                            'maximum_cores':'half'}},
    #                        subdir='stage1')