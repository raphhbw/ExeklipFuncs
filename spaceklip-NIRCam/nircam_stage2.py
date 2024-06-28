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
    
    # input_dir = '/Users/rbendahan/Astronomy/data/JWST/NIRCam-444/HD195627/spaceklip/stage1_v3/'
    input_dir = '/Users/rbendahan/Astronomy/data/JWST/NIRCam-200/HD195627/spaceklip/stage1/'
    output_dir = '/Users/rbendahan/Astronomy/data/JWST/NIRCam-200/HD195627/spaceklip/'

    skip_stage2 = False
    SpT = 'F0V'

    if skip_stage2:
        # idir = input_dir.replace('stage1_v3', 'stage2_v3')
        idir = input_dir.replace('stage1', 'stage2')
        Database = database.Database(output_dir=output_dir)
        fitsfiles = sorted([idir + f for f in os.listdir(idir) if f.endswith('_calints.fits')])
        Database.read_jwst_s012_data(datapaths=fitsfiles, bgpaths=None)
        Database.summarize()

    else:
        Database = database.Database(output_dir=output_dir)
        fitsfiles = sorted([input_dir + f for f in os.listdir(input_dir) if f.endswith('_rateints.fits')])
        Database.read_jwst_s012_data(datapaths=fitsfiles, bgpaths=None)
        Database.summarize()
        coron2pipeline.run_obs(database=Database,
                                steps={'outlier_detection': {'skip': False}},
                                # subdir='stage2_v3')
                                subdir='stage2')
        
    ImageTools = imagetools.ImageTools(database=Database)

    ImageTools.subtract_median(types=['SCI', 'SCI_TA', 'SCI_BG',
                                      'REF', 'REF_TA', 'REF_BG'],
                            #    subdir='medsub_v2')
                               subdir='medsub')
    
    bpmap = np.zeros((320, 320))
    bpmap[153, 153] = 1
    bpmap[159, 142] = 1
    bpmap[159, 143] = 1
    bpmap[159, 157] = 1
    bpmap[159, 158] = 1
    bpmap[159, 164] = 1
    bpmap[164, 141] = 1
    bpmap[165, 146] = 1
    bpmap[169, 140] = 1
    bpmap[177, 120] = 1
    custom_kwargs = {}
    custom_kwargs['JWST_NIRCAM_NRCA2_F200W_MASKRND_MASKA335R_SUB320A335R'] = bpmap # added
    custom_kwargs['JWST_NIRCAM_NRCALONG_F250M_MASKRND_MASKA335R_SUB320A335R'] = bpmap
    custom_kwargs['JWST_NIRCAM_NRCALONG_F300M_MASKRND_MASKA335R_SUB320A335R'] = bpmap
    custom_kwargs['JWST_NIRCAM_NRCALONG_F356W_MASKRND_MASKA335R_SUB320A335R'] = bpmap
    custom_kwargs['JWST_NIRCAM_NRCALONG_F410M_MASKRND_MASKA335R_SUB320A335R'] = bpmap
    custom_kwargs['JWST_NIRCAM_NRCALONG_F444W_MASKRND_MASKA335R_SUB320A335R'] = bpmap
    ImageTools.fix_bad_pixels(method='bpclean+custom+timemed+dqmed+medfilt',
                              bpclean_kwargs={'sigclip': 5,
                                              'shift_x': [-1, 0, 1],
                                              'shift_y': [-1, 0, 1]},
                              custom_kwargs=custom_kwargs,
                              timemed_kwargs={},
                              dqmed_kwargs={'shift_x': [-1, 0, 1],
                                            'shift_y': [-1, 0, 1]},
                              medfilt_kwargs={'size': 4},
                              subdir='bpcleaned')
                            #   subdir='bpcleaned_v3')
    

    ImageTools.replace_nans(cval=0.,
                            types=['SCI', 'SCI_BG', 'REF', 'REF_BG'],
                            subdir='nanreplaced')
                            # subdir='nanreplaced_v3')
    
    # Blur frames.
    fact = {}
    fact['JWST_NIRCAM_NRCA2_F200W_MASKRND_MASKA335R_SUB320A335R'] = ['auto'] * 11 # added
    fact['JWST_NIRCAM_NRCALONG_F250M_MASKRND_MASKA335R_SUB320A335R'] = ['auto'] * 11
    fact['JWST_NIRCAM_NRCALONG_F300M_MASKRND_MASKA335R_SUB320A335R'] = [None] * 11
    fact['JWST_NIRCAM_NRCALONG_F356W_MASKRND_MASKA335R_SUB320A335R'] = [None] * 11
    fact['JWST_NIRCAM_NRCALONG_F410M_MASKRND_MASKA335R_SUB320A335R'] = [None] * 11
    fact['JWST_NIRCAM_NRCALONG_F444W_MASKRND_MASKA335R_SUB320A335R'] = [None] * 11
    ImageTools.blur_frames(fact=fact,
                           types=['SCI', 'SCI_BG', 'REF', 'REF_BG'],
                           subdir='blurred')
                        #    subdir='blurred_v3')
    
    # Recenter frames. Before, update NIRCam coronagraphic mask centers to
    # on-sky values measured by Jarron. Might not be required for simulated
    # data!
    ImageTools.update_nircam_centers()
    ImageTools.recenter_frames(method='fourier',
                               subpix_first_sci_only=False,
                               spectral_type= SpT,
                               kwargs={},
                               subdir='recentered')
                            #    subdir='recentered_v3')
    
    # Use image registration to align all frames in concatenation to first
    # science frame in that concatenation.
    ImageTools.align_frames(method='fourier',
                            kwargs={},
                            subdir='aligned')
                            # subdir='aligned_v3')
    
    # Pad all frames.
    ImageTools.pad_frames(npix=160,
                          cval=np.nan,
                          types=['SCI', 'SCI_BG', 'REF', 'REF_BG'],
                          subdir='padded')
                        #   subdir='padded_v3')
    
