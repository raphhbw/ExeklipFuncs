# =============================================================================
# IMPORTS
# =============================================================================

import os

import numpy as np

from spaceKLIP import database, coron2pipeline, imagetools


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    
    # Set the input and output directories and grab the input FITS files.
    # star = 'HD21997'

    os.nice(10) # change niceness

    idir = '/data/rb941/JWST/MIRI-1140/HD92945/spaceklip/stage1/'
    odir = '/data/rb941/JWST/MIRI-1140/HD92945/spaceklip/'

    skip_stage2 = False

    if skip_stage2:
        idir = idir.replace('stage1', 'stage2')
        Database = database.Database(output_dir=odir)
        fitsfiles = sorted([idir + f for f in os.listdir(idir) if f.endswith('_calints.fits')])
        # bgfitsfiles =  sorted([f for f in fitsfiles if '8011' in f or '8015' in f]) # HD206893
        # bgfitsfiles =  sorted([f for f in fitsfiles if '8001' in f or '8005' in f]) # HD107146
        bgfitsfiles =  sorted([f for f in fitsfiles if '8006' in f or '8010' in f]) # HD92945
        # bgfitsfiles =  sorted([f for f in fitsfiles if '8011' in f or '8011' in f]) # HD21997
        Database.read_jwst_s012_data(datapaths=fitsfiles,
                                #  psflibpaths=None,
                                 bgpaths=bgfitsfiles)
    else:
        Database = database.Database(output_dir=odir)
        fitsfiles = sorted([idir + f for f in os.listdir(idir) if f.endswith('_rateints.fits')])
        print(fitsfiles, len(fitsfiles))
        # bgfitsfiles = sorted([f for f in fitsfiles if '8011' in f or '8015' in f])# HD206893
        # bgfitsfiles =  sorted([f for f in fitsfiles if '8001' in f or '8005' in f]) # HD107146
        bgfitsfiles =  sorted([f for f in fitsfiles if '8006' in f or '8010' in f]) # HD92945
        # bgfitsfiles =  sorted([f for f in fitsfiles if '8011' in f or '8011' in f]) # HD21997
        print(bgfitsfiles)
        Database.read_jwst_s012_data(datapaths=fitsfiles,
                                #  psflibpaths=None,
                                 bgpaths=bgfitsfiles)
        Database.summarize()
        coron2pipeline.run_obs(database=Database,
                                steps={'outlier_detection': {'skip': False}},
                                subdir='stage2')
        
    ImageTools = imagetools.ImageTools(database=Database)

    # Remove the first frame due to reset switch charge delay. Only required
    # for MIRI.
    ImageTools.remove_frames(index=[0],
                            #  frame=None, # old version had this param
                             types=['SCI', 'SCI_BG', 'REF', 'REF_BG'],
                             subdir='removed')
    
    # Fix bad pixels using custom spaceKLIP routines. Multiple routines can be
    # combined in custom order by joining them with + signs.
    # - bpclean: use sigma clipping to find additional bad pixels.
    # - custom: use custom bad pixel map.
    # - timemed: replace pixels which are only bad in some frames with their
    #            median value from good frames.
    # - dqmed:   replace bad pixels with median of surrounding good pixels.
    # - medfilt: replace bad pixels with image plane median filter.
    bpmap_F1140C = np.zeros((224, 288))
    bpmap_F1140C[ 73, 42] = 1
    bpmap_F1140C[ 218, 24] = 1
    bpmap_F1140C[ 178, 14] = 1
    bpmap_F1140C[ 207, 24] = 1
    bpmap_F1140C[ 190, 36] = 1
    bpmap_F1140C[ 160, 35] = 1
    bpmap_F1140C[ 162, 46] = 1
    bpmap_F1140C[ 174, 56] = 1
    bpmap_F1140C[ 133, 67] = 1
    bpmap_F1140C[ 208, 94] = 1
    bpmap_F1140C[ 186, 101] = 1
    bpmap_F1140C[ 178, 98] = 1
    bpmap_F1140C[ 168, 106] = 1
    bpmap_F1140C[ 161, 119] = 1
    bpmap_F1140C[ 136, 115] = 1
    bpmap_F1140C[ 113, 18] = 1
    bpmap_F1140C[ 56, 56] = 1
    bpmap_F1140C[ 73, 60] = 1
    bpmap_F1140C[ 43, 58] = 1
    bpmap_F1140C[ 40, 84] = 1
    bpmap_F1140C[ 14, 16] = 1
    bpmap_F1140C[ 35, 104] = 1
    bpmap_F1140C[ 109, 103] = 1
    bpmap_F1140C[ 141, 166] = 1
    bpmap_F1140C[ 122, 174] = 1
    bpmap_F1140C[ 140, 213] = 1
    bpmap_F1140C[ 57, 58] = 1
    bpmap_F1140C[ 103, 110] = 1
    bpmap_F1140C[ 78, 127] = 1
    bpmap_F1140C[ 67, 143] = 1
    bpmap_F1140C[ 36, 153] = 1
    bpmap_F1140C[ 79, 149] = 1
    bpmap_F1140C[ 76, 151] = 1
    bpmap_F1140C[ 75, 150] = 1
    bpmap_F1140C[ 81, 167] = 1
    bpmap_F1140C[ 40, 139] = 1
    bpmap_F1140C[ 37, 154] = 1
    bpmap_F1140C[ 29, 171] = 1
    bpmap_F1140C[ 87, 186] = 1
    bpmap_F1140C[ 29, 171] = 1
    bpmap_F1140C[ 39, 195] = 1
    bpmap_F1140C[ 41, 213] = 1
    bpmap_F1140C[ 53, 215] = 1
    bpmap_F1140C[ 53, 200] = 1
    bpmap_F1140C[ 60, 73] = 1
    bpmap_F1140C[ 191, 126] = 1
    bpmap_F1140C[ 139, 24] = 1
    bpmap_F1140C[ 97, 48] = 1
    bpmap_F1140C[ 114, 48] = 1

    custom_kwargs = {}
    custom_kwargs['JWST_MIRI_MIRIMAGE_F1140C_NONE_4QPM_1140_MASK1140'] = bpmap_F1140C
    ImageTools.fix_bad_pixels(method='custom+timemed+dqmed',
                              bpclean_kwargs={'sigclip': 5,
                                              'shift_x': [-1, 0, 1],
                                              'shift_y': [-1, 0, 1]},
                              custom_kwargs=custom_kwargs,
                              timemed_kwargs={},
                              dqmed_kwargs={'shift_x': [-1, 0, 1],
                                            'shift_y': [-1, 0, 1]},
                              medfilt_kwargs={'size': 4},
                              subdir='bpcleaned')
    
    # Crop all frames.
    # npix defined manually
    ImageTools.crop_frames(npix=[13, 61, 7, 8],
                           types=['SCI', 'SCI_BG', 'REF', 'REF_BG'],
                           subdir='cropped')
    
    # Replace all nans.
    ImageTools.replace_nans(cval=0.,
                            types=['SCI', 'SCI_BG', 'REF', 'REF_BG'],
                            subdir='nanreplaced')
    
    # # Perform background subtraction to remove MIRI glowstick. Only required for MIRI.
    # ImageTools.subtract_background(nsplit=1, old calcon branch
    ImageTools.subtract_background(subdir='bgsub')
    
    # Align frames
    # ImageTools.align_frames(method='fourier',
    #                         align_algo='leastsq',
    #                         kwargs={},
    #                         subdir='aligned')

    # Pad all frames.
    ImageTools.pad_frames(npix=160,
                          cval=np.nan,
                          types=['SCI', 'SCI_BG', 'REF', 'REF_BG'],
                          subdir='padded')
    
    # Coadd frames.
    # ImageTools.coadd_frames(nframes=1,
    #                         types=['SCI', 'REF'],
    #                         subdir='coadded')