import os
from astropy.io import fits
from spaceKLIP import database, coron1pipeline, coron2pipeline,coron3pipeline, imagetools, pyklippipeline

import numpy as np

class Sklip():
    def __init__(self, star, bkg, instr, subtraction, KLmodes, subsect,annuli, uncal_dir, 
                odir, SpT, photo_path, companions, stage='stage1', final_stage='end'):
        
        self.uncal_dir = uncal_dir
        self.odir = odir
        
        self.bkg = bkg

        self.stage = stage
        self.final_stage = final_stage

        self.KLmodes = KLmodes
        self.subsections = subsect
        self.annuli = annuli
        self.subtraction = subtraction
        
        self.star = star
        self.SpT = SpT
        self.instr = instr
        self.ramp_cores = 'all'

        self.Database = None
        self.AnalysisTools = None

        # self.load_database() # initialise database

    def load_database(self):
        ''' Loads data into database based on the processing stage
        different possible stages:
        - stage1
        - stage2
        - img_proc
        - sub
        Last 2 stages are addressed separately in rawcon() and calcon()
        - rawcon
        - calcon
        '''
        Database = database.Database(output_dir=self.odir)

        if self.stage == 'stage1' :
            fitsfiles = sorted([self.uncal_dir + f for f in os.listdir(self.uncal_dir) if f.endswith('_uncal.fits')])
            if self.instr == 'miri':
                bgfitsfiles =  sorted([f for f in fitsfiles if self.bkg[0] in f or self.bkg[1] in f]) # make note of which files are bkg
            elif self.instr == 'nircam': # NIRCam implementation
                bgfitsfiles = None
            Database.read_jwst_s012_data(datapaths=fitsfiles, bgpaths=bgfitsfiles)
            Database.summarize()
        elif self.stage == 'stage2':
            # new_idir = self.odir + 'stage1/'
            # fitsfiles = sorted([new_idir + f for f in os.listdir(new_idir) if f.endswith('_rateints.fits')])
            if self.instr == 'miri':
                new_idir = self.odir + 'stage1_masked/' # default masking groups for MIRI
                fitsfiles = sorted([new_idir + f for f in os.listdir(new_idir) if f.endswith('_rateints.fits')])
                bgfitsfiles =  sorted([f for f in fitsfiles if self.bkg[0] in f or self.bkg[1] in f]) # make note of which files are bkg
            elif self.instr == 'nircam': # NIRCam implementation
                new_idir = self.odir + 'stage1/'
                fitsfiles = sorted([new_idir + f for f in os.listdir(new_idir) if f.endswith('_rateints.fits')])
                bgfitsfiles = None
            Database.read_jwst_s012_data(datapaths=fitsfiles, bgpaths=bgfitsfiles)
            Database.summarize()
        elif self.stage == 'img_proc':
            new_idir = self.odir + 'stage2/'
            fitsfiles = sorted([new_idir + f for f in os.listdir(new_idir) if f.endswith('_calints.fits')])
            if self.instr == 'miri':
                bgfitsfiles =  sorted([f for f in fitsfiles if self.bkg[0] in f or self.bkg[1] in f]) # make note of which files are bkg
            elif self.instr == 'nircam': # NIRCam implementation
                bgfitsfiles = None
            Database.read_jwst_s012_data(datapaths=fitsfiles, bgpaths=bgfitsfiles)
            Database.summarize()
        elif self.stage == 'sub':
            # if self.instr == 'miri':
            #     # new_idir = self.odir + 'bgsub/' # TODO: check if image processing is done right -- Original 2023
            #     new_idir = self.odir + 'padded/' # new updated
            # elif self.instr == 'nircam':
            # previous version
            # new_idir = self.odir + 'padded/' # TODO: check if image processing is done right
            # new
            new_idir = self.odir + 'coadded/' # TODO: check if image processing is done right

            # could have option to manually start from whatever folder you want
                
            fitsfiles = sorted([new_idir + f for f in os.listdir(new_idir) if f.endswith('_calints.fits')])
            if self.instr == 'miri':
                bgfitsfiles =  sorted([f for f in fitsfiles if self.bkg[0] in f or self.bkg[1] in f]) # make note of which files are bkg
            elif self.instr == 'nircam': # NIRCam implementation
                bgfitsfiles = None
            Database.read_jwst_s012_data(datapaths=fitsfiles, bgpaths=bgfitsfiles)
            Database.summarize()
        # TODO: add loading database for rawcon and calcon -- requires knowing datapaths

        self.Database = Database

    def stage1_v0(self, subdir='stage1'):
        """ Ramp fitting stage -- Original 2023"""
        coron1pipeline.run_obs(database=self.Database,
                           steps={'saturation': {'n_pix_grow_sat': 1,
                                                 'grow_diagonal': False},
                                  'refpix': {'odd_even_columns': True,
                                             'odd_even_rows': True,
                                             'nlower': 0,
                                             'nupper': 0,
                                             'nleft': 0,
                                             'nright': 0,
                                             'nrow_off': 0,
                                             'ncol_off': 0},
                                  'dark_current': {'skip': True},
                                  'jump': {'rejection_threshold': 8.,
                                           'three_group_rejection_threshold': 8.,
                                           'four_group_rejection_threshold': 8.,
                                           'maximum_cores':self.ramp_cores},
                                  'ramp_fit': {'save_calibrated_ramp': False,
                                                'maximum_cores':self.ramp_cores}},
                           subdir=subdir)
    
    def stage1(self, subdir='stage1'): #changed subdir default
        """ Ramp fitting stage -- Adapted Original to 2024 version"""
        if self.instr == 'miri':
            jump_reject = 8.
        else:
            jump_reject = 4.
        coron1pipeline.run_obs(database=self.Database,
                           steps={'saturation': {'n_pix_grow_sat': 1,
                                                 'grow_diagonal': False},
                                  'refpix': {'odd_even_rows': True},
                                  'dark_current': {'skip': True},
                                  'ipc': {'skip': True},
                                  'jump': {'rejection_threshold': jump_reject,
                                           'maximum_cores':self.ramp_cores},
                                  'ramp_fit': {'save_calibrated_ramp': False,
                                                'maximum_cores':self.ramp_cores}},
                           subdir=subdir)
        
    def stage1_mask(self, groups_to_mask, subdir='stage1_masked'): #changed subdir default
        """ Ramp fitting stage -- miri_impr
        Includes masking of specific groups """
        coron1pipeline.run_obs(database=self.Database,
                           steps={'mask_groups': {'groups_to_mask': groups_to_mask,
                                                  'types':['REF', 'REF_BG']},
                                  'saturation': {'n_pix_grow_sat': 1,
                                                 'grow_diagonal': False},
                                  'refpix': {'odd_even_rows': True},
                                  'dark_current': {'skip': False},
                                #   'dark_current': {'skip': True},
                                  'ipc': {'skip': True},
                                  'jump': {'rejection_threshold': 8., 'maximum_cores':self.ramp_cores}, # ERS recommend threshold = 4 NIRCam, = 8 MIRI
                                  'ramp_fit': {'save_calibrated_ramp': False, 'maximum_cores':self.ramp_cores}},
                           subdir=subdir,
                           verbose = True)
    
    def stage2(self, subdir='stage2'):
        coron2pipeline.run_obs(database=self.Database,
                                steps={'outlier_detection': {'skip': False}},
                                subdir=subdir)
    
    def stage2_nircam(self, subdir='stage2'):
        coron2pipeline.run_obs(database=self.Database,
                           steps={'outlier_detection': {'skip': False}},
                           subdir=subdir)
        coron3pipeline.run_obs(database=self.Database,
                            steps={'klip': {'truncate': 100}},
                            subdir='stage3')

    def miri_image_processing_v0(self):
        ImageTools = imagetools.ImageTools(database=self.Database)

        # Remove the first frame due to reset switch charge delay. Only required for MIRI.
        ImageTools.remove_frames(index=[0],
                                #  frame=None, # old version had this param
                                types=['SCI', 'SCI_BG', 'REF', 'REF_BG'],
                                subdir='removed')

        # custom bad pixel map miri1140c
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
    
    def miri_image_processing(self):
        ImageTools = imagetools.ImageTools(database=self.Database)

        # Remove the first frame due to reset switch charge delay. Only required for MIRI.
        ImageTools.remove_frames(index=[0],
                                #  frame=None, # old version had this param
                                types=['SCI', 'SCI_BG', 'REF', 'REF_BG'],
                                subdir='removed')

        # custom bad pixel map miri1140c
        bpmap_F1140C = np.zeros((224, 288))
        bpmap_F1140C[115, 126] = 1
        bpmap_F1140C[111, 125] = 1
        bpmap_F1140C[95, 116] = 1
        # bpmap_F1140C[ 73, 42] = 1
        # bpmap_F1140C[ 218, 24] = 1
        # bpmap_F1140C[ 178, 14] = 1
        # bpmap_F1140C[ 207, 24] = 1
        # bpmap_F1140C[ 190, 36] = 1
        # bpmap_F1140C[ 160, 35] = 1
        # bpmap_F1140C[ 162, 46] = 1
        # bpmap_F1140C[ 174, 56] = 1
        # bpmap_F1140C[ 133, 67] = 1
        # bpmap_F1140C[ 208, 94] = 1
        # bpmap_F1140C[ 186, 101] = 1
        # bpmap_F1140C[ 178, 98] = 1
        # bpmap_F1140C[ 168, 106] = 1
        # bpmap_F1140C[ 161, 119] = 1
        # bpmap_F1140C[ 136, 115] = 1
        # bpmap_F1140C[ 113, 18] = 1
        # bpmap_F1140C[ 56, 56] = 1
        # bpmap_F1140C[ 73, 60] = 1
        # bpmap_F1140C[ 43, 58] = 1
        # bpmap_F1140C[ 40, 84] = 1
        # bpmap_F1140C[ 14, 16] = 1
        # bpmap_F1140C[ 35, 104] = 1
        # bpmap_F1140C[ 109, 103] = 1
        # bpmap_F1140C[ 141, 166] = 1
        # bpmap_F1140C[ 122, 174] = 1
        # bpmap_F1140C[ 140, 213] = 1
        # bpmap_F1140C[ 57, 58] = 1
        # bpmap_F1140C[ 103, 110] = 1
        # bpmap_F1140C[ 78, 127] = 1
        # bpmap_F1140C[ 67, 143] = 1
        # bpmap_F1140C[ 36, 153] = 1
        # bpmap_F1140C[ 79, 149] = 1
        # bpmap_F1140C[ 76, 151] = 1
        # bpmap_F1140C[ 75, 150] = 1
        # bpmap_F1140C[ 81, 167] = 1
        # bpmap_F1140C[ 40, 139] = 1
        # bpmap_F1140C[ 37, 154] = 1
        # bpmap_F1140C[ 29, 171] = 1
        # bpmap_F1140C[ 87, 186] = 1
        # bpmap_F1140C[ 29, 171] = 1
        # bpmap_F1140C[ 39, 195] = 1
        # bpmap_F1140C[ 41, 213] = 1
        # bpmap_F1140C[ 53, 215] = 1
        # bpmap_F1140C[ 53, 200] = 1
        # bpmap_F1140C[ 60, 73] = 1
        # bpmap_F1140C[ 191, 126] = 1
        # bpmap_F1140C[ 139, 24] = 1
        # bpmap_F1140C[ 97, 48] = 1
        # bpmap_F1140C[ 114, 48] = 1

        custom_kwargs = {}
        custom_kwargs['JWST_MIRI_MIRIMAGE_F1140C_NONE_4QPM_1140_MASK1140'] = bpmap_F1140C
        ImageTools.fix_bad_pixels(method='bpclean+custom+dqmed+timemed',
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
        ImageTools.crop_frames(npix=[13, 60, 7, 7],
                            #    npix=[13, 61, 7, 8],
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
        ImageTools.pad_frames(npix=80,
                            cval=np.nan,
                            types=['SCI', 'SCI_BG', 'REF', 'REF_BG'],
                            subdir='padded')
        
        # Coadd frames.
        # ImageTools.coadd_frames(nframes=1,
        #                         types=['SCI', 'REF'],
        #                         subdir='coadded')

    def miri_image_processing_godoy(self):
        ImageTools = imagetools.ImageTools(database=self.Database)

        # Coadd frames.
        # ImageTools.coadd_frames(nframes=1,
        #                         types=['SCI', 'REF'],
        #                         subdir='coadded')

        # TODO Currently only custom for HD107146 - should just be wavelength depednent but could update the custom list
        bpmap_F1140C = np.zeros((224, 288))
        bpmap_F1140C[115, 126] = 1
        bpmap_F1140C[111, 125] = 1
        bpmap_F1140C[95, 116] = 1
        custom_kwargs = {}
        custom_kwargs['JWST_MIRI_MIRIMAGE_F1140C_NONE_4QPM_1140_MASK1140'] = bpmap_F1140C
        ImageTools.fix_bad_pixels(method='bpclean+custom+dqmed+timemed', # new bpclean method
                                bpclean_kwargs={'sigclip': 5,
                                                'shift_x': [-1, 0, 1],
                                                'shift_y': [-1, 0, 1]},
                                custom_kwargs=custom_kwargs,
                                timemed_kwargs={},
                                dqmed_kwargs={'shift_x': [-1, 0, 1],
                                                'shift_y': [-1, 0, 1]},
                                medfilt_kwargs={'size': 4},
                                subdir='bpcleaned')
        
        # # Crop all frames.
        ImageTools.crop_frames(npix=[13, 60, 7, 7], # Numbers have changed. Why? TODO
                            types=['SCI', 'SCI_BG', 'REF', 'REF_BG'],
                            subdir='cropped')

        # Replace all nans. Same has not changed
        ImageTools.replace_nans(cval=0.,
                                types=['SCI', 'SCI_BG', 'REF', 'REF_BG'],
                                subdir='nanreplaced')
        
        # New bkg subtraction method
        ImageTools.subtract_background_godoy(subdir='bgsub_godoy')

        # Align frames
        # ImageTools.align_frames(method='fourier',
        #                         align_algo='leastsq',
        #                         kwargs={},
        #                         subdir='aligned')

        # Pad all frames.
        # ImageTools.pad_frames(npix=80,
        #                     cval=np.nan,
        #                     types=['SCI', 'SCI_BG', 'REF', 'REF_BG'],
        #                     subdir='padded')

        # Coadd frames.
        ImageTools.coadd_frames(nframes=None, #None
                                types=['SCI', 'REF'],
                                subdir='coadded')
    
    def nircam_image_processing(self):
        ImageTools = imagetools.ImageTools(database=self.Database)

        ImageTools.subtract_median(types=['SCI', 'SCI_TA', 'SCI_BG',
                                        'REF', 'REF_TA', 'REF_BG'],
                                subdir='medsub')
        
        bpmap = np.zeros((320, 320))
        bpmap[90, 122] = 1
        bpmap[120, 90] = 1
        bpmap[121, 90] = 1
        bpmap[122, 90] = 1
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
        custom_kwargs['JWST_NIRCAM_NRCA2_F200W_MASKRND_MASKA335R_SUB320A335R'] = bpmap
        custom_kwargs['JWST_NIRCAM_NRCALONG_F250M_MASKRND_MASKA335R_SUB320A335R'] = bpmap
        custom_kwargs['JWST_NIRCAM_NRCALONG_F300M_MASKRND_MASKA335R_SUB320A335R'] = bpmap
        custom_kwargs['JWST_NIRCAM_NRCALONG_F356W_MASKRND_MASKA335R_SUB320A335R'] = bpmap
        custom_kwargs['JWST_NIRCAM_NRCALONG_F410M_MASKRND_MASKA335R_SUB320A335R'] = bpmap
        custom_kwargs['JWST_NIRCAM_NRCALONG_F444W_MASKRND_MASKA335R_SUB320A335R'] = bpmap
        # ImageTools.fix_bad_pixels(method='bpclean+timemed+dqmed+medfilt',
        ImageTools.fix_bad_pixels(method='timemed+dqmed+medfilt',
                              bpclean_kwargs={'sigclip': 3,
                                              'shift_x': [-1, 0, 1],
                                              'shift_y': [-1, 0, 1]},
                              custom_kwargs=custom_kwargs,
                              timemed_kwargs={},
                              dqmed_kwargs={'shift_x': [-1, 0, 1],
                                            'shift_y': [-1, 0, 1]},
                              medfilt_kwargs={'size': 4},
                              subdir='bpcleaned_test')
        # return
        
        # Replace nans.
        ImageTools.replace_nans(cval=0.,
                            types=['SCI', 'SCI_BG', 'REF', 'REF_BG'],
                            subdir='nanreplaced')
        
        ImageTools.align_frames(method='fourier',
                            kwargs={},
                            subdir='aligned')
        
        ImageTools.coadd_frames(nframes=None, 
                                types=['SCI', 'REF'], 
                                subdir='coadded')
        
        # # Pad all frames.
        # ImageTools.pad_frames(npix=160,
        #                   cval=np.nan,
        #                   types=['SCI', 'SCI_BG', 'REF', 'REF_BG'],
        #                   subdir='padded')

    def nircam_image_processing_v0(self):
        ImageTools = imagetools.ImageTools(database=self.Database)

        ImageTools.subtract_median(types=['SCI', 'SCI_TA', 'SCI_BG',
                                        'REF', 'REF_TA', 'REF_BG'],
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
        custom_kwargs['JWST_NIRCAM_NRCA2_F200W_MASKRND_MASKA335R_SUB320A335R'] = bpmap
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
        

        ImageTools.replace_nans(cval=0.,
                                types=['SCI', 'SCI_BG', 'REF', 'REF_BG'],
                                subdir='nanreplaced')
        
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
        
        # Recenter frames. Before, update NIRCam coronagraphic mask centers to
        # on-sky values measured by Jarron. Might not be required for simulated
        # data!
        ImageTools.update_nircam_centers()
        ImageTools.recenter_frames(method='fourier',
                                subpix_first_sci_only=False,
                                spectral_type= self.SpT,
                                kwargs={},
                                subdir='recentered')
        
        # Use image registration to align all frames in concatenation to first
        # science frame in that concatenation.
        ImageTools.align_frames(method='fourier',
                                kwargs={},
                                subdir='aligned')
        
        # Pad all frames.
        ImageTools.pad_frames(npix=160,
                            cval=np.nan,
                            types=['SCI', 'SCI_BG', 'REF', 'REF_BG'],
                            subdir='padded')

    def klip_subtraction(self, subdir='klipsub_coadded'):
        pyklippipeline.run_obs(database=self.Database,
                        kwargs={'mode': self.subtraction,
                                 'annuli': self.annuli,
                                 'subsections': self.subsections,
                                 'numbasis': self.KLmodes,
                                 'algo': 'klip',
                                 'save_rolls': True},
                        subdir=subdir)

    def run_pipeline_v0(self):
        while self.stage != self.final_stage:
            if self.stage == 'stage1':
                # TODO: add print statements to know what steps we are at
                self.load_database()
                self.stage1()
                self.stage = 'stage2'

            elif self.stage == 'stage2':
                print('\n', 'Performing Stage2 for {}'.format(self.star), '\n')
                self.load_database()
                self.stage2()
                self.stage = 'img_proc'

            elif self.stage == 'img_proc':
                print('\n', 'Performing Image Processing for {}'.format(self.star), '\n')
                self.load_database()
                if self.instr == 'miri':
                    self.miri_image_processing()
                elif self.instr == 'nircam':
                    self.nircam_image_processing()
                self.stage = 'sub'

            elif self.stage == 'sub':
                print('\n', 'Performing Image Processing for {}'.format(self.star), '\n')
                self.load_database()
                self.klip_subtraction()
                self.stage = 'rawcon'
            
            # if self.stage == 'rawcon':
                
    def SwitchRefSci(self, hdul_Ref,hdul_Sci,Cx,Cy,KernelPix,Method):
        SciCube=hdul_Sci["SCI"].data[0,-1]
        RefCube=hdul_Ref["SCI"].data[0,-1]
        SciCubeCrop=SciCube[(Cy-KernelPix):(Cy+KernelPix),(Cx-KernelPix):(Cx+KernelPix)]
        RefCubeCrop=RefCube[(Cy-KernelPix):(Cy+KernelPix),(Cx-KernelPix):(Cx+KernelPix)]
        if Method=="MaxPixel":
            
            if np.nanmax(RefCubeCrop)>np.nanmax(SciCubeCrop):
                return False
            else:
                return True
        elif "Summed" in Method:
            if np.nansum(RefCubeCrop)>np.nansum(SciCubeCrop):
                return False
            else:
                return True
    
    def optimise_groups(self, RefPath,SciPath,KernelPix=15,Method="Summed+Nan"):
        hdul_Ref=fits.open(RefPath)
        hdul_Sci=fits.open(SciPath)
        integration_Ref=0 #Just chose the first integration. THis could be changed?
        integration_Sci=0
        
        group_Ref=int(hdul_Ref["SCI"].header["NAXIS3"]) #How many total groups are there?
        group_Sci=int(hdul_Sci["SCI"].header["NAXIS3"])
        print(f"the total groups available are: {group_Ref}")
        Cx=int(hdul_Ref["SCI"].header["CRPIX1"]) #Centre pix location for image (centred on PSF not frame)
        Cy=int(hdul_Ref["SCI"].header["CRPIX2"])

        if self.SwitchRefSci(hdul_Ref,hdul_Sci,Cx,Cy,KernelPix,Method):
            print("The Science target is brighter than the reference images!")
            group_Ref,group_Sci=group_Sci,group_Ref
            hdul_Ref,hdul_Sci=hdul_Sci,hdul_Ref
            
        
        SciCube=hdul_Sci["SCI"].data[integration_Sci,group_Sci-1] #Load the Science Cube
        
        SciCubeCrop=SciCube[(Cy-KernelPix):(Cy+KernelPix),(Cx-KernelPix):(Cx+KernelPix)] #Crop it around the central PSF
        if "Nan" in Method:
            SciCubeCrop=np.where(SciCubeCrop>np.max(SciCubeCrop)//2,SciCubeCrop,np.nan) #Get rid of data that isnt at least half the max pixel
                                                                                    #this is to isolate the PSF better?
        Nminimized=np.inf #Large number used for minimisation.
        
        for i in range(group_Ref): 
            RefCube=hdul_Ref["SCI"].data[integration_Ref,i] #load cube
            RefCubeCrop=RefCube[(Cy-KernelPix):(Cy+KernelPix),(Cx-KernelPix):(Cx+KernelPix)] #crop it
            if "Summed" in Method:
                if "Nan" in Method:
                    RefCubeCrop=np.where(RefCubeCrop>np.max(RefCubeCrop)//2,RefCubeCrop,np.nan) #delete unrelated data
                SubtractedCube=SciCubeCrop-RefCubeCrop #take reference layer from Science layer
                summed=np.nansum(SubtractedCube) #add up total flux
                if summed<Nminimized and summed>0: #minimize total flux. making sure total flux is still above 0 (could have over subtraction issues still)
                    Nminimized=summed
                    minimizedGroups=i
                    # minimized=SubtractedCube
            elif Method=="MaxPixel":
                if np.nanmax(SciCubeCrop)<np.nanmax(RefCubeCrop): #when any pixel's counts is larger than any pixels count, then return that image
                    print(f"Optimized Groups using {Method} is: {i}")
                    return np.arange(i, group_Ref)
        print(f"Optimized Groups using {Method} is: {minimizedGroups}")
        return np.arange(minimizedGroups, group_Ref) #return frame that minimizes the total flux
        
    
    def run_pipeline(self, godoy=False):
        while self.stage != self.final_stage:
            if self.stage == 'stage1':
                print('\n', 'Performing Stage1 for {}'.format(self.star), '\n')
                self.load_database()
                if godoy:
                    data = self.Database
                    for i, key in enumerate(data.obs.keys()): # get the path of a reference and science target
                        # print(data.obs[key].colnames)
                        science = data.obs[key][data.obs[key]['TYPE'] == 'SCI']['FITSFILE'][0]
                        ref = data.obs[key][data.obs[key]['TYPE'] == 'REF']['FITSFILE'][0]
                    
                    groups_to_mask = self.optimise_groups(RefPath=ref, SciPath= science)
                    # groups_to_mask = np.arange(ref_groups // 8, ref_groups) ## TODO optimise the number of groups to consider
                    self.stage1_mask(groups_to_mask=groups_to_mask)
                else:
                    self.stage1()
                self.stage = 'stage2'

            elif self.stage == 'stage2':
                print('\n', 'Performing Stage2 for {}'.format(self.star), '\n')
                self.load_database()
                if self.instr == 'nircam':
                    self.stage2_nircam()
                else:
                    self.stage2()
                self.stage = 'img_proc'

            elif self.stage == 'img_proc':
                print('\n', 'Performing Image Processing for {}'.format(self.star), '\n')
                self.load_database()
                if self.instr == 'miri':
                    if godoy:
                        self.miri_image_processing_godoy()
                    else:
                        self.miri_image_processing()
                elif self.instr == 'nircam':
                    self.nircam_image_processing()
                self.stage = 'sub'

            elif self.stage == 'sub':
                print('\n', 'Performing Image Processing for {}'.format(self.star), '\n')
                self.load_database()
                self.klip_subtraction()
                self.stage = 'rawcon'
            
            # if self.stage == 'rawcon':
                
            
