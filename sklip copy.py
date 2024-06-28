import os, glob, sys
from astropy.io import fits
from spaceKLIP import database, coron1pipeline, coron2pipeline,coron3pipeline, imagetools, pyklippipeline, analysistools

import numpy as np

class Sklip():
    def __init__(self, star, bkg, instr, subtraction, KLmodes, subsect,annuli, uncal_dir, 
                odir, SpT, photo_path, companions, stage='stage1', final_stage='end', skip_1_f = False):
        
        self.uncal_dir = uncal_dir
        self.odir = odir
        
        self.bkg = bkg

        self.stage = stage
        self.final_stage = final_stage
        self.stage1_subdir = 'stage1'

        self.KLmodes = KLmodes
        self.subsections = subsect
        self.annuli = annuli
        self.subtraction = subtraction
        self.photo_path = photo_path
        self.companions = companions
        self.contrast_path = None # default to None

        self.godoy = False # default to False
        self.skip_1_f = skip_1_f # default to False
        
        self.star = star
        self.SpT = SpT
        self.instr = instr
        self.ramp_cores = 'all'

        self.Database = None
        self.AnalysisTools = None

        # self.load_database() # initialise database

    def load_database(self, stage1_path='stage1', subpath='padded'):
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
            new_idir = self.odir + '{}/'.format(stage1_path)
            fitsfiles = sorted([new_idir + f for f in os.listdir(new_idir) if f.endswith('_rateints.fits')])
            if self.instr == 'miri':
                # fitsfiles = sorted([new_idir + f for f in os.listdir(new_idir) if f.endswith('_rateints.fits')])
                bgfitsfiles =  sorted([f for f in fitsfiles if self.bkg[0] in f or self.bkg[1] in f]) # make note of which files are bkg
            elif self.instr == 'nircam': # NIRCam implementation
                # fitsfiles = sorted([new_idir + f for f in os.listdir(new_idir) if f.endswith('_rateints.fits')])
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
            new_idir = self.odir + '{}/'.format(subpath)
                
            fitsfiles = sorted([new_idir + f for f in os.listdir(new_idir) if f.endswith('_calints.fits')])
            if self.instr == 'miri':
                bgfitsfiles =  sorted([f for f in fitsfiles if self.bkg[0] in f or self.bkg[1] in f]) # make note of which files are bkg
            elif self.instr == 'nircam': # NIRCam implementation
                bgfitsfiles = None
            Database.read_jwst_s012_data(datapaths=fitsfiles, bgpaths=bgfitsfiles)
            Database.summarize()
        
        elif (self.stage == 'rawcon') or (self.stage == 'calcon'):
            Database.read_jwst_s3_data(self.contrast_path)
            Database.summarize()

        self.Database = Database
    
    def stage1(self, subdir='stage1'): #changed subdir default
        """ Ramp fitting stage -- Adapted Original to 2024 version"""
        if self.instr == 'miri':
            jump_reject = 8.
        else:
            jump_reject = 4.
            coron1pipeline.run_obs(database=self.Database,
                                skip_fnoise=self.skip_1_f,
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
                                    'persistence': {'skip': True},
                                    'jump': {'rejection_threshold': jump_reject,
                                       'three_group_rejection_threshold': 4.,
                                       'four_group_rejection_threshold': 4.,
                                            'maximum_cores':self.ramp_cores},
                                    'ramp_fit': {'save_calibrated_ramp': False,
                                                    'maximum_cores':self.ramp_cores}},
                            subdir=subdir)
        
    def stage1_mask(self, subdir='stage1_masked'): #changed subdir default
        """ Ramp fitting stage -- miri_impr
        Includes masking of specific groups """
        coron1pipeline.run_obs(database=self.Database,
                   steps={'mask_groups': {'skip': False,
                                          'mask_method': 'advanced', # new method from aarynn
                                          'save_results':False, 
                                          'types':['REF', 'REF_BG']},
                          'experimental_jumpramp': {'use': True, 'nproc':4}},
                   subdir=subdir,
                   verbose=True)
        
    def stage1_mask_v0(self, groups_to_mask, jump_threshold, subdir='stage1_masked'): #changed subdir default
        """ Ramp fitting stage -- miri_impr
        Includes masking of specific groups """
        print(type(jump_threshold))
        coron1pipeline.run_obs(database=self.Database,
                           steps={'mask_groups': {'groups_to_mask': groups_to_mask,
                                                  'types':['REF', 'REF_BG']},
                                  'saturation': {'n_pix_grow_sat': 1,
                                                 'grow_diagonal': False},
                                  'refpix': {'odd_even_rows': True},
                                  'dark_current': {'skip': False},
                                #   'dark_current': {'skip': True},
                                  'ipc': {'skip': True},
                                  'jump': {'rejection_threshold': jump_threshold, 'maximum_cores':self.ramp_cores}, # ERS recommend threshold = 4 NIRCam, = 8 MIRI
                                  'ramp_fit': {'save_calibrated_ramp': False, 'maximum_cores':self.ramp_cores}},
                           subdir=subdir,
                           verbose = True)
    
    def stage2(self, subdir='stage2'):
        coron2pipeline.run_obs(database=self.Database,
                                steps={'outlier_detection': {'skip': False}},
                                subdir=subdir)
    
    def nircam_jwst_sub(self, subdir='jwst'):
        """ TODO: rename this for jwst pipeline subtraction. Maybe separate step from stage2? """
        coron2pipeline.run_obs(database=self.Database,
                           steps={'outlier_detection': {'skip': False}},
                           subdir=subdir)
        coron3pipeline.run_obs(database=self.Database, #performs jwst pipeline subtraction
                            steps={'klip': {'truncate': 100}},
                            subdir='stage3')

    def miri_image_processing(self):
        ImageTools = imagetools.ImageTools(database=self.Database)

        if self.godoy == False:
            # Remove the first frame due to reset switch charge delay. Only required for MIRI.
            ImageTools.remove_frames(index=[0],
                                    #  frame=None, # old version had this param
                                    types=['SCI', 'SCI_BG', 'REF', 'REF_BG'],
                                    subdir='removed')

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
        ImageTools.crop_frames(npix=[13, 60, 7, 7], # [left, right, bottom, top]
                            types=['SCI', 'SCI_BG', 'REF', 'REF_BG'],
                            subdir='cropped')

        # Replace all nans. Same has not changed
        ImageTools.replace_nans(cval=0.,
                                types=['SCI', 'SCI_BG', 'REF', 'REF_BG'],
                                subdir='nanreplaced')
        
        if self.godoy:
            # New bkg subtraction method
            ImageTools.subtract_background_godoy(subdir='bgsub_godoy')
        else:
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
        
    def nircam_image_processing(self):
        ImageTools = imagetools.ImageTools(database=self.Database)

        # Subtract the median from each frame. 
        # Clip everything brighter than 5- sigma from the background before computing the median.
        ImageTools.subtract_median(types=['SCI', 'SCI_TA', 'SCI_BG',
                                        'REF', 'REF_TA', 'REF_BG'],
                                subdir='medsub')
        
        # bad pixel cleaning
        # bpmap = np.zeros((320, 320))
        # bpmap[153, 153] = 1
        # bpmap[159, 142] = 1
        # bpmap[159, 143] = 1
        # bpmap[159, 157] = 1
        # bpmap[159, 158] = 1
        # bpmap[159, 164] = 1
        # bpmap[164, 141] = 1
        # bpmap[165, 146] = 1
        # bpmap[169, 140] = 1
        # bpmap[177, 120] = 1
        # custom_kwargs = {}
        # custom_kwargs['JWST_NIRCAM_NRCA2_F200W_MASKRND_MASKA335R_SUB320A335R'] = bpmap
        # custom_kwargs['JWST_NIRCAM_NRCALONG_F250M_MASKRND_MASKA335R_SUB320A335R'] = bpmap
        # custom_kwargs['JWST_NIRCAM_NRCALONG_F300M_MASKRND_MASKA335R_SUB320A335R'] = bpmap
        # custom_kwargs['JWST_NIRCAM_NRCALONG_F356W_MASKRND_MASKA335R_SUB320A335R'] = bpmap
        # custom_kwargs['JWST_NIRCAM_NRCALONG_F410M_MASKRND_MASKA335R_SUB320A335R'] = bpmap
        # custom_kwargs['JWST_NIRCAM_NRCALONG_F444W_MASKRND_MASKA335R_SUB320A335R'] = bpmap
        # ImageTools.fix_bad_pixels(method='bpclean+custom+timemed+dqmed+medfilt', # default:timemed+dqmed+medfilt

        ImageTools.fix_bad_pixels(method='bpclean+timemed+dqmed+medfilt',
                              bpclean_kwargs={'sigclip': 3,
                                              'shift_x': [-1, 0, 1],
                                              'shift_y': [-1, 0, 1]},
                              custom_kwargs={},
                              timemed_kwargs={},
                              dqmed_kwargs={'shift_x': [-1, 0, 1],
                                            'shift_y': [-1, 0, 1]},
                              medfilt_kwargs={'size': 4},
                              subdir='bpcleaned')

        # ImageTools.fix_bad_pixels(method='bpclean+timemed+dqmed+medfilt', # default:timemed+dqmed+medfilt
        #                         bpclean_kwargs={'sigclip': 5,
        #                                         'shift_x': [-1, 0, 1],
        #                                         'shift_y': [-1, 0, 1]},
        #                         custom_kwargs=custom_kwargs,
        #                         timemed_kwargs={},
        #                         dqmed_kwargs={'shift_x': [-1, 0, 1],
        #                                         'shift_y': [-1, 0, 1]},
        #                         medfilt_kwargs={'size': 4},
        #                         subdir='bpcleaned')
        
        # Replace all nans in the data with a constant value.
        ImageTools.replace_nans(cval=0.,
                                types=['SCI', 'SCI_BG', 'REF', 'REF_BG'],
                                subdir='nanreplaced')
        
        # Recenter frames. Before, update NIRCam coronagraphic mask centers to
        # on-sky values measured by Jarron. Might not be required for simulated
        # data!
        ImageTools.update_nircam_centers()

        # Recenter frames so that the host star position is data.shape // 2. 
        # For NIRCam coronagraphy, use a WebbPSF model to determine the star position 
        # behind the coronagraphic mask for the first SCI frame. 
        # Then, shift all other SCI and REF frames by the same amount.
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
        
    def klip_subtraction(self, subdir='klipsub'):
        pyklippipeline.run_obs(database=self.Database,
                        kwargs={'mode': self.subtraction,
                                 'annuli': self.annuli,
                                 'subsections': self.subsections,
                                 'numbasis': self.KLmodes,
                                 'algo': 'klip',
                                 'save_rolls': True},
                        subdir=subdir)

    def raw_contrast(self, subdir='rawcon'):
        # Initialize spaceKLIP analysis tools class.
        AnalysisTools = analysistools.AnalysisTools(self.Database)

        # Compute raw contrast.
        if self.companions == False:
            AnalysisTools.raw_contrast(starfile=self.photo_path,
                                        spectral_type=self.SpT,
                                    #    companions=[[0.431, -0.717, 2.]],
                                        subdir='rawcon')
        else:
            AnalysisTools.raw_contrast(starfile=self.photo_path,
                                        spectral_type=self.SpT,
                                        companions=self.companions,
                                        subdir='rawcon')
            
    def calibrated_contrast(self, subdir='calcon'):
        # Initialize spaceKLIP analysis tools class.
        AnalysisTools = analysistools.AnalysisTools(self.Database)

        # Compute raw contrast.
        if self.companions == False:
            AnalysisTools.calibrate_contrast(rawcon_subdir='rawcon',
                                 #  companions=[[0.431, -0.717, 2.]],
                                    use_saved=False )
        else:
            AnalysisTools.calibrate_contrast(rawcon_subdir='rawcon',
                                    companions=self.companions,
                                    use_saved=False )

    def find_contrast_paths(self):
        contrast_paths = []
        for sub_method in self.subtraction:
            for ann in self.annuli:
                for subsec in self.subsections:
                    contrast_paths.append(glob.glob(self.odir + 'klipsub/{}_NANNU{}_NSUBS{}*all.fits'.format(sub_method, ann, subsec))[0])

        return contrast_paths
        
    def run_pipeline(self, godoy=False, sub_path='padded', jump_threshold=8.):
        self.godoy = godoy

        while self.stage != self.final_stage:
            if self.stage == 'stage1':
                print('\n', 'Performing Stage1 for {}'.format(self.star), '\n')
                self.load_database()
                if godoy:
                    self.stage1_mask(subdir=self.stage1_subdir)
                else:
                    self.stage1_subdir = 'stage1_th{}'.format(jump_threshold)
                    self.stage1(subdir=self.stage1_subdir)
                self.stage = 'stage2'

            elif self.stage == 'stage2':
                print('\n', 'Performing Stage2 for {}'.format(self.star), '\n')
                self.load_database(stage1_path=self.stage1_subdir)
                self.stage2()

                # if self.instr == 'nircam':
                #     self.stage2_nircam()
                # else:
                #     self.stage2()
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
                self.load_database(subpath=sub_path)
                self.klip_subtraction()
                self.stage = 'rawcon'
            
            elif self.stage == 'rawcon':
                print('\n', 'Performing Raw contrast for {}'.format(self.star), '\n')
                self.contrast_path = self.find_contrast_paths()
                self.load_database()
                self.raw_contrast()
                self.stage = 'calcon'
            
            elif self.stage == 'calcon':
                print('\n', 'Performing calibrated contrast for {}'.format(self.star), '\n')
                if self.contrast_path == None:
                    self.contrast_path = self.find_contrast_paths()
                self.load_database()
                self.calibrated_contrast()
                self.stage = 'end'
                
            
