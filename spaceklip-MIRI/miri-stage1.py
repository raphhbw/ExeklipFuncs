import os, time
from spaceKLIP import database, coron1pipeline

if __name__ == '__main__':

    # Directories
    os.nice(10) # change niceness

    t0 = time.time()
    star = 'HD206893'

    idir = '/data/rb941/JWST/MIRI-1140/{}/uncal/'.format(star)
    odir = '/data/rb941/JWST/MIRI-1140/{}/spaceklip/'.format(star)

    # Initialise Database
    Database = database.Database(output_dir=odir)

    fitsfiles = sorted([idir + f for f in os.listdir(idir) if f.endswith('_uncal.fits')])
    # bgfitsfiles =  sorted([idir + f for f in os.listdir(idir) if '8011' in f or '8011' in f]) # HD21997
    # bgfitsfiles =  sorted([idir + f for f in os.listdir(idir) if '8006' in f or '8010' in f]) # HD92945
    # bgfitsfiles =  sorted([idir + f for f in os.listdir(idir) if '8001' in f or '8005' in f]) # GJ14

    Database.read_jwst_s012_data(datapaths=fitsfiles, 
                                #  bgpaths=bgfitsfiles
                                 ) 

    Database.summarize()

    # coron1pipeline.run_obs(database=Database,
    #                        steps={'saturation': {'n_pix_grow_sat': 1,
    #                                              'grow_diagonal': False},
    #                               'refpix': {'odd_even_columns': True,
    #                                          'odd_even_rows': True,
    #                                          'nlower': 0,
    #                                          'nupper': 0,
    #                                          'nleft': 0,
    #                                          'nright': 0,
    #                                          'nrow_off': 0,
    #                                          'ncol_off': 0},
    #                               'dark_current': {'skip': True},
    #                               'jump': {'rejection_threshold': 8.,
    #                                        'three_group_rejection_threshold': 8.,
    #                                        'four_group_rejection_threshold': 8.,
    #                                        'maximum_cores':'all'},
    #                               'ramp_fit': {'save_calibrated_ramp': False,
    #                                             'maximum_cores':'all'}},
    #                        subdir='stage1')

    
    t1 = time.time()
    print("Total execution time: {}".format(t1-t0))