# =============================================================================
# IMPORTS
# =============================================================================

import os

from spaceKLIP import database, pyklippipeline


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    
   # Set the input and output directories and grab the input FITS files.
   os.nice(10) # change niceness

   idir = '/data/rb941/JWST/MIRI-1140/HD92945/spaceklip/stage1/'
   odir = '/data/rb941/JWST/MIRI-1140/HD92945/spaceklip/'

   idir = idir.replace('stage1', 'bgsub')
   Database = database.Database(output_dir=odir)
   fitsfiles = sorted([idir + f for f in os.listdir(idir) if f.endswith('_calints.fits')])
#  bgfitsfiles =  sorted([f for f in fitsfiles if '8011' in f or '8015' in f]) # HD206893
#  bgfitsfiles =  sorted([f for f in fitsfiles if '8001' in f or '8005' in f]) # HD107146
   bgfitsfiles =  sorted([f for f in fitsfiles if '8006' in f or '8010' in f]) # HD92945
   # bgfitsfiles =  sorted([f for f in fitsfiles if '8011' in f or '8011' in f]) # HD21997
   Database.read_jwst_s012_data(datapaths=fitsfiles,
                              bgpaths=bgfitsfiles)
   
   pyklippipeline.run_obs(database=Database,
                     kwargs={'mode': ['ADI', 'RDI', 'ADI+RDI'],
                              'annuli': [1, 6],
                           # annuli & subsections defines the different sections where 
                           # klip does subtraction
                              'subsections': [1, 6],
                              'numbasis': [10],
                              # numbasis = klmode
                                 'algo': 'klip',
                                 'save_rolls': True},
                        subdir='klipsub')