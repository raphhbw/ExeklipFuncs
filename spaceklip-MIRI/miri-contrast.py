import os, argparse, sys

from spaceKLIP import database, analysistools

if __name__ == '__main__':
   parser = argparse.ArgumentParser(description='Help')
   parser.add_argument('--nice', metavar='niceness', default=10, help='niceness of job (default: 10)')
   parser.add_argument('--star', metavar='star', help='input name of star')
   args = parser.parse_args()

   os.nice(int(args.nice))
   star = args.star
   # print(star, type(star))
   # sys.exit()

   if star == 'HD92945':
      SpT = 'K1V'
      print(star, SpT)
   elif star == 'HD107146':
      SpT = 'G2V'
      print(star, SpT)
   elif star == 'HD206893':
      SpT = 'F5V'
      print(star, SpT)

   #  odir = '/data/rb941/JWST/MIRI-1140/{}/spaceklip/'.format(star)
   odir = '/data/rb941/JWST/MIRI-1140/{}/opt_groups/'.format(star)

   Database = database.Database(output_dir=odir)
   datapaths = [odir + 'klipsub/ADI_NANNU1_NSUBS1_JWST_MIRI_MIRIMAGE_F1140C_NONE_4QPM_1140_MASK1140-KLmodes-all.fits',
               odir + 'klipsub/RDI_NANNU1_NSUBS1_JWST_MIRI_MIRIMAGE_F1140C_NONE_4QPM_1140_MASK1140-KLmodes-all.fits',
               odir + 'klipsub/ADI+RDI_NANNU1_NSUBS1_JWST_MIRI_MIRIMAGE_F1140C_NONE_4QPM_1140_MASK1140-KLmodes-all.fits']
   Database.read_jwst_s3_data(datapaths)

   # Initialize spaceKLIP analysis tools class.
   AnalysisTools = analysistools.AnalysisTools(Database)

   # Compute raw contrast.
   AnalysisTools.raw_contrast(starfile='/data/rb941/photometry/{}_phoenix_sol+modbb_disk_r_.txt'.format(star.lower()),
                              spectral_type=SpT,
                           #    companions=[[0.431, -0.717, 2.]],
                              subdir='rawcon')

   AnalysisTools.calibrate_contrast(rawcon_subdir='rawcon',
                                 #  companions=[[0.431, -0.717, 2.]],
                                    use_saved=False )