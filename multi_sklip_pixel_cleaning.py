""" Code updated
9-01-24
"""

import argparse, os, json, sys
from sklip_pixel_cleaning import Sklip
import numpy as np

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Help')
    parser.add_argument('--p', metavar='param.json', help='parameter file name (default:param.json)')
    parser.add_argument('--nice', metavar='niceness', default=10, help='niceness of job (default: 10)')
    parser.add_argument('--cores', metavar='cores', default='all', 
                        help='cores used in the ramp fitting (default: all) - options: quarter, half, all')
    args = parser.parse_args()

    os.nice(int(args.nice))

    with open(args.p) as paramfile:
            param = json.load(paramfile)

    stages = ['stage1', 'stage2', 'img_proc', 'sub', 'rawcon', 'calcon', 'end']
    # jump_thresholds = range(2,9)

    bp_file = '/Users/rbendahan/Astronomy/roach/MIRI-1140/HD92945/test_custom_pb_clean/klipsub/ADI+RDI_NANNU1_NSUBS1_JWST_MIRI_MIRIMAGE_F1140C_NONE_4QPM_1140_MASK1140-KLmodes-all_roll2_init.fits'

    for target in param.keys():

            print('##########################################')
            print(target, '----------')
            print(param[target])
            star_param = param[target]
            print('##########################################')

            # if jump_th == 1:
            #     stage = 'rawcon'
            # else:
            stage = star_param["start_stage"]
            final_stage = star_param["final_stage"]

            updated_odir = star_param["odir"]
            # print(updated_odir)

            Pipeline = Sklip(star=star_param["star"], 
                        bkg=star_param["bkg"], 
                        instr=star_param["instr"],
                        subtraction=star_param["subtraction"], 
                        KLmodes=star_param["KLmodes"],
                        subsect=star_param["subsect"], 
                        annuli=star_param["annuli"], 
                        uncal_dir=star_param["uncal_dir"], 
                        odir=updated_odir,
                        SpT=star_param["SpT"],
                        photo_path=star_param["photo_path"],
                        companions=star_param["companions"],
                        stage=stage,
                        final_stage=final_stage
                        )
            print('--------------')
            # if jump_th == 1:
            #     Pipeline.stage1_subdir = 'stage1_masked_thr1'
            # Pipeline.load_database()
            
            Pipeline.run_pipeline(godoy=star_param["godoy"], sub_path='bgsub_godoy', bp_file=False)

            # Pipeline.find_contrast_paths()
            # print(Pipeline.contrast_path)

            print('\n' + '##########################################')
            print('Finished with {}'.format(target), '----------')
            print('##########################################')

