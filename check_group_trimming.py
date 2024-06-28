""" Code updated
9-01-24
"""
if __name__ == '__main__':
    import argparse, os, json, sys
    from sklip_v1 import Sklip
    import numpy as np
    # from sklip_v0 import Sklip


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

    for target in param.keys():

        print('##########################################')
        print(target, '----------')
        print(param[target])
        star_param = param[target]
        print('##########################################')

        stage = 'stage1'
        final_stage = star_param["final_stage"]

        updated_odir = star_param["odir"]

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
        
        Pipeline.run_pipeline(godoy=star_param["godoy"])

        # Pipeline.find_contrast_paths()
        # print(Pipeline.contrast_path)

        print('\n' + '##########################################')
        print('Finished with {}'.format(target))
        print('##########################################')

