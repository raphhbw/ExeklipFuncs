if __name__ == '__main__':
    import argparse, os, json
    from sklip import Sklip

    # setting argparse options
    parser = argparse.ArgumentParser(description='Parameter option for multi_sklip.py')
    parser.add_argument('--p', metavar='param.json', help='parameter file name (default:param.json)')
    parser.add_argument('--nice', metavar='niceness', default=10, help='niceness of job (default: 10)')
    parser.add_argument('--thr', metavar='jump_threshold', default=4., help='jump threshold in ramp fitting (default: 4.)')
    parser.add_argument('--cores', metavar='cores', default='all', 
                    help='cores used in the ramp fitting (default: all) - options: quarter, half, all')
    args = parser.parse_args()

    # Set niceness
    os.nice(int(args.nice))

    # Load params
    with open(args.p) as paramfile:
        param = json.load(paramfile)

    # Stage option reminder
    stages = ['stage1', 'stage2', 'img_proc', 'sub', 'rawcon', 'calcon', 'end']

    # Looping over the different stars in the param file
    for target in param.keys():
        print('##########################################')
        print(target, '----------')
        print(param[target])
        star_param = param[target]
        print('##########################################')

        # Set up the Sklip pipeline with info from param file
        Pipeline = Sklip(star=star_param["star"], 
                    bkg=star_param["bkg"], 
                    instr=star_param["instr"],
                    subtraction=star_param["subtraction"], 
                    KLmodes=star_param["KLmodes"],
                    subsect=star_param["subsect"], 
                    annuli=star_param["annuli"], 
                    uncal_dir=star_param["uncal_dir"], 
                    odir=star_param["odir"],
                    SpT=star_param["SpT"],
                    photo_path=star_param["photo_path"],
                    companions=star_param["companions"],
                    stage=star_param["start_stage"],
                    final_stage=star_param["final_stage"],
                    skip_1_f=star_param['skip_1_f']
                    )
        print('--------------')
        
        # Run pipeline from sklip.py
        # sub_path defines the files used to perform the PSF subtraction on (default: padded)
        Pipeline.run_pipeline(godoy=star_param["godoy"], sub_path='padded', jump_threshold=args.thr)

        print('\n' + '##########################################')
        print('Finished with {}'.format(target))
        print('##########################################')

