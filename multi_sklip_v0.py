""" Code from 2023 """

import argparse, os, json
from sklip import Sklip

parser = argparse.ArgumentParser(description='Help')
parser.add_argument('--p', metavar='param.json', help='parameter file name (default:param.json)')
parser.add_argument('--nice', metavar='niceness', default=10, help='niceness of job (default: 10)')
parser.add_argument('--cores', metavar='cores', default='all', 
                    help='cores used in the ramp fitting (default: all) - options: quarter, half, all')
args = parser.parse_args()

os.nice(int(args.nice))

with open(args.p) as paramfile:
        param = json.load(paramfile)

stages = ['stage1', 'stage2', 'img_proc', 'sub', 'rawcon', 'calcon']

for target in param.keys():

#     if target != 'HD21997':
#           continue
    
    stage = 'stage1'
    final_stage = 'rawcon'

    print('##########################################')
    print(target)
    print(param[target])
    star_param = param[target]
    print('##########################################')

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
                  stage=stage,
                  final_stage=final_stage
                  )

    Pipeline.run_pipeline()
    print('\n' + '##########################################')
    print('Finished with {}'.format(target))
    print('##########################################')

