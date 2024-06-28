""" Code updated
9-01-24
"""

import argparse, os, json, sys
from sklip_test import Sklip

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
    
#     stage = 'stage1'
#     final_stage = 'rawcon'
#     stage = param["start_stage"]
#     final_stage = param["final_stage"]

    print('##########################################')
    print(target)
    print(param[target])
    star_param = param[target]
    print('##########################################')

    stage = star_param["start_stage"]
    final_stage = star_param["final_stage"]

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
    print('--------------')
#     Pipeline.load_database()
#     data = Pipeline.Database
# #     data.summarize()
# #     df = Pipeline.Database.to_pandas()
#     for i, key in enumerate(data.obs.keys()):
#         # print(data.obs[key].colnames)
#         science = data.obs[key][data.obs[key]['TYPE'] == 'SCI']['FITSFILE'][0]
#         ref = data.obs[key][data.obs[key]['TYPE'] == 'REF']['FITSFILE'][0]
#     # sys.exit()
#     Pipeline.optimise_groups(RefPath=ref, SciPath= science)
    Pipeline.run_pipeline(godoy=star_param["godoy"])
    print('\n' + '##########################################')
    print('Finished with {}'.format(target))
    print('##########################################')

