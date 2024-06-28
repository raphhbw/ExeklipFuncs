# ExeklipFuncs

Functions used to run [SpaceKLIP](https://github.com/kammerje/spaceKLIP) for MIRI and NIRCam JWST data.

These functions have been built using the [ExeterSpaceKLIP](https://github.com/raphhbw/ExeterSpaceKLIP) package. For NIRCam reductions use the `develop` branch and for MIRI use the `miri_improvements` branch of SpaceKLIP.

### Run code
To run the code make sure you are using the correct branch of SpaceKLIP for the desired type of subtraction. Run `python multi_sklip.py --help` to check the different options available when running the code.
```
python multi_sklip.py --p params/params.json
```

### Files
- `params/params.json`: example list of parameters that can be given to `multi_sklip.py`. Parameters include:
                - `star`: name of star
                - `instr`: instrumnet (nircam or miri)
                - `bkg`: None for nircam, e.g. ["8005","8006"] for miri
                - `subtraction`: subtraction mode e.g. ["ADI", "RDI", "ADI+RDI"]
                - `KLmodes`: number of KLmodes for PSF subtraction
                - `subsect` and `annuli`: number of subsections and annuli for the PSF subtraction (e.g. [1,2,3])
                - `uncal_dir`: path to the directory with the uncal files
                - `odir`: path to the directory where the SpaceKLIP outputs will generate files
                - `SpT`: spectral type of the star
                - `photo_path`: path to the sed file for the star. Used for contrast calculations
                - `companions`: companion fitting parameters in SpaceKLIP (default: 0, i.e. False)
                - `start_stage` and `final_stage`: what stage of the reduction do you want to start from and what stage do you want to **finish before**. Options: 'stage1', 'stage2', 'img_proc', 'sub', 'rawcon', 'calcon', 'end'
                - `godoy`: godoy method for background subtraction in miri. Flag turned on or off
                - `skip_1_f`: flag to skip 1/f correction in nircam. Can be useful in some cases with extended structures for example
- `multi_sklip.py`: Main code to run. Will call the correct function from `sklip.py` where all the different SpaceKLIP stages are coded. Variable parameters can be seen using the help:
                - `--p param_file.json`: param file used for the reduction
                - `--nice niceness`: useful when working on server
                - `thr jump_threshold`: jump threshold in the ramp fitting stage1. Used for nircam only
                - `--cores cores`: set the number of cores to use for the ramp fitting stage
- `sklip.py`: All the SpaceKLIP functions are called in this file.