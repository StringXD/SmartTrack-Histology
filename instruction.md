# Histology tools

Some code to work with the Allen Inst Mouse Brain CCF data, specifically the 10µm voxel 2017 version.
In this repository, we added some new features in SHARP-Track to make it easier to use.


## Requirements
- A computer mouse with a scroll wheel
- MATLAB (R2017 on a Windows computer was used for all testing)
- Fiji (optional) (http://fiji.sc/)
- This repository (add all folders and subfolders to your MATLAB path)
- [The npy-matlab repository](http://github.com/kwikteam/npy-matlab)
- [The Allen Mouse Brain Atlas volume and annotations](http://data.cortexlab.net/allenCCF/) (download all 4 files from this link) 
If you have access, you could also locate the data files at //zserver/Lab/Atlas/allenCCF. Alternatively, see setup_utils to download and preprocess the files yourself. See also https://alleninstitute.github.io/AllenSDK/reference_space.html to access the data directly via the Allen Inst python API.


## Format convertion

To use SHARP-Track, images have to be converted into .tif format. For myself, I use olympus imaging system the most and their data format is .vsi, so I made a macro to convert image format and also do some preprocessing like divide or merge channels, adjust brightness and contrast, and compute micron-pixel ratio.(See NPpreprocess.txt)

## Multiboxing

