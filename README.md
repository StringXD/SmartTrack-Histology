# SmartTrack Histology

## What does it do ?

The basis of this repository is [allenCCF](https://github.com/cortex-lab/allenCCF), we thank very much to its developers for providing such a fantastic tool for validate electrode probe track in histology image. To make it easier to use, we added some new functions like image format conversion and designed track preview. These new functions are listed below:

## Image format conversion

*For now, only conversion from olympus .vsi file to .tif format is supported.*

Before using allenCCF, images need to be converted to .tif format and the micron pixel ratio needs to be specified. But for our own daily use, what we get after brain slice imaging with Olympus imaging system are .vsi files. So we are using [FiJi](http://fiji.sc/) to  convert image format, roughly adjust brightness and contrast, and compute micron pixel ratio for batches of images with our script.

## Design probe tracks

We care about how many brain areas does the designed track go through and what they are. The [Display_Designed_Track.m](https://github.com/StringXD/SmartTrack-Histology/blob/master/SHARP-Track/Display_Designed_Track.m) and [Display_Designed_Track3.m](https://github.com/StringXD/SmartTrack-Histology/blob/master/SHARP-Track/Display_Designed_Track3.m) provide a WYSIWYG way to preview your track given the implantion angle, depth and coordinates.

## Using images with different contrast for registration and picking ROI

This feature is used when there is difficulty both differentiating different brain areas when registration and picking ROIs after registration. In this condition, user can use two set of images, one for registration, the other for picking ROIs.

** * Note this two set of images must be identical except color contrast**

## And some optomization

