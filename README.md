# FORS2 - Photopolarimetry

## General Information
The calibration scripts (here made available) were designed to process raw ESO VLT-FORS2/IPOL fits files, as such they expect a specific header structure.
These scripts, as well as the reduction scripts, purpose is to produce polarimetric descriptions for extended sources.
These scripts were tested in a system with the following specs 
* Processor: Intel® Core™ i7-8750H CPU @ 2.20GHz × 12,
* RAN: 31.2 GB,
* OS: Ubuntu 20.04.5 LTS 64-bit
* GNOME: 3.36.8
* R: 4.2.1
* Python: 3.8

Nevertheless, all scripts except for "process_BIAS.R" require less than 8GB of RAM at any given time. "process_BIAS.R" was designed to adapt itself to the available RAM.

## Setting up
First, you'll need to edit the "home_folder" variable, in each script. This variable should point to folder that holds all (or at least most of) the files/folders that will be called be the scripts.

Second. the calibration and reduction scripts will call functions within the libraries contained in the "Libraries/" folder. You'll need to edit the variable "libs_folder", in each script, so that it points to the folder holding those libraries. 

Third, within the library "process_lib.R" you'll also need to properly point to the folder that will hold the "sexyDarce.py" library, and to the your preferred python version (all functions were tested within python 3.8, some were also tested within 3.10).

All other paths will be setup via prompt. I recommend however accepting the folder structure I planned, even if it is not be optimal, otherwise it may be harder to troubleshoot.

## Running the calibration scripts

After setting up, if you have calibration data you should run the scripts in the following order:

1- "process_BIAS.R", to generate master bias files.

2- "process_FLAT.R", to generate individual CCD Chips master flat files (requires outputs of "process_BIAS.R").

- "process_ACQ.R", to generate joined CCD Chips, calibrated Acquisition files.

- "process_STD.R", to generate joined CCD Chips, calibrated standard files.

- "process_OBJ.R", to generate joined CCD Chips, calibrated science files.

These last three scripts can be run in any order. 
The outputs of "process_ACQ.R" will be used to generate masks and astrometry files in the reduction scripts.
The "process_OBJ.R" will pre-process science images into ordinary and extraordinary beam maps that will serve as input for the reduction scripts.

## Running the reduction scripts

The calibration data you should run the scripts in the following order:

1- "process_BEAM.R", to produce individual offset polarimetric maps.

2- "merge_OFFSET.R", to produce joined offsets polarimetric maps.

3- "crop_MERGED.R", to crop out the areas of the image that are not of interest.

## List of dependencies

### System
Astrometry.net (http://astrometry.net/use.html)

### R libraries
##### (some of these may no longer be required but were not pruned from the code in due time)
colorspace, dipsaus, EBImage, fields, FITSio, ggplot2, graphics, jjb, MASS, plot3D, plotrix, profmem, reticulate, spam, stringr 

### Python libraries
astropy,  astroquery, lacosmic, math, matplotlib, numpy, photutils, scipy, sep, random, warnings
