rm(list=ls())
require(FITSio)
require(stringr)
require(jjb)
require(colorspace)
require(graphics)
require(ggplot2)
require(plot3D)
require(fields)
require(plotrix)
require(reticulate)
require(dipsaus)
require(MASS)
require(spam)
library(EBImage)
options(scipen = 999)
options(digits = 8)

home_folder <- "/media/joaomfras/GalaxyPol/Pol-Gal/"
libs_folder <- paste0(home_folder, 
                      "Polarimetric-Imaging-Reduction-Scripts-(FORS2)/Commit/")
lib_head_path <- paste0(libs_folder, "fetch_fits_header_info.R")
lib_pro_path <- paste0(libs_folder, "process_lib.R")
source(lib_head_path)
source(lib_pro_path)
rm(lib_head_path, lib_pro_path)

################################################################################
################################# User Prompts #################################
################################################################################
og_obj <- "NGC-3351"                                                           #
print("Provide the name of the target: ")                                      #
og_obj <- edit(og_obj)                                                         #
                                                                               #
psf_est <- 5                                                                   #
print(paste0("Provide an estimation for the PSF radius, this will be used to ",#
             "estimate the background flux: "))                                #
psf_est <- edit(psf_est)                                                       #
while(psf_est < 1){                                                            #
  print(paste0("This value must be >= 1: "))                                   #
  psf_est <- edit(psf_est)                                                     #
}                                                                              #
                                                                               #
q_prob <- .85                                                                  #
print(paste0("After the beams are background subtracted, a gaussian will be ", #
             "fitted to the negative values and mirrored to the positive side",#
             ". When plotting a histogram of these values and the gaussian ",  #
             "distribution on top of it we can see that at some point, there ",#
             "are more positive values than the distribution estimates. ",     #
             "Those are likely to be signal of intereset. Define a minimum ",  #
             "ratio between the distribution frequency estimation and the ",   #
             "actual frequency of values in a given positive interval, so ",   #
             "that that interval is still removed from the signal: "))         #
q_prob <- edit(q_prob)                                                         #
while(q_prob <= 0 | q_prob >= 1){                                              #
  print(paste0("This value must be larger than 0 and smaller than 1: "))       #
  q_prob <- edit(q_prob)                                                       #
}                                                                              #
                                                                               #
binS <- 5                                                                      #
print("Define a valid bin size (integer > 0) for the maps to be produced:")    #
binS <- edit(binS)                                                             #
                                                                               #
while(binS < 1 || binS %% 1 != 0){                                             #
  print("Invalid input. Please input a positive integer.")                     #
  binS <- edit(binS)                                                           #
}                                                                              #
                                                                               #
og_folder <- paste0(home_folder, og_obj, "/Processed_OBJ/")                    #
                                                                               #
MWstars_folder <- paste0(home_folder, "Milky_Way_stars/")                      #
print("Provide an output path for Milky Way star polarimetry: ")               #
MWstars_folder <- edit(MWstars_folder)                                         #
if(!dir.exists(MWstars_folder)){                                               #
  mkdir(MWstars_folder)                                                        #
}                                                                              #
                                                                               #
acq_folder <- paste0(home_folder, og_obj, "/Processed_ACQ/Merged_ACQ/")        #
                                                                               #
beam_folder <- paste0(og_folder, "Beams/")                                     #
input_folder <- paste0(beam_folder, "Obs/")                                    #
stopifnot(dir.exists(input_folder))                                            #
                                                                               #
bin_fold <- paste0("bin", binS,"/")                                            #
                                                                               #
instV <- 1                                                                     #
print(paste0("Chose which version of FORS2 polarization files to use: 0 (GGM2",#
             "019) or 1 (JRS2023)"))                                           #
instV <- edit(instV)                                                           #
                                                                               #
while(instV != 0 && instV != 1){                                               #
  print("Invalid input. Please input either 0 or 1.")                          #
  instV <- edit(instV)                                                         #
}                                                                              #
                                                                               #
inst_folder <- paste0(home_folder, "FORS2_QU_CORR/inla/")                      #
print("Provide the folder path for FORS2 polarization files:")                 #
inst_folder <- edit(inst_folder)                                               #
                                                                               #
if(!instV){                                                                    #
  use_inst_bin <- 0                                                            #
  save_JRS2023 <- 1                                                            #
  print(paste0("This script will adjust the GGM2019 files dimensions in order",#
               " to match those of the beams files. This will result in the J",#
               "RS2023 version of the FORS2 polarization files. Do you wish t",#
               "o save these for later use? 0 (No) or 1 (Yes)"))               #
  save_JRS2023 <- edit(save_JRS2023)                                           #
                                                                               #
  while(save_JRS2023 != 0 && save_JRS2023 != 1){                               #
    print("Invalid input. Please input either 0 or 1.")                        #
    save_JRS2023 <- edit(save_JRS2023)                                         #
  }                                                                            #
                                                                               #
  if(binS > 1){                                                                #
    save_inst_bin <- 1                                                         #
    print(paste0("This script will bin the JRS2023 files pixels in order to p",#
                 "erform some operations along with other files later on. The",#
                 "se files will be needed if you run 'merge_OFFSETS_v3.R'. Do",#
                 " you wish to save the resulting files for later use? 0 (NO)",#
                 " or 1 (Yes)"))                                               #
    save_inst_bin <- edit(save_inst_bin)                                       #
                                                                               #
    while(save_inst_bin != 0 && save_inst_bin != 1){                           #
      print("Invalid input. Please input either 0 or 1.")                      #
      save_inst_bin <- edit(save_inst_bin)                                     #
    }                                                                          #
  }else{                                                                       #
    save_inst_bin <- 0                                                         #
  }                                                                            #
                                                                               #
  if(save_JRS2023 || save_inst_bin){                                           #
    inst_new_folder <- paste0(inst_folder, "vJRS2023/")                        #
    print("Provide the folder path for FORS2 JRS2023 polarization files:")     #
    inst_new_folder <- edit(inst_new_folder)                                   #
                                                                               #
    if(!dir.exists(inst_new_folder)){                                          #
      mkdir(inst_new_folder)                                                   #
    }                                                                          #
                                                                               #
    if(save_inst_bin){                                                         #
      inst_bin_folder <- paste0(inst_new_folder, "Binned/", bin_fold)          #
      if(!dir.exists(inst_bin_folder)){                                        #
        mkdir(inst_bin_folder)                                                 #
      }                                                                        #
    }                                                                          #
  }                                                                            #
}else{                                                                         #
  if(binS > 1){                                                                #
    use_inst_bin <- 1                                                          #
    print(paste0("Do you wish to input a previously binned version of FORS2 J",#
                 "RS2023 files? 0 (No) or 1 (Yes)"))                           #
    use_inst_bin <- edit(use_inst_bin)                                         #
                                                                               #
    while(use_inst_bin != 0 && use_inst_bin != 1){                             #
      print("Invalid input. Please input either 0 or 1.")                      #
      use_inst_bin <- edit(use_inst_bin)                                       #
    }                                                                          #
  }else{                                                                       #
    use_inst_bin <- 0                                                          #
  }                                                                            #
                                                                               #
  inst_folder <- paste0(inst_folder, "vJRS2023/")                              #
  print("Provide the folder path for FORS2 JRS2023 files you wish to use:")    #
  inst_folder <- edit(inst_folder)                                             #
                                                                               #
  inst_bin_folder <- paste0(inst_folder, "Binned/", bin_fold)                  #
  print(paste0("Provide the folder path for binned FORS2 JRS2023 files you wi",#
               "sh to use:"))                                                  #
  inst_folder <- edit(inst_folder)                                             #
                                                                               #
  if(!use_inst_bin){                                                           #
    if(binS > 1){                                                              #
      save_inst_bin <- 1                                                       #
      print(paste0("This script will bin the JRS2023 files pixels in order to",#
                   " perform some operations along with other files later on.",#
                   " Do you wish to save the resulting files for later use? 0",#
                   " (No) or 1 (Yes)"))                                        #
      save_inst_bin <- edit(save_inst_bin)                                     #
                                                                               #
      while(save_inst_bin != 0 && save_inst_bin != 1){                         #
        print("Invalid input. Please input either 0 or 1.")                    #
        save_inst_bin <- edit(save_inst_bin)                                   #
      }                                                                        #
                                                                               #
      if(save_inst_bin){                                                       #
        if(!dir.exists(inst_bin_folder)){                                      #
          mkdir(inst_bin_folder)                                               #
        }                                                                      #
      }                                                                        #
    }                                                                          #
  }                                                                            #
}                                                                              #
                                                                               #
obs_bin_folder <- paste0(input_folder, "Binned/", bin_fold)                    #
if(!dir.exists(obs_bin_folder)){                                               #
  mkdir(obs_bin_folder)                                                        #
}                                                                              #
obs_stat_folder <- paste0(obs_bin_folder, "Stats/")                            #
if(!dir.exists(obs_stat_folder)){                                              #
  mkdir(obs_stat_folder)                                                       #
}                                                                              #
                                                                               #
noSky_beam_folder <- paste0(beam_folder, "No_Sky/")                            #
if(!dir.exists(noSky_beam_folder)){                                            #
  mkdir(noSky_beam_folder)                                                     #
}                                                                              #
trgt_noSky_folder <- paste0(noSky_beam_folder, "Target_noSky/")                #
if(!dir.exists(trgt_noSky_folder)){                                            #
  mkdir(trgt_noSky_folder)                                                     #
}                                                                              #
mw_noSky_folder <- paste0(noSky_beam_folder, "MW_noSky/")                      #
if(!dir.exists(mw_noSky_folder)){                                              #
  mkdir(mw_noSky_folder)                                                       #
}                                                                              #
                                                                               #
trgt_noSky_bin_folder <- paste0(trgt_noSky_folder, "Binned/", bin_fold)        #
if(!dir.exists(trgt_noSky_bin_folder)){                                        #
  mkdir(trgt_noSky_bin_folder)                                                 #
}                                                                              #
                                                                               #
trgt_noSky_stat_folder <- paste0(trgt_noSky_bin_folder, "Stats/")              #
if(!dir.exists(trgt_noSky_stat_folder)){                                       #
  mkdir(trgt_noSky_stat_folder)                                                #
}                                                                              #
                                                                               #
beam_sky_folder <- paste0(beam_folder, "Sky/")                                 #
if(!dir.exists(beam_sky_folder)){                                              #
  mkdir(beam_sky_folder)                                                       #
}                                                                              #
trgt_sky_folder <- paste0(beam_sky_folder, "Target_Sky/")                      #
if(!dir.exists(trgt_sky_folder)){                                              #
  mkdir(trgt_sky_folder)                                                       #
}                                                                              #
mw_sky_folder <- paste0(beam_sky_folder, "MW_Sky/")                            #
if(!dir.exists(mw_sky_folder)){                                                #
  mkdir(mw_sky_folder)                                                         #
}                                                                              #
                                                                               #
trgt_mask_folder <- paste0(trgt_sky_folder, "Masks/")                          #
if(!dir.exists(trgt_mask_folder)){                                             #
  mkdir(trgt_mask_folder)                                                      #
}                                                                              #
mw_mask_folder <- paste0(mw_sky_folder, "Masks/")                              #
if(!dir.exists(mw_mask_folder)){                                               #
  mkdir(mw_mask_folder)                                                        #
}                                                                              #
                                                                               #
trgt_sky_bin_folder <- paste0(trgt_sky_folder, "Binned/", bin_fold)            #
if(!dir.exists(trgt_sky_bin_folder)){                                          #
  mkdir(trgt_sky_bin_folder)                                                   #
}                                                                              #
                                                                               #
trgt_sky_stat_folder <- paste0(trgt_sky_bin_folder, "Stats/")                  #
if(!dir.exists(trgt_sky_stat_folder)){                                         #
  mkdir(trgt_sky_stat_folder)                                                  #
}                                                                              #
                                                                               #
offset_path <- paste0(og_folder, "Processed_Offsets/")                         #
if(!dir.exists(offset_path)){                                                  #
  mkdir(offset_path)                                                           #
}                                                                              #
                                                                               #
obs_Iflux_folder <- paste0(offset_path, "I_flux/Obs/")                         #
if(!dir.exists(obs_Iflux_folder)){                                             #
  mkdir(obs_Iflux_folder)                                                      #
}                                                                              #
noSky_Iflux_folder <- paste0(offset_path, "I_flux/Sky_C/")                     #
if(!dir.exists(noSky_Iflux_folder)){                                           #
  mkdir(noSky_Iflux_folder)                                                    #
}                                                                              #
sky_Iflux_folder <- paste0(offset_path, "I_flux/Sky/")                         #
if(!dir.exists(sky_Iflux_folder)){                                             #
  mkdir(sky_Iflux_folder)                                                      #
}                                                                              #
                                                                               #
pol_path <- paste0(offset_path, "Stokes/")                                     #
if(!dir.exists(pol_path)){                                                     #
  mkdir(pol_path)                                                              #
}                                                                              #
rm(og_folder)                                                                  #
                                                                               #
out_path <- paste0(pol_path, bin_fold)                                         #
print("Provide an output path for this session: ")                             #
out_path <- edit(out_path)                                                     #
if(!dir.exists(out_path)){                                                     #
  mkdir(out_path)                                                              #
}                                                                              #
                                                                               #
# Path for files without any polarization correction                           #
obsPol_path <- paste0(out_path, "Obs/")                                        #
if(!dir.exists(obsPol_path)){                                                  #
  mkdir(obsPol_path)                                                           #
}                                                                              #
                                                                               #
# Path for files related to Milky Way polarization                             #
mwPol_path <- paste0(pol_path, "MW/")                                          #
if(!dir.exists(mwPol_path)){                                                   #
  mkdir(mwPol_path)                                                            #
}                                                                              #
rm(pol_path)                                                                   #
                                                                               #
mwSrcACP_path <- paste0(mwPol_path, "Lists_of_Sources/Accepted/")              #
if(!dir.exists(mwSrcACP_path)){                                                #
  mkdir(mwSrcACP_path)                                                         #
}                                                                              #
mwSrcAPT_path <- paste0(mwSrcACP_path, "Optimal_Aperture_Hist/")               #
if(!dir.exists(mwSrcAPT_path)){                                                #
  mkdir(mwSrcAPT_path)                                                         #
}                                                                              #
mwSrcFLX_path <- paste0(mwSrcACP_path, "Fluxes/")                              #
if(!dir.exists(mwSrcFLX_path)){                                                #
  mkdir(mwSrcFLX_path)                                                         #
}                                                                              #
mwSrcPOL_path <- paste0(mwSrcACP_path, "Pol/")                                 #
if(!dir.exists(mwSrcPOL_path)){                                                #
  mkdir(mwSrcPOL_path)                                                         #
}                                                                              #
mwSrcREJ_path <- paste0(mwPol_path, "Lists_of_Sources/Rejected/")              #
if(!dir.exists(mwSrcREJ_path)){                                                #
  mkdir(mwSrcREJ_path)                                                         #
}                                                                              #
seg_path <- paste0(mwPol_path, "Segmentation_Maps/")                           #
if(!dir.exists(seg_path)){                                                     #
  mkdir(seg_path)                                                              #
}                                                                              #
mwPolCSV_path <- paste0(mwPol_path, "Pol/CSV/")                                #
if(!dir.exists(mwPolCSV_path)){                                                #
  mkdir(mwPolCSV_path)                                                         #
}                                                                              #
mwPolPDF_path <- paste0(mwPol_path, "Pol/PDF/")                                #
if(!dir.exists(mwPolPDF_path)){                                                #
  mkdir(mwPolPDF_path)                                                         #
}                                                                              #
                                                                               #
# Path for files related to background polarization                            #
skyPol_path <- paste0(out_path, "Sky/")                                        #
if(!dir.exists(skyPol_path)){                                                  #
  mkdir(skyPol_path)                                                           #
}                                                                              #
                                                                               #
# Path for files with background polarization correction                       #
sky_instCPol_path <- paste0(out_path, "Sky&Inst_C/")                           #
if(!dir.exists(sky_instCPol_path)){                                            #
  mkdir(sky_instCPol_path)                                                     #
}                                                                              #
                                                                               #
# Path for files with all polarization corrections                             #
allCPol_path <- paste0(out_path, "All_C/")                                     #
if(!dir.exists(allCPol_path)){                                                 #
  mkdir(allCPol_path)                                                          #
}                                                                              #
rm(out_path)                                                                   #
                                                                               #
pnt_src_t <- 10                                                                #
print("Define the threshold value for point source masking:")                  #
pnt_src_t <- edit(pnt_src_t)                                                   #
                                                                               #
src_sig <- 10                                                                  #
print("Define the sigma value for point source extraction:")                   #
src_sig <- edit(src_sig)                                                       #
################################################################################

##_##################################_###_###################################_##
#/ \--------------------------------/ \|/ \---------------------------------/ \#
#\_/--------------------------------\_/|\_/---------------------------------\_/#
################################################################################

################################################################################
########### Setting up files lists, confirming bands and HWP angles ############
################################################################################
print(paste0("##### Starting to process beam files of ", og_obj, " to bin siz",#
             "e ", binS, "x", binS, " #####"))                                 #
                                                                               #
# List of output pdf labels                                                    #
pdf_labels <- c("obs", 'sky & inst C', "all C", "sky")                         #
                                                                               #
# List of output FITS file keywords                                            #
file_type_1 <- c('OBSERV', 'SKY_INST C', 'ALL C', 'SKY')                       #
                                                                               #
# List of output file labels                                                   #
file_label <- c('obs', 'sky_instC', 'allC', 'sky')                             #
                                                                               #
# List of output paths                                                         #
folder_t <- c(obsPol_path, sky_instCPol_path, allCPol_path, skyPol_path)       #
rm(obsPol_path, sky_instCPol_path, allCPol_path)                               #
                                                                               #
# Listof paths to, and names of, files in the specified folder                 #
pathList <- list.files(input_folder, full.names = T)                           #
fileList <- list.files(input_folder, full.names = F)                           #
instList <- list.files(inst_folder, full.names = T)                            #
instNames <- list.files(inst_folder, full.names = F)                           #
acqList <- list.files(acq_folder, full.names = T)                              #
acqNames <- list.files(acq_folder, full.names = F)                             #
rm(input_folder, acq_folder)                                                   #
                                                                               #
if(use_inst_bin){                                                              #
  instBinList <- list.files(inst_bin_folder, full.names = T)                   #
  instBinNames <- list.files(inst_bin_folder, full.names = F)                  #
}                                                                              #
                                                                               #
# Exclusion of non-fits files from the previous lists                          #
is_fits <- grep(".fits", pathList)                                             #
is_acq <- grep("Acq-no-CRs.fits", acqList)                                     #
is_instfits <- grep(".fits", instList)                                         #
is_instBinfits <- grep(".fits", instBinList)                                   #
                                                                               #
if(length(is_fits) != 0){                                                      #
  pathList <- pathList[is_fits]                                                #
  fileList <- fileList[is_fits]                                                #
}                                                                              #
rm(is_fits)                                                                    #
                                                                               #
if(length(is_instfits) != 0){                                                  #
  instList <- instList[is_instfits]                                            #
  instNames <- instNames[is_instfits]                                          #
}                                                                              #
rm(is_instfits)                                                                #
                                                                               #
if(length(is_instBinfits) != 0){                                               #
  instBinList <- instBinList[is_instBinfits]                                   #
  instBinNames <- instBinNames[is_instBinfits]                                 #
}                                                                              #
rm(is_instBinfits)                                                             #
                                                                               #
if(length(is_acq) != 0){                                                       #
  acqList <- acqList[is_acq]                                                   #
  acqNames <- acqNames[is_acq]                                                 #
}                                                                              #
rm(is_acq)                                                                     #
                                                                               #
acqN <- length(acqList)                                                        #
                                                                               #
# This assumes the following filters were used, change if necessary:           ###
# "u_HIGH", "b_HIGH", "v_HIGH", "R_SPECIAL", "I_BESS"                          #
band_str <- c("u_HIGH", "b_HIGH", "v_HIGH", "R_SPECIAL", "I_BESS")             #
nB <- length(band_str)                                                         #
bands <- 1:nB                                                                  #
                                                                               #
inst_ind <- array(NA, dim = c(2, nB, 2),                                       #
                  dimnames = list(c("Q", "U"), band_str, c("Data", "Unc")))    #
                                                                               #
Q_ind <- grep("_inlaQ", instList)                                              #
U_ind <- grep("_inlaU", instList)                                              #
eQ_ind <- grep("_erinlaQ", instList)                                           #
eU_ind <- grep("_erinlaU", instList)                                           #
                                                                               #
if(use_inst_bin){                                                              #
  inst_bin_ind <- array(NA, dim = c(2, nB),                                    #
                        dimnames = list(c("Q","U"), band_str))                 #
                                                                               #
  Qb_ind <- grep("_inlaQ", instBinList)                                        #
  Ub_ind <- grep("_inlaU", instBinList)                                        #
}                                                                              #
                                                                               #
for(b in bands[-1]){                                                           #
  bs <- band_str[b]                                                            #
                                                                               #
  instband_ind <- grep(bs, instList)                                           #
                                                                               #
  inst_ind["Q", b, "Data"] <- intersect(instband_ind, Q_ind)                   #
  inst_ind["U", b, "Data"] <- intersect(instband_ind, U_ind)                   #
  inst_ind["Q", b, "Unc"] <- intersect(instband_ind, eQ_ind)                   #
  inst_ind["U", b, "Unc"] <- intersect(instband_ind, eU_ind)                   #
                                                                               #
  if(use_inst_bin){                                                            #
    instband_ind <- grep(bs, instBinList)                                      #
                                                                               #
    inst_bin_ind["Q", b] <- intersect(instband_ind, Qb_ind)                    #
    inst_bin_ind["U", b] <- intersect(instband_ind, Ub_ind)                    #
  }                                                                            #
}                                                                              #
rm(Q_ind, U_ind, instband_ind, b, bs, eQ_ind, eU_ind)                          #
                                                                               #
if(use_inst_bin){                                                              #
  rm(Qb_ind, Ub_ind)                                                           #
}                                                                              #
                                                                               #
fitsCount <- length(fileList)                                                  #
                                                                               #
# Keywords of interest in fits files' headers                                  #
REFX <- "CRPIX1"                                                               #
REFY <- "CRPIX2"                                                               #
VALX <- "CRVAL1"                                                               #
VALY <- "CRVAL2"                                                               #
TYPEX <- "CTYPE1"                                                              #
TYPEY <- "CTYPE2"                                                              #
DATE <- "DATE-OBS"                                                             #
FILTER <- "HIERARCH ESO INS FILT1 NAME"                                        #
HWP_ANG <- "HIERARCH ESO INS RETA2 POSANG"                                     #
AUTHOR <- "AUTHOR"                                                             #
TYP1 <- "FILETYP1"                                                             #
TYP2 <- "FILETYP2"                                                             #
TARGET <- "OBJECT  ="                                                          #
REFRA <- "RA     "                                                             #
REFDEC <- "DEC    "                                                            #
PIXSCALE <- "HIERARCH ESO INS PIXSCALE"                                        #
BINX <- "HIERARCH ESO DET WIN1 BINX"                                           #
BINY <- "HIERARCH ESO DET WIN1 BINY"                                           #
GAIN <- "HIERARCH ESO DET OUT1 GAIN"                                           #
                                                                               #
# Determine type of fits files present                                         #
# Possible types of fits files:                                                #
# "PRODUCT", "UNCERT", "NOT EXP"                                               #
types_arr <- get_fits_header_list_str(pathList, TYP1)                          #
# Determine type of fits files present                                         #
# Possible types of fits files:                                                #
# "E_BEAM", "O_BEAM", "NOT BEAM"                                               #
beams_arr <- get_fits_header_list_str(pathList, TYP2)                          #
# Determine type of fits files present                                         #
# Possible types of fits files:                                                #
# "JRS2023", "NOT EXP"                                                         #
auth_arr <- get_fits_header_list_str(pathList, AUTHOR)                         #
rm(TYP1, TYP2, AUTHOR)                                                         #
                                                                               #
print("Checking for unexpected type files...")                                 #
                                                                               #
dimsObj <- NULL                                                                #
                                                                               #
for(f in 1:fitsCount){                                                         #
                                                                               #
  if((types_arr[f] != "PRODUCT" && types_arr[f] != "UNCERT")){                 #
    types_arr[f] <- "NOT EXP"                                                  #
  }                                                                            #
  if((beams_arr[f] != "O_BEAM" && beams_arr[f] != "E_BEAM")){                  #
    beams_arr[f] <- "NOT BEAM"                                                 #
  }                                                                            #
                                                                               #
  if(auth_arr[f] != "JRS2023"){                                                #
    auth_arr[f] <- "NOT EXP"                                                   #
  }                                                                            #
                                                                               #
  if(is.null(dimsObj) && types_arr[f] != "NOT EXP" &&                          #
     beams_arr[f] != "NOT BEAM" && auth_arr[f] != "NOT EXP"){                  #
    dimsObj <- dim(readFITS(pathList[f])$imDat)[1:2]                           #
    delt <- dimsObj[1] * .8 / dimsObj[2]                                       #
  }                                                                            #
}                                                                              #
rm(fitsCount, f)                                                               #
                                                                               #
print("Indexing beam files...")                                                #
                                                                               #
# Create list of pathList indexes for E_BEAMs fits                             #
E_ind <- as.vector(which(beams_arr == "E_BEAM" & auth_arr == "JRS2023" &       #
                           types_arr == "PRODUCT", arr.ind = TRUE))            #
eCount <- length(E_ind)                                                        #
O_ind <- as.vector(which(beams_arr == "O_BEAM" & auth_arr == "JRS2023" &       #
                           types_arr == "PRODUCT", arr.ind = TRUE))            #
oCount <- length(O_ind)                                                        #
rm(beams_arr, auth_arr, types_arr)                                             #
                                                                               #
print("Checking if beam data and uncertainty files are paired...")             #
                                                                               #
err_no_pair <- "At least one file is not paired."                              #
stopifnot(err_no_pair = eCount == oCount)                                      #
print("E and O beams are properly grouped.")                                   #
tCount <- eCount                                                               #
rm(eCount, oCount, err_no_pair)                                                #
                                                                               #
# Determine bands in which beams are available                                 #
bands <- 1:nB                                                                  #
ang_str <- c("0", "22.5", "45", "67.5")                                        #
angs <- 1:length(ang_str)                                                      #
                                                                               #
print("Tagging beam files by band filter...")                                  #
tagDim <- c(2, tCount)                                                         #
tagDimNames <- list(c("E", "O"), NULL)                                         #
                                                                               #
bands_arr <- array(NA, dim = tagDim, dimnames = tagDimNames)                   #
bands_arr["E",] <- get_fits_header_list_str(pathList[E_ind], FILTER)           #
bands_arr["O",] <- get_fits_header_list_str(pathList[O_ind], FILTER)           #
rm(FILTER)                                                                     #
                                                                               #
print("Tagging beam files by HWP angle...")                                    #
angs_num_arr <- array(NA, dim = tagDim, dimnames = tagDimNames)                #
angs_num_arr["E",] <- get_fits_header_list_num(pathList[E_ind], HWP_ANG)       #
angs_num_arr["O",] <- get_fits_header_list_num(pathList[O_ind], HWP_ANG)       #
rm(HWP_ANG)                                                                    #
                                                                               #
angs_str_arr <- array(NA, dim = tagDim, dimnames = tagDimNames)                #
numbs_arr <- array(0, dim = tagDim, dimnames = tagDimNames)                    #
rm(tagDim, tagDimNames)                                                        #
                                                                               #
# To assure there are 4 pol ang per band                                       #
ang_count <- array(0, dim = c(2, 4), dimnames = list(c("E", "O"), ang_str))    #
                                                                               #
print(paste0("Listing beam data and uncertainty files while checking they wer",#
             "e measured at four different HWP angles..."))                    #
for(t in 1:tCount){                                                            #
  name_e <- fileList[E_ind[t]]                                                 #
  name_o <- fileList[O_ind[t]]                                                 #
                                                                               #
  numbs_arr["E", t] <- as.numeric(str_extract(str_extract(name_e, "#[0-9]+"),  #
                                              "[0-9]+"))                       #
  numbs_arr["O", t] <- as.numeric(str_extract(str_extract(name_o, "#[0-9]+"),  #
                                              "[0-9]+"))                       #
                                                                               #
  if(-2 < angs_num_arr["E", t] && angs_num_arr["E", t] < 2){                   #
    ang_count["E", "0"] <- ang_count["E", "0"] + 1                             #
    angs_str_arr["E", t] <- ang_str[1]                                         #
  }                                                                            #
  if((22.5 - 2) < angs_num_arr["E", t] && angs_num_arr["E", t] < (22.5 + 2)){  #
    ang_count["E", "22.5"] <- ang_count["E", "22.5"] + 1                       #
    angs_str_arr["E", t] <- ang_str[2]                                         #
  }                                                                            #
  if((45 - 2) < angs_num_arr["E", t] && angs_num_arr["E", t] < (45 + 2)){      #
    ang_count["E", "45"] <- ang_count["E", "45"] + 1                           #
    angs_str_arr["E", t] <- ang_str[3]                                         #
  }                                                                            #
  if((67.5 - 2) < angs_num_arr["E", t] && angs_num_arr["E", t] < (67.5 + 2)){  #
    ang_count["E", "67.5"] <- ang_count["E", "67.5"] + 1                       #
    angs_str_arr["E", t] <- ang_str[4]                                         #
  }                                                                            #
                                                                               #
  if(-2 < angs_num_arr["O", t] && angs_num_arr["O", t] < 2){                   #
    ang_count["O", "0"] <- ang_count["O", "0"] + 1                             #
    angs_str_arr["O", t] <- ang_str[1]                                         #
  }                                                                            #
  if((22.5 - 2) < angs_num_arr["O", t] && angs_num_arr["O", t] < (22.5 + 2)){  #
    ang_count["O", "22.5"] <- ang_count["O", "22.5"] + 1                       #
    angs_str_arr["O", t] <- ang_str[2]                                         #
  }                                                                            #
  if((45 - 2) < angs_num_arr["O", t] && angs_num_arr["O", t] < (45 + 2)){      #
    ang_count["O", "45"] <- ang_count["O", "45"] + 1                           #
    angs_str_arr["O", t] <- ang_str[3]                                         #
  }                                                                            #
  if((67.5 - 2) < angs_num_arr["O", t] && angs_num_arr["O", t] < (67.5 + 2)){  #
    ang_count["O", "67.5"] <- ang_count["O", "67.5"] + 1                       #
    angs_str_arr["O", t] <- ang_str[4]                                         #
  }                                                                            #
                                                                               #
  print(paste0("O beam #", t, " is in ", bands_arr["O", t], " band, observe",  #
               "d with the HWP at an angle of ", angs_str_arr["O", t], "º."))  #
  print(paste0("E beam #", t, " is in ", bands_arr["E", t], " band, observe",  #
               "d with the HWP at an angle of ", angs_str_arr["E", t], "º."))  #
                                                                               #
  if(t == tCount){                                                             #
                                                                               #
    error_HWP_1 <- paste0("ERROR: the number of files at each HWP angle is ei",#
                          "ther not the same, or zero. There are ",            #
                          ang_count["O", "0"], " files at 0º, ",               #
                          ang_count["O", "22.5"], " files at 22.5º, ",         #
                          ang_count["O", "45"], " files at 45º and ",          #
                          ang_count["O", "67.5"], " files at 67.5º for O beam",#
                          ". And there are ", ang_count["E", "0"], " files at",#
                          " 0º, ", ang_count["E", "22.5"], " files at 22.5º, ",#
                          ang_count["E", "45"], " files at 45º and ",          #
                          ang_count["E", "67.5"], " files at 67.5º for E beam",#
                          ".")                                                 #
    stopifnot(error_HWP_1  = ang_count["O", "0"] == ang_count["O", "22.5"],    #
              ang_count["O", "0"] == ang_count["O", "45"],                     #
              ang_count["O", "0"] == ang_count["O", "67.5"],                   #
              ang_count["O", "0"] == ang_count["E", "0"],                      #
              ang_count["O", "0"] == ang_count["E", "22.5"],                   #
              ang_count["O", "0"] == ang_count["E", "45"],                     #
              ang_count["O", "0"] == ang_count["E", "67.5"],                   #
              ang_count["O", "0"] != 0)                                        #
    rm(error_HWP_1)                                                            #
  }                                                                            #
  rm(name_e, name_o)                                                           #
}                                                                              #
rm(t, tCount, ang_count)                                                       #
                                                                               #
print("Creating O and E beam data files band filter indexes...")               #
max_N <- array(0, dim = nB, dimnames = list(band_str))                         #
for(b in bands){                                                               #
  bs <- band_str[b]                                                            #
                                                                               #
  t_O_ind <- which(bands_arr["O",] == bs, arr.ind = TRUE)                      #
  t_E_ind <- which(bands_arr["E",] == bs, arr.ind = TRUE)                      #
                                                                               #
  switch(b,                                                                    #
         c(u_ind_O <- t_O_ind, u_ind_E <- t_E_ind),                            #
         c(b_ind_O <- t_O_ind, b_ind_E <- t_E_ind),                            #
         c(v_ind_O <- t_O_ind, v_ind_E <- t_E_ind),                            #
         c(r_ind_O <- t_O_ind, r_ind_E <- t_E_ind),                            #
         c(i_ind_O <- t_O_ind, i_ind_E <- t_E_ind)                             #
  )                                                                            #
                                                                               #
  print(paste0("Checking how many measurement sets there are in band ", bs,    #
               "..."))                                                         #
  max_N[bs] <- max(numbs_arr["O", t_O_ind], na.rm=T)                           #
                                                                               #
  if(is.infinite(max_N[bs]) || max_N[bs] == 0){                                #
    max_N[bs] <- 0                                                             #
    bands <- bands[-which(bands == b, arr.ind = T)]                            #
                                                                               #
    switch(b,                                                                  #
           rm(u_ind_O, u_ind_E), rm(b_ind_O, b_ind_E), rm(v_ind_O, v_ind_E),   #
           rm(r_ind_O, r_ind_E), rm(i_ind_O, i_ind_E)                          #
    )                                                                          #
  }                                                                            #
}                                                                              #
rm(b, bands_arr)                                                               #
                                                                               #
print("Creating O and E beam data files HWP angle indexes...")                 #
for(ang in angs){                                                              #
  as <- ang_str[ang]                                                           #
                                                                               #
  t_O_ind <- which(angs_str_arr["O",] == as, arr.ind = TRUE)                   #
  t_E_ind <- which(angs_str_arr["E",] == as, arr.ind = TRUE)                   #
                                                                               #
  switch(ang,                                                                  #
         c(ind_0_O <- t_O_ind, ind_0_E <- t_E_ind),                            #
         c(ind_22_O <- t_O_ind, ind_22_E <- t_E_ind),                          #
         c(ind_45_O <- t_O_ind, ind_45_E <- t_E_ind),                          #
         c(ind_68_O <- t_O_ind, ind_68_E <- t_E_ind)                           #
  )                                                                            #
}                                                                              #
rm(ang, t_O_ind, t_E_ind, angs_num_arr, angs_str_arr)                          #
################################################################################

##_##################################_###_###################################_##
#/ \--------------------------------/ \|/ \---------------------------------/ \#
#\_/--------------------------------\_/|\_/---------------------------------\_/#
################################################################################

################################################################################
#################### Setting up Python required libraries ######################
################################################################################
print("Loading relevant python libraries...")                                  #
py_run_string("import numpy as np")                                            #
py_run_string("import sep")                                                    #
py_run_string("import scipy")                                                  #
py_run_string("from photutils import MedianBackground")                        #
py_run_string("from photutils.segmentation import make_source_mask")           #
py_run_string("from photutils.isophote import EllipseGeometry")                #
py_run_string("from photutils.isophote import Ellipse")                        #
py_run_string("from astropy.io import fits")                                   #
py_run_string("from astropy.stats import SigmaClip")                           #
py_run_string("sig_Clip = SigmaClip(sigma=3.0)")                               #
py_run_string("bkg_Est = MedianBackground()")                                  #
py_run_string("sep.set_extract_pixstack(600000)")                              #
rm(libs_folder)                                                                #
                                                                               #
sat <- 2^16 * .9                                                               #
max_ecc <- 0.6                                                                 #
max_A <- max(c(binS^2, 30^2))                                                  #
min_A <- 9                                                                     #
rm(load_mod_cmd, mypy_lib)                                                     #
################################################################################

################################################################################
#-_----------------------------------_---_-----------------------------------_-#
#/ \--------------------------------/ \|/ \---------------------------------/ \#
#\_/--------------------------------\_/|\_/---------------------------------\_/#
#------------------------------------------------------------------------------#
################################################################################

################################################################################
################################ BEAM PROCESSOR ################################
################################################################################
run_date <-  paste0("_", gsub(" ", "_", Sys.time()))                           #
                                                                               #
singDim <- c(dimsObj, 2)                                                       #
multDim <- c(dimsObj, 2, 4)                                                    #
pxDim <- c(dimsObj, 3)                                                         #
singNames <- list(NULL, NULL, c("Data", "Unc"))                                #
multNames <- list(NULL, NULL, c("Data", "Unc"), ang_str)                       #
pxNames <- list(NULL, NULL, c("Data", "Unc", "Unc(%)"))                        #
                                                                               #
# From FORS2 User Manual Tab. 4.7, values in deg                               #
X_err <- array(NA, dim = nB, dimnames = list(band_str))                        #
X_err["u_HIGH"] <- -2.07 * pi / 180                                            #
X_err["b_HIGH"] <- +1.54 * pi / 180                                            #
X_err["v_HIGH"] <- +1.80 * pi / 180                                            #
X_err["R_SPECIAL"] <- -1.19 * pi / 180                                         #
X_err["I_BESS"] <- -2.89 * pi / 180                                            #
uX_err <- 0.005 * pi / 180                                                     #
                                                                               #
N_max <- max(max_N, na.rm = T)                                                 #
                                                                               #
apt_R_flag <- file.exists(paste0(mwSrcPOL_path, "optimal_apt_radius_tab.csv")) #
                                                                               #
if(!apt_R_flag){                                                               #
  apt_rad_arr <- array(NA, dim = c(nB, N_max),                                 #
                       dimnames = list(band_str, paste0("#", 1:N_max)))        #
}                                                                              #
                                                                               #
for(b in bands){                                                               #
  bs <- band_str[b]                                                            #
  print(paste0("Starting to process files in band ", bs, "..."))               #
                                                                               #
  print("Loading relevant band filter file indexes...")                        #
  switch(b,                                                                    #
         c(band_ind_E <- u_ind_E, band_ind_O <- u_ind_O,                       #
           rm(u_ind_E, u_ind_O)),                                              #
         c(band_ind_E <- b_ind_E, band_ind_O <- b_ind_O,                       #
           rm(b_ind_E, b_ind_O)),                                              #
         c(band_ind_E <- v_ind_E, band_ind_O <- v_ind_O,                       #
           rm(v_ind_E, v_ind_O)),                                              #
         c(band_ind_E <- r_ind_E, band_ind_O <- r_ind_O,                       #
           rm(r_ind_E, r_ind_O)),                                              #
         c(band_ind_E <- i_ind_E, band_ind_O <- i_ind_O,                       #
           rm(i_ind_E, i_ind_O))                                               #
  )                                                                            #
                                                                               #
  # Polarization degree instrumental corrections to be applied later on        #
  c2X <- cos(2 * X_err[b])                                                     #
  s2X <- sin(2 * X_err[b])                                                     #
  u_c2X <- sqrt((2 * s2X * uX_err)^2)                                          #
  u_s2X <- sqrt((2 * c2X * uX_err)^2)                                          #
                                                                               #
  # Load the instrumental polarization correction maps                         #
  print("Loading instrumental polarization correction maps...")                #
                                                                               #
  Qinst <- array(NA, dim = singDim, dimnames = singNames)                      #
  Uinst <- array(NA, dim = singDim, dimnames = singNames)                      #
                                                                               #
  if(b > 1){                                                                   #
    instQ <- readFITS(instList[inst_ind["Q", b, "Data"]])$imDat                #
    instU <- readFITS(instList[inst_ind["U", b, "Data"]])$imDat                #
    unc_instQ <- readFITS(instList[inst_ind["Q", b, "Unc"]])$imDat             #
    unc_instU <- readFITS(instList[inst_ind["U", b, "Unc"]])$imDat             #
  }                                                                            #
                                                                               #
  ############################################################################ #
  ########## Fixing instrumental correction Q and U map dimensions ########### #
  #### Assumes all instrumental correction maps have the same dimensions ##### #
  ##############################################################################
  if(b > 1){                                                                 # #
    if(!instV){                                                              # #
      print("Adjusting instrumental correction files dimensions...")         # #
                                                                             # #
      corr_dims <- dim(instQ)                                                # #
      dif_dims <- dimsObj - corr_dims                                        # #
      rm(corr_dims)                                                          # #
                                                                             # #
      dif_x <- dif_dims[1] %/% 2                                             # #
      x_ini <- 1 + dif_x                                                     # #
      if(2 * dif_x < dif_dims[1]){                                           # #
        dif_x <- dif_x + 1                                                   # #
      }                                                                      # #
      x_fin <- dimsObj[1] - dif_x                                            # #
      x_CD <- x_ini:x_fin                                                    # #
      rm(x_fin, x_ini, dif_x)                                                # #
                                                                             # #
      dif_y <- dif_dims[2] %/% 2                                             # #
      y_ini <- 1 + dif_y                                                     # #
      if(2 * dif_y < dif_dims[2]){                                           # #
        dif_y <- dif_y + 1                                                   # #
      }                                                                      # #
      y_fin <- dimsObj[2] - dif_y                                            # #
      y_CD <- y_ini:y_fin                                                    # #
      rm(y_fin, y_ini, dif_y)                                                # #
                                                                             # #
      Qinst[x_CD, y_CD, "Data"] <- instQ                                     # #
      Qinst[x_CD, y_CD, "Unc"] <- unc_instQ                                  # #
      Uinst[x_CD, y_CD, "Data"] <- instU                                     # #
      Uinst[x_CD, y_CD, "Unc"] <- unc_instU                                  # #
      rm(dif_dims, x_CD, y_CD, instQ, unc_instQ, instU, unc_instU)           # #
                                                                             # #
      if(save_JRS2023){                                                      # #
        print(paste0("* Saving FORS2 polarization maps version JRS2023, for",# #
                     " band ", bs, " *"))                                    # #
                                                                             # #
        writeFITSim(Qinst[,,"Data"],                                         # #
                    paste0(inst_new_folder,                                  # #
                           instNames[inst_ind["Q", b, "Data"]]))             # #
        writeFITSim(Uinst[,,"Data"],                                         # #
                    paste0(inst_new_folder,                                  # #
                           instNames[inst_ind["U", b, "Data"]]))             # #
        writeFITSim(Qinst[,,"Unc"],                                          # #
                    paste0(inst_new_folder,                                  # #
                           instNames[inst_ind["Q", b, "Unc"]]))              # #
        writeFITSim(Uinst[,,"Unc"],                                          # #
                    paste0(inst_new_folder,                                  # #
                           instNames[inst_ind["U", b, "Unc"]]))              # #
      }                                                                      # #
    }else{                                                                   # #
        Qinst[,,"Data"] <- instQ                                             # #
        Qinst[,, "Unc"] <- unc_instQ                                         # #
        Uinst[,, "Data"] <- instU                                            # #
        Uinst[,, "Unc"] <- unc_instU                                         # #
        rm(instQ, unc_instQ, instU, unc_instU)                               # #
    }                                                                        # #
  }else{                                                                     # #
    Qinst[,,] <- 0                                                           # #
    Uinst[,,] <- 0                                                           # #
  }                                                                          # #
  ##############################################################################
                                                                               #
  #--------------------------------------------------------------------------# #
                                                                               #
  ######################### Binning Instrumental Data ######################## #
  ##############################################################################
  if(!use_inst_bin){                                                         # #
    if(binS > 1){                                                            # #
      print(paste0("* Binning instrumental Q and U maps pixels by their med",# #
                   "ian within ", binS, "x", binS, " pixel windows, of band",# #
                   " ", bs, " *"))                                           # #
                                                                             # #
      cr_px <- c(get_fits_header_num(readFITS(acqList[1])$header, REFX),     # #
                 get_fits_header_num(readFITS(acqList[1])$header, REFY))     # #
                                                                             # #
      bQinst <- bin_code_map(Qinst, binS, cr_px)                             # #
      bUinst <- bin_code_map(Uinst, binS, cr_px)                             # #
      rm(cr_px)                                                              # #
                                                                             # #
      binDim <- dim(bQinst)                                                  # #
                                                                             # #
      if(save_inst_bin){                                                     # #
        print(paste0("* Saving FORS2 polarization maps binned by ", binS,    # #
                     "x", binS, " pixel windows, of band ", bs, " *"))       # #
                                                                             # #
        tmp_split <- strsplit(instNames[inst_ind["Q", b, "Data"]], ".fits")  # #
        tmp_name <- paste0(tmp_split[[1]][1], "_bin=", binS, ".fits")        # #
        writeFITSim(bQinst, paste0(inst_bin_folder, tmp_name))               # #
                                                                             # #
        tmp_split <- strsplit(instNames[inst_ind["U", b, "Data"]], ".fits")  # #
        tmp_name <- paste0(tmp_split[[1]][1], "_bin=", binS, ".fits")        # #
        writeFITSim(bUinst, paste0(inst_bin_folder, tmp_name))               # #
        rm(tmp_name, tmp_split)                                              # #
      }                                                                      # #
    }else{                                                                   # #
      binDim <- c(dimsObj, 2)                                                # #
    }                                                                        # #
  }else{                                                                     # #
    binDim <- dim(readFITS(instBinList[1])$imDat)                            # #
    bQinst <- array(0, dim = binDim, dimnames = singNames)                   # #
    bUinst <- array(0, dim = binDim, dimnames = singNames)                   # #
                                                                             # #
    if(b > 1){                                                               # #
      temp <- readFITS(instBinList[inst_bin_ind["Q",b]])$imDat               # #
      bQinst[,,"Data"] <- temp[,,1]                                          # #
      bQinst[,,"Unc"] <- temp[,,2]                                           # #
                                                                             # #
      temp <- readFITS(instBinList[inst_bin_ind["U",b]])$imDat               # #
      bUinst[,,"Data"] <- temp[,,1]                                          # #
      bUinst[,,"Unc"] <- temp[,,2]                                           # #
      rm(temp)                                                               # #
    }                                                                        # #
  }                                                                          # #
                                                                             # #
  bin_dims <- c(binDim, 4)                                                   # #
  ##############################################################################
                                                                               #
  N <- max_N[bs]                                                               #
                                                                               #
  # Iterating over the different offset runs at the same band filter           #
  for(n in 1:N){                                                               #
    print(paste0("*** Starting to process file set #", n, " of band ", bs,     #
                 " ***"))                                                      #
                                                                               #
    print("* Loading relevant set number file indexes *")                      #
    numb_ind_E <- which(numbs_arr["E",] == n, arr.ind = TRUE)                  #
    numb_ind_O <- which(numbs_arr["O",] == n, arr.ind = TRUE)                  #
                                                                               #
    # Intersect indexes of Band with Set_Numb                                  #
    print(paste0("* Getting relevant data files indexes from intersection of ",#
                 "set number indexes with band filter indexes *"))             #
    temp_ind_E <- intersect(band_ind_E, numb_ind_E)                            #
    length_ind_E <- length(temp_ind_E)                                         #
    temp_ind_O <- intersect(band_ind_O, numb_ind_O)                            #
    length_ind_O <- length(temp_ind_O)                                         #
                                                                               #
    print(paste0("* Checking that are 4 different E and O beam data and uncer",#
                 "tainty files *"))                                            #
    if(length_ind_E != 4 || length_ind_O != 4){                                #
      print(paste0("#SKIP! Expected 8 data files, 4 on O beam and 4 on E beam",#
                   ", ", length_ind_O, " were found for O beam, and ",         #
                   length_ind_E, "were found for E beam instead for the ",     #
                   "filter & number subset: band = ", bs, ", number = ", n))   #
      next                                                                     #
    }                                                                          #
    rm(length_ind_E, length_ind_O)                                             #
                                                                               #
    # Holder arrays                                                            #
    Es <- array(NA, dim = multDim, dimnames = multNames)                       #
    Os <- array(NA, dim = multDim, dimnames = multNames)                       #
                                                                               #
    ######################### Loading Relevant Files ######################### #
    ############################################################################
    # Intersect Beam indexes with HWP angle indexes                          # #
    # Save file name to later fill header                                    # #
    # Load file to holder array                                              # #
    print(paste0("* Loading E beam data from set #", n, " of band ", bs,     # #
                 " *"))                                                      # #
                                                                             # #
    i <- intersect(temp_ind_E, ind_0_E)                                      # #
    or1 <- fileList[E_ind[i]]                                                # #
    Es[,,,ang_str[1]] <- readFITS(pathList[E_ind[i]])$imDat                  # #
                                                                             # #
    i <- intersect(temp_ind_E, ind_22_E)                                     # #
    or2 <- fileList[E_ind[i]]                                                # #
    Es[,,,ang_str[2]] <- readFITS(pathList[E_ind[i]])$imDat                  # #
                                                                             # #
    i <- intersect(temp_ind_E, ind_45_E)                                     # #
    or3 <- fileList[E_ind[i]]                                                # #
    Es[,,,ang_str[3]] <- readFITS(pathList[E_ind[i]])$imDat                  # #
                                                                             # #
    i <- intersect(temp_ind_E, ind_68_E)                                     # #
    or4 <- fileList[E_ind[i]]                                                # #
    Es[,,,ang_str[4]] <- readFITS(pathList[E_ind[i]])$imDat                  # #
    rm(temp_ind_E)                                                           # #
                                                                             # #
    print(paste0("* Loading O beam data from set #", n, " of band ", bs,     # #
                 " *"))                                                      # #
    i <- intersect(temp_ind_O, ind_0_O)                                      # #
    or5 <- fileList[O_ind[i]]                                                # #
    Os[,,,ang_str[1]] <- readFITS(pathList[O_ind[i]])$imDat                  # #
                                                                             # #
    i <- intersect(temp_ind_O, ind_22_O)                                     # #
    or6 <- fileList[O_ind[i]]                                                # #
    Os[,,,ang_str[2]] <- readFITS(pathList[O_ind[i]])$imDat                  # #
                                                                             # #
    i <- intersect(temp_ind_O, ind_45_O)                                     # #
    or7 <- fileList[O_ind[i]]                                                # #
    Os[,,,ang_str[3]] <- readFITS(pathList[O_ind[i]])$imDat                  # #
                                                                             # #
    i <- intersect(temp_ind_O, ind_68_O)                                     # #
    or8 <- fileList[O_ind[i]]                                                # #
    Os[,,,ang_str[4]] <- readFITS(pathList[O_ind[i]])$imDat                  # #
    rm(temp_ind_O)                                                           # #
    ############################################################################
                                                                               #
    #------------------------------------------------------------------------# #
                                                                               #
    ############## Finding and Saving Coords of Stripes Borders ############## #
    ############################################################################
                                                                             # #
    # Finding rows and cols mostly filled with NA pixels                     # #
    # to determine stripe boundaries at each band                            # #
    E_NArow_bol <- apply(Es[,,"Data",], 2:3,                                 # #
                         function(x)                                         # #
                           (length(which(is.na(x))) >= (dim(Es)[1] * .5)))   # #
    O_NArow_bol <- apply(Os[,,"Data",], 2:3,                                 # #
                         function(x)                                         # #
                           (length(which(is.na(x))) >= (dim(Os)[1] * .5)))   # #
    E_NAcol_bol <- apply(Es[,,"Data",], c(1,3),                              # #
                         function(x)                                         # #
                           (length(which(is.na(x))) >= (dim(Es)[2] * .75)))  # #
    O_NAcol_bol <- apply(Os[,,"Data",], c(1,3),                              # #
                         function(x)                                         # #
                           (length(which(is.na(x))) >= (dim(Os)[2] * .75)))  # #
                                                                             # #
    E_y_lims <- apply(E_NArow_bol, 2, find_0to1_and_1to0)                    # #
    O_y_lims <- apply(O_NArow_bol, 2, find_0to1_and_1to0)                    # #
    E_x_lims <- apply(E_NAcol_bol, 2, find_0to1_and_1to0)                    # #
    O_x_lims <- apply(O_NAcol_bol, 2, find_0to1_and_1to0)                    # #
    rm(E_NArow_bol, O_NArow_bol, E_NAcol_bol, O_NAcol_bol)                   # #
    ############################################################################
                                                                               #
    #------------------------------------------------------------------------# #
                                                                               #
    ##################### Setting Up Output FITS Headers ##################### #
    ############################################################################
    print(paste0("* Preparing output data files' headers for set #", n, " o",# #
                 "f band ", bs, " *"))                                       # #
                                                                             # #
    # Extract header information                                             # #
    temp_fits <- readFITS(pathList[O_ind[i]])                                # #
    temp_head <- temp_fits$header[-(8:10)]                                   # #
    temp_head <- temp_head[-grep("WOLL POSANG",temp_head)]                   # #
    temp_head <- temp_head[-grep("RETA2 ROT",temp_head)]                     # #
    temp_head <- temp_head[-grep("RETA2 POSANG",temp_head)]                  # #
    rm(temp_fits)                                                            # #
                                                                             # #
    # Preserve attribute order across the different file types               # #
    # generated by this script series                                        # #
    temp_head <- delKwv('CAL-BIAS', temp_head)                               # #
    temp_head <- delKwv('CAL-FLAT', temp_head)                               # #
    temp_head <- delKwv('CAL-COSM', temp_head)                               # #
    temp_head <- delKwv('COSM-SIG', temp_head)                               # #
    temp_head <- delKwv('COSM-CTR', temp_head)                               # #
    temp_head <- delKwv('FILETYP1', temp_head)                               # #
    temp_head <- delKwv('FILETYP2', temp_head)                               # #
    temp_head <- delKwv('SOURCE', temp_head)                                 # #
    temp_head <- delKwv('AUTHOR', temp_head)                                 # #
    temp_head <- delKwv('ORIGIN_1', temp_head)                               # #
    temp_head <- delKwv('ORIGIN_2', temp_head)                               # #
    temp_head <- modVal('EXTNAME', 'OFFSET', "Extension name", temp_head)    # #
    temp_head <- addKwv('SRC_1', or1, "E at 0º", temp_head)                  # #
    temp_head <- addKwv('SRC_2', or2, "E at 22.5º", temp_head)               # #
    temp_head <- addKwv('SRC_3', or3, "E at 45º", temp_head)                 # #
    temp_head <- addKwv('SRC_4', or4, "E at 67.5º", temp_head)               # #
                                                                             # #
    skyE_head <- temp_head                                                   # #
                                                                             # #
    temp_head <- addKwv('SRC_5', or5, "O at 0º", temp_head)                  # #
    temp_head <- addKwv('SRC_6', or6, "O at 22.5º", temp_head)               # #
    temp_head <- addKwv('SRC_7', or7, "O at 45º", temp_head)                 # #
    temp_head <- addKwv('SRC_8', or8, "O at 67.5º", temp_head)               # #
    temp_head <- addKwv('FILETYP1', 'DATA & UNC', "File type", temp_head)    # #
    temp_head <- addKwv('FILETYP2', 'STOKES Q', "File sub-type", temp_head)  # #
    temp_head <- addKwv('METHOD', 'BAGNULO',                                 # #
                        "Method used to calc. Stoke param.", temp_head)      # #
    sky_head <- temp_head                                                    # #
    temp_head <- addKwv('CORR-SKY', 'NO', "Was sky pol. subtracted?",        # #
                        temp_head)                                           # #
    temp_head <- addKwv('CORR-MW', 'NO', "Was MW pol. subtracted?",          # #
                        temp_head)                                           # #
    sky_head <- addKwv('CORR-MW', 'YES', "Was MW pol. subtracted?",          # #
                       sky_head)                                             # #
    temp_head <- addKwv('CORR-INST', 'NO', "Was instr. pol. subtracted?",    # #
                        temp_head)                                           # #
    sky_head <- addKwv('CORR-INST', 'YES', "Was instr. pol. subtracted?",    # #
                       sky_head)                                             # #
    temp_head <- addKwv('CORR-CHRM', 'NO', "Was pol. angle corrected?",      # #
                        temp_head)                                           # #
    sky_head <- addKwv('CORR-CHRM', 'YES', "Was pol. angle corrected?",      # #
                       sky_head)                                             # #
    temp_head <- addKwv('CORR-PBIAS', 'NO', "Was pol. bias corrected?",      # #
                        temp_head)                                           # #
    sky_head <- addKwv('CORR-PBIAS', 'NO', "Was pol. bias corrected?",       # #
                       sky_head)                                             # #
    temp_head <- addKwv('BINNING', toString(binS), "Bin length", temp_head)  # #
    temp_head <- addKwv('SOURCE', 'process_Beam.R',                          # #
                        "Script used to generate this file", temp_head)      # #
    temp_head <- addKwv('AUTHOR', 'JRS2023', "Script creator", temp_head)    # #
                                                                             # #
    skyE_head <- addKwv('SOURCE', 'process_Beam.R',                          # #
                        "Script used to generate this file", skyE_head)      # #
    skyE_head <- addKwv('AUTHOR', 'JRS2023', "Script creator", skyE_head)    # #
                                                                             # #
    skyO_head <- skyE_head                                                   # #
    skyO_head <- modVal('SRC_1', or5, "O at 0º", skyO_head)                  # #
    skyO_head <- modVal('SRC_2', or6, "O at 22.5º", skyO_head)               # #
    skyO_head <- modVal('SRC_3', or7, "O at 45º", skyO_head)                 # #
    skyO_head <- modVal('SRC_4', or8, "O at 67.5º", skyO_head)               # #
    rm(or1, or2, or3, or4, or5, or6, or7, or8)                               # #
                                                                             # #
    sky_head <- addKwv('BINNING', toString(binS), "Bin length", sky_head)    # #
    sky_head <- addKwv('SOURCE', 'process_Beam.R',                           # #
                       "Script used to generate this file", sky_head)        # #
    sky_head <- addKwv('AUTHOR', 'JRS2023', "Script creator", sky_head)      # #
    I_head <- sky_head                                                       # #
    I_head <- modVal('FILETYP2', 'I Flux', "File sub-type", I_head)          # #
                                                                             # #
    # Prepare attributes for output files' header                            # #
    print(paste0("* Preparing attributes for output files' headers for set ",# #
                 "#", n, " of band ", bs, " *"))                             # #
                                                                             # #
    refx <- get_fits_header_num(temp_head, REFX)                             # #
    refy <- get_fits_header_num(temp_head, REFY)                             # #
    refz <- 1                                                                # #
    refz2 <- 1                                                               # #
    valx <- get_fits_header_num(temp_head, VALX)                             # #
    valy <- get_fits_header_num(temp_head, VALY)                             # #
    valz <- NA                                                               # #
    valz2 <- NA                                                              # #
    typex <- "  'PIXEL     '              / Coordinate system of x-axis "    # #
    typey <- "  'PIXEL     '              / Coordinate system of y-axis "    # #
    typez <- "  'DATA & UNC'              / Coordinate system of z-axis "    # #
    typez2 <-"  'HWP ANGLE  '              / Coordinate system of z2-axis"   # #
    ref_coords <- c(get_fits_header_num(temp_head, REFRA),                   # #
                    get_fits_header_num(temp_head, REFDEC))                  # #
    crpix <- c(refx, refy, refz)                                             # #
    crval <- c(valx, valy, valz)                                             # #
    ctype <- c(typex, typey, typez)                                          # #
    crpix2 <- c(refx, refy, refz, refz2)                                     # #
    crval2 <- c(valx, valy, valz, valz2)                                     # #
    ctype2 <- c(typex, typey, typez, typez2)                                 # #
    rm(refx, refy, refz, refz2, valx, valy, valz, valz2, typex, typey,       # #
       typez, typez2)                                                        # #
                                                                             # #
    if(n == 1){                                                              # #
      ref_crpix <- crpix[1:2]                                                # #
    }                                                                        # #
                                                                             # #
    date_name <- get_fits_header_str(temp_head, DATE)                        # #
    acq_date <- strsplit(date_name, "T")[[1]][1]                             # #
    acq_year <- as.POSIXct(acq_date, tz = "UTC")                             # #
    target_name <- get_fits_header_str(temp_head, TARGET)                    # #
    fits_fname <- paste0("FORS2_", target_name, "_", date_name, "_", bs,     # #
                         "_#", n)                                            # #
    fits_name <- paste0(fits_fname, run_date)                                # #
                                                                             # #
    # Calculate scale relationship between pixel unit and fov angular unit   # #
    binx <- get_fits_header_num(temp_head, BINX)                             # #
    biny <- get_fits_header_num(temp_head, BINY)                             # #
    pxscl <- get_fits_header_num(temp_head, PIXSCALE)                        # #
    scl_x <- pxscl * binx / 3600                                             # #
    scl_y <- pxscl * biny / 3600                                             # #
    gain <- get_fits_header_num(temp_head, GAIN)                             # #
    ascl <- pxscl * mean(c(binx, biny))                                      # #
                                                                             # #
    py_run_string("gain_m = r.gain")                                         # #
                                                                             # #
    # Get Coordinates and Diameter of object                                 # #
    tgt_coords <- get_NED_coords(target_name, 6)                             # #
    temp_a <- get_NED_median_diam(target_name)                               # #
    tgt_sma <- round(temp_a / (2 * ascl))                                    # #
    b_over_a <- get_NED_largest_axis_ratio(target_name)                      # #
    tgt_pa <- (get_NED_newest_posang(target_name) - 90) / 180 * pi           # #
    rm(target_name)                                                          # #
                                                                             # #
    # Difference between celestial coordinates of image reference and object # #
    dif_coords <- ref_coords - tgt_coords                                    # #
    cosd <- cos(mean(c(ref_coords[2], tgt_coords[2])) * pi / 180)            # #
                                                                             # #
    # Pixel coordinates of object within image                               # #
    tgt_xy <- round(crpix[1:2] - c(dif_coords[1] * cosd / scl_x,             # #
                                   dif_coords[2] / scl_y))                   # #
                                                                             # #
    # Find Acquisition image that matches in date and offset with            # #
    # polarimetric observations                                              # #
    found_acq <- FALSE                                                       # #
    aN <- 1                                                                  # #
    acq_ind <- NULL                                                          # #
                                                                             # #
    while(!found_acq && aN <= acqN){                                         # #
      t_acq <- acqList[aN]                                                   # #
      t_date <- strsplit(strsplit(t_acq, "FORS2_")[[1]][2], "T")[[1]][1]     # #
                                                                             # #
      if(t_date == acq_date){                                                # #
        acq_ind <- aN                                                        # #
        found_acq <- TRUE                                                    # #
      }                                                                      # #
                                                                             # #
      aN <- aN + 1                                                           # #
    }                                                                        # #
    rm(t_acq, t_date, aN)                                                    # #
                                                                             # #
    ############################################################################
                                                                               #
    ########################## Creating Target Mask ########################## #
    ############################################################################
    print(paste0("* Creating target object mask for set #", n, " of band ",  # #
                 bs, " *"))                                                  # #
                                                                             # #
    # If there's not necessary astrometric information                       # #
    if(is.na(b_over_a) | is.na(tgt_pa)){                                     # #
      # Nor is there an aquisition image on which to do isophote analysis    # #
      if(is.null(acq_ind)){                                                  # #
        maskGal <- drawCircle(img = array(0, dim = dimsObj), x = tgt_xy[1],  # #
                              y = tgt_xy[2], radius = tgt_sma + 1, fill = T, # #
                              col = 1)                                       # #
      }else{# But there is an acquisition file to do isophote analysis       # #
        maskGal <- array(0, dim = dimsObj)                                   # #
                                                                             # #
        acq_head <- readFITS(acqList[acq_ind])$header                        # #
                                                                             # #
        acq_name <- get_fits_header_str(acq_head, TARGET)                    # #
                                                                             # #
        # Get reference coords for acq                                       # #
        acq_refx <- get_fits_header_num(acq_head, REFX)                      # #
        acq_refy <- get_fits_header_num(acq_head, REFY)                      # #
        acq_crpix <- c(acq_refx, acq_refy)                                   # #
        acq_ref_coords <- c(get_fits_header_num(acq_head, REFRA),            # #
                            get_fits_header_num(acq_head, REFDEC))           # #
        rm(acq_refx, acq_refy)                                               # #
                                                                             # #
        acq_dif_coords <- acq_ref_coords - tgt_coords                        # #
        acq_cosd <- cos(mean(c(acq_ref_coords[2], tgt_coords[2])) * pi / 180)# #
        acq_tgt_xy <- round(acq_crpix[1:2] - c(acq_dif_coords[1] *           # #
                                                 acq_cosd / scl_x,           # #
                                               acq_dif_coords[2] / scl_y))   # #
        rm(acq_dif_coords, acq_cosd)                                         # #
                                                                             # #
        # Estimate pixel coordinates difference between acq ref & object     # #
        dif_tgt_xy <- acq_tgt_xy - tgt_xy                                    # #
                                                                             # #
        # Load Acquisition image and target info to python                   # #
        load_cmd <- paste0("acqF = fits.open('", acqList[acq_ind], "')")     # #
        py_run_string(load_cmd)                                              # #
        rm(load_cmd)                                                         # #
                                                                             # #
        py_run_string("img_dat = acqF[0].data[0]")                           # #
        py_run_string("acqF.close()")                                        # #
        py_run_string("X = r.acq_tgt_xy[0]")                                 # #
        py_run_string("Y = r.acq_tgt_xy[1]")                                 # #
        py_run_string("SMA = r.tgt_sma")                                     # #
        rm(acq_tgt_xy)                                                       # #
                                                                             # #
        # Python commands to get shape of target mask                        # #
        get_cmd <- paste0("iso = sexD.get_closest_iso_for_target(img_dat, X",# #
                          ", Y, SMA, 0, 90)")                                # #
        py_run_string(get_cmd)                                               # #
        rm(get_cmd)                                                          # #
                                                                             # #
        # Get isophotal info from Python to R                                # #
        tgt_x <- py$iso$x0 - dif_tgt_xy[1]                                   # #
        tgt_y <- py$iso$y0 - dif_tgt_xy[2]                                   # #
        tgt_eps <- py$iso$eps                                                # #
        tgt_pa <- py$iso$pa                                                  # #
                                                                             # #
        tgt_a <- py$iso$sma                                                  # #
        tgt_b <- tgt_a * sqrt(1 - tgt_eps)                                   # #
        py_run_string("del(X, Y, SMA, iso, img_dat)")                        # #
        rm(tgt_eps)                                                          # #
                                                                             # #
        for(xa in 1:dimsObj[1]){                                             # #
          for(ya in 1:dimsObj[2]){                                           # #
            rx <- ((xa - tgt_x) * cos(tgt_pa) + (ya - tgt_y) *               # #
                     sin(tgt_pa))^2 / tgt_a^2                                # #
            ry <- ((xa - tgt_x) * sin(tgt_pa) - (ya - tgt_y) *               # #
                     cos(tgt_pa))^2 / tgt_b^2                                # #
                                                                             # #
            if(rx + ry <= 1){                                                # #
              maskGal[xa, ya] <- 1                                           # #
            }                                                                # #
          }                                                                  # #
        }                                                                    # #
      }                                                                      # #
    }else{# If there is enough astrometry to draw the ellipse                # #
      maskGal <- array(0, dim = dimsObj)                                     # #
                                                                             # #
      tgt_a <- round((temp_a + sqrt(temp_a)) / (2 * ascl))                   # #
      tgt_b <- b_over_a * tgt_a                                              # #
                                                                             # #
      for(xa in 1:dimsObj[1]){                                               # #
        for(ya in 1:dimsObj[2]){                                             # #
          rx <- ((xa - tgt_xy[1]) * cos(tgt_pa) + (ya - tgt_xy[2]) *         # #
                   sin(tgt_pa))^2 / tgt_a^2                                  # #
          ry <- ((xa - tgt_xy[1]) * sin(tgt_pa) - (ya - tgt_xy[2]) *         # #
                   cos(tgt_pa))^2 / tgt_b^2                                  # #
                                                                             # #
          if(rx + ry <= 1){                                                  # #
            maskGal[xa, ya] <- 1                                             # #
          }                                                                  # #
        }                                                                    # #
      }                                                                      # #
      rm(xa, ya, rx, ry)                                                     # #
    }                                                                        # #
    ############################################################################
                                                                               #
    #------------------------------------------------------------------------# #
                                                                               #
    ########################## Backing Up NA Coords ########################## #
    ############################################################################
    print(paste0("* Removing Infs and saturated pixels from beams of set #", # #
                 n, " of band ", bs, " *"))                                  # #
                                                                             # #
    Os <- remove_infs(Os, NA)                                                # #
    Es <- remove_infs(Es, NA)                                                # #
    Os[,,"Data",] <- replace_qtts(Os[,,"Data",], sat, '>', NA)               # #
    Es[,,"Data",] <- replace_qtts(Es[,,"Data",], sat, '>', NA)               # #
    Os[,,"Unc",] <- bij_match_val(Os[,,"Data",], Os[,,"Unc",], NA, NA)       # #
    Es[,,"Unc",] <- bij_match_val(Es[,,"Data",], Es[,,"Unc",], NA, NA)       # #
                                                                             # #
    out_I_fits <- paste0(obs_Iflux_folder, fits_fname, "-I_obs.fits")        # #
                                                                             # #
    Iobs <- array(NA, dim = singDim, dimnames = singNames)                   # #
                                                                             # #
    if(!file.exists(out_I_fits)){                                            # #
      print(paste0("* Calculating observed total flux for set #", n,         # #
                   " of band ", bs, " *"))                                   # #
                                                                             # #
      Is <- array(NA, dim = multDim, dimnames = multNames)                   # #
                                                                             # #
      Is[,,"Data",] <- Es[,,"Data",] + Os[,,"Data",]                         # #
      Is[,,"Unc",] <- unc_add(Es[,,"Unc",], Os[,,"Unc",])                    # #
      Iobs[,,"Data"] <- apply(Is[,,"Data",], 1:2, median, na.rm = T)         # #
      Iobs[,,"Unc"] <- unc_median(Is[,,"Data",], Is[,,"Unc",], 1:2)          # #
      rm(Is)                                                                 # #
                                                                             # #
      writeFITSim(Iobs, file = out_I_fits, crvaln = crval, crpixn = crpix,   # #
                  ctypen = ctype, header = I_head)                           # #
    }else{                                                                   # #
      print(paste0("* Loading observed total flux for set #", n, " of band ",# #
                   bs, " *"))                                                # #
                                                                             # #
      Iobs[,,] <- readFITS(out_I_fits)$imDat                                 # #
    }                                                                        # #
    rm(out_I_fits)                                                           # #
                                                                             # #
    print(paste0("* Backing up pixel coordinates for instrumentaly masked p",# #
                 "ixels for set #", n, " of band ", bs, " *"))               # #
                                                                             # #
    NAs_ind <- which(is.na(Iobs[,,"Data"]), arr.ind = TRUE)                  # #
    E_na_inds <- which(is.na(Es[,,"Data",]), arr.ind = TRUE)                 # #
    O_na_inds <- which(is.na(Os[,,"Data",]), arr.ind = TRUE)                 # #
                                                                             # #
    Es_pix_mask <- array(FALSE, dim = c(dimsObj, length(ang_str)))           # #
    Os_pix_mask <- array(FALSE, dim = c(dimsObj, length(ang_str)))           # #
                                                                             # #
    Es_pix_mask[E_na_inds] <- TRUE                                           # #
    Os_pix_mask[O_na_inds] <- TRUE                                           # #
    ############################################################################
                                                                               #
    e_ext_mask_out_fits <- paste0(trgt_mask_folder, fits_fname,                #
                                  "-src_mask_4_bkg_Eb.fits")                   #
    o_ext_mask_out_fits <- paste0(trgt_mask_folder, fits_fname,                #
                                  "-src_mask_4_bkg_Ob.fits")                   #
    e_pnt_mask_out_fits <- paste0(trgt_mask_folder, fits_fname,                #
                                  "-src_mask_4_bkg_E_MW.fits")                 #
    o_pnt_mask_out_fits <- paste0(trgt_mask_folder, fits_fname,                #
                                  "-src_mask_4_bkg_O_MW.fits")                 #
    e_sky_out_fits <- paste0(trgt_sky_folder, fits_fname, "-Es_Sky.fits")      #
    o_sky_out_fits <- paste0(trgt_sky_folder, fits_fname, "-Os_Sky.fits")      #
    e_mw_sky_out_fits <- paste0(mw_sky_folder, fits_fname, "-Es_MW_Sky.fits")  #
    o_mw_sky_out_fits <- paste0(mw_sky_folder, fits_fname, "-Os_MW_Sky.fits")  #
                                                                               #
    exist_flag <- prod(file.exists(c(e_ext_mask_out_fits, e_sky_out_fits,      #
                                     o_ext_mask_out_fits, o_sky_out_fits,      #
                                     e_pnt_mask_out_fits, e_mw_sky_out_fits,   #
                                     o_pnt_mask_out_fits, o_mw_sky_out_fits))) #
                                                                               #
    Es_sky <- array(NA, dim = multDim, dimnames = multNames)                   #
    Os_sky <- array(NA, dim = multDim, dimnames = multNames)                   #
    Es_MW_sky <- array(NA, dim = multDim, dimnames = multNames)                #
    Os_MW_sky <- array(NA, dim = multDim, dimnames = multNames)                #
                                                                               #
    if(!exist_flag){                                                           #
      ########################## Creating Source Masks ####################### #
      ##########################################################################
      print(paste0("* Creating source masks for beams for set #", n,         # #
                   " of band ", bs, " *"))                                   # #
                                                                             # #
      Es_pnt_src_mask <- array(FALSE, dim = c(dimsObj, length(ang_str)))     # #
      Os_pnt_src_mask <- array(FALSE, dim = c(dimsObj, length(ang_str)))     # #
      Es_ext_src_mask <- array(FALSE, dim = c(dimsObj, length(ang_str)))     # #
      Os_ext_src_mask <- array(FALSE, dim = c(dimsObj, length(ang_str)))     # #
                                                                             # #
      E_pnt_mask_cmd <- paste0("E_pnt_mask = make_source_mask(Ebeam, nsigma",# #
                               " =", pnt_src_t, ", npixels = 5, dilate_size",# #
                               " = 3)")                                      # #
      O_pnt_mask_cmd <- paste0("O_pnt_mask = make_source_mask(Obeam, nsigma",# #
                               " =", pnt_src_t, ", npixels = 5, dilate_size",# #
                               " = 3)")                                      # #
                                                                             # #
      for(ang in angs){                                                      # #
        as <- ang_str[ang]                                                   # #
                                                                             # #
        print(paste0("--> Loading, to Python Env., E and O beam data of HWP",# #
                     " angle = ", as, "º..."))                               # #
                                                                             # #
        temp <- Es[,, "Data", ang]                                           # #
        py_run_string("Ebeam = r.temp")                                      # #
        temp <- Os[,, "Data", ang]                                           # #
        py_run_string("Obeam = r.temp")                                      # #
        rm(temp)                                                             # #
                                                                             # #
        print(paste0("--> Creating preliminary point source mask for E beam",# #
                     " data of HWP angle = ", as, "º..."))                   # #
                                                                             # #
        py_run_string(E_pnt_mask_cmd)                                        # #
                                                                             # #
        print(paste0("--> Creating preliminary point source mask for O beam",# #
                     " data of HWP angle = ", as, "º..."))                   # #
                                                                             # #
        py_run_string(O_pnt_mask_cmd)                                        # #
                                                                             # #
        py_run_string("del(Ebeam, Obeam)")                                   # #
                                                                             # #
        print(paste0("--> Saving masks created for beam data of HWP angle =",# #
                     " ", as,  "º to R Env. objects..."))                    # #
                                                                             # #
        Es_pnt_src_mask[,,ang] <- py$E_pnt_mask                              # #
        Os_pnt_src_mask[,,ang] <- py$O_pnt_mask                              # #
      }                                                                      # #
      py_run_string("del(E_pnt_mask, O_pnt_mask)")                           # #
      rm(E_pnt_mask_cmd, O_pnt_mask_cmc)                                     # #
                                                                             # #
      Es_pnt_src_mask[E_na_inds] <- TRUE                                     # #
      Os_pnt_src_mask[O_na_inds] <- TRUE                                     # #
      rm(E_na_inds, O_na_inds)                                               # #
                                                                             # #
      print(paste0("* Masking galaxy in E_beams and O_beams for set #", n,   # #
                   " of band ", bs, " *"))                                   # #
      for(ang in angs){                                                      # #
        Es_ext_src_mask[,,ang] <- Es_pnt_src_mask[,,ang] | maskGal           # #
        Os_ext_src_mask[,,ang] <- Os_pnt_src_mask[,,ang] | maskGal           # #
      }                                                                      # #
      ##########################################################################
                                                                               #
      #----------------------------------------------------------------------# #
                                                                               #
      ########################## Estimating Sky Beams ######################## #
      ##########################################################################
      nl_E <- dim(E_y_lims)[1]                                               # #
      nl_O <- dim(O_y_lims)[1]                                               # #
                                                                             # #
      E_ext_bkg_cmd <- paste0("E_ext_bkg = sep.Background(Ebeam, E_ext_mask",# #
                              ", maskthresh = 0, bw = 64, bh = 64, fw = ",   # #
                              psf_est, ", fh = ", psf_est,")")               # #
      E_pnt_bkg_cmd <- paste0("E_pnt_bkg = sep.Background(Ebeam, E_pnt_mask",# #
                              ", maskthresh = 0, bw = 64, bh = 64, fw = ",   # #
                              psf_est, ", fh = ", psf_est,")")               # #
      O_ext_bkg_cmd <- paste0("O_ext_bkg = sep.Background(Obeam, O_ext_mask",# #
                              ", maskthresh = 0, bw = 64, bh = 64, fw = ",   # #
                              psf_est, ", fh = ", psf_est,")")               # #
      O_pnt_bkg_cmd <- paste0("O_pnt_bkg = sep.Background(Obeam, O_pnt_mask",# #
                              ", maskthresh = 0, bw = 64, bh = 64, fw = ",   # #
                              psf_est, ", fh = ", psf_est,")")               # #
                                                                             # #
      print(paste0("* Calculating backgrounds for E and O beams for set #",  # #
                   n, " of band ", bs, " *"))                                # #
                                                                             # #
      for(ang in angs){                                                      # #
        as <- ang_str[ang]                                                   # #
                                                                             # #
        xs_E <- E_x_lims[1, ang]:E_x_lims[2, ang]                            # #
        xs_O <- O_x_lims[1, ang]:O_x_lims[2, ang]                            # #
                                                                             # #
        for(nl in seq(1, nl_E, 2)){                                          # #
          sn <- (nl + 1) / 2                                                 # #
          ys_E <- E_y_lims[nl, ang]:E_y_lims[nl + 1, ang]                    # #
                                                                             # #
          print(paste0("--> Loading, to Python Env., E beam data of HWP ang",# #
                       "le = ", as, "º: stripe ", sn, "..."))                # #
          temp <- Es[xs_E, ys_E, "Data", ang]                                # #
          py_run_string("Ebeam = r.temp.copy(order='C')")                    # #
                                                                             # #
          print(paste0("--> Loading, to Python Env., E beam extended source",# #
                       " mask of HWP angle = ", as, "º: stripe ", sn, "..."))# #
          temp <- Es_ext_src_mask[xs_E, ys_E, ang]                           # #
          py_run_string("E_ext_mask = r.temp.copy(order='C')")               # #
                                                                             # #
          print(paste0("--> Loading, to Python Env., E beam point source ma",# #
                       "sk of HWP angle = ", as, "º: stripe ", sn, "..."))   # #
          temp <- Es_pnt_src_mask[xs_E, ys_E, ang]                           # #
          py_run_string("E_pnt_mask = r.temp.copy(order='C')")               # #
          rm(temp)                                                           # #
                                                                             # #
          print(paste0("--> Estimating SKY level for E beam data of HWP ang",# #
                       "le = ", as, "º: stripe ", sn, "..."))                # #
          py_run_string(E_ext_bkg_cmd)                                       # #
                                                                             # #
          print(paste0("--> Estimating MW SKY level for E beam data of HWP ",# #
                       "angle = ", as, "º: stripe ", sn, "..."))             # #
          py_run_string(E_pnt_bkg_cmd)                                       # #
                                                                             # #
          py_run_string("del(Ebeam, E_ext_mask, E_pnt_mask)")                # #
                                                                             # #
          py_run_string("E_ext_sky = E_ext_bkg.back()")                      # #
          py_run_string("E_ext_sky_U = E_ext_bkg.rms()")                     # #
          py_run_string("E_pnt_sky = E_pnt_bkg.back()")                      # #
          py_run_string("E_pnt_sky_U = E_pnt_bkg.rms()")                     # #
          py_run_string("del(E_ext_bkg, E_pnt_bkg)")                         # #
                                                                             # #
          print(paste0("--> Saving SKY level maps for E beam data of HWP an",# #
                       "gle = ", as,  "º to R Env. objects: stripe ", sn,    # #
                       "..."))                                               # #
                                                                             # #
          # Place background stripes on maps                                 # #
          Es_sky[xs_E, ys_E, "Data", ang] <- py$E_ext_sky                    # #
          Es_sky[xs_E, ys_E, "Unc", ang] <- py$E_ext_sky_U                   # #
          Es_MW_sky[xs_E, ys_E, "Data", ang] <- py$E_pnt_sky                 # #
          Es_MW_sky[xs_E, ys_E, "Unc", ang] <- py$E_pnt_sky_U                # #
                                                                             # #
          py_run_string("del(E_ext_sky, E_ext_sky_U,E_pnt_sky,E_pnt_sky_U)") # #
        }                                                                    # #
                                                                             # #
        for(nl in seq(1, nl_O, 2)){                                          # #
          sn <- (nl + 1) / 2                                                 # #
          ys_O <- O_y_lims[nl, ang]:O_y_lims[nl + 1, ang]                    # #
                                                                             # #
          print(paste0("--> Loading, to Python Env., O beam data of HWP ang",# #
                       "le = ", as, "º: stripe ", sn, "..."))                # #
          temp <- Os[xs_O, ys_O, "Data", ang]                                # #
          py_run_string("Obeam = r.temp.copy(order='C')")                    # #
                                                                             # #
          print(paste0("--> Loading, to Python Env., O beam extended source",# #
                       " mask of HWP angle = ", as, "º: stripe ", sn, "..."))# #
          temp <- Os_ext_src_mask[xs_O, ys_O, ang]                           # #
          py_run_string("O_ext_mask = r.temp.copy(order='C')")               # #
                                                                             # #
          print(paste0("--> Loading, to Python Env., O beam point source ma",# #
                       "sk of HWP angle = ", as, "º: stripe ", sn, "..."))   # #
          temp <- Os_pnt_src_mask[xs_O, ys_O, ang]                           # #
          py_run_string("O_pnt_mask = r.temp.copy(order='C')")               # #
          rm(temp)                                                           # #
                                                                             # #
          print(paste0("--> Estimating SKY level for O beam data of HWP ang",# #
                       "le = ", as, "º: stripe ", sn, "..."))                # #
          py_run_string(O_ext_bkg_cmd)                                       # #
                                                                             # #
          print(paste0("--> Estimating MW SKY level for O beam data of HWP ",# #
                       "angle = ", as, "º: stripe ", sn, "..."))             # #
          py_run_string(O_pnt_bkg_cmd)                                       # #
                                                                             # #
          py_run_string("del(Obeam, O_ext_mask, O_pnt_mask)")                # #
                                                                             # #
          py_run_string("O_ext_sky = O_ext_bkg.back()")                      # #
          py_run_string("O_ext_sky_U = O_ext_bkg.rms()")                     # #
          py_run_string("O_pnt_sky = O_pnt_bkg.back()")                      # #
          py_run_string("O_pnt_sky_U = O_pnt_bkg.rms()")                     # #
          py_run_string("del(O_ext_bkg, O_pnt_bkg)")                         # #
                                                                             # #
          print(paste0("--> Saving SKY level maps for O beam data of HWP an",# #
                       "gle = ", as,  "º to R Env. objects: stripe ", sn,    # #
                       "..."))                                               # #
                                                                             # #
          # Place background stripes on maps                                 # #
          Os_sky[xs_O, ys_O, "Data", ang] <- py$O_ext_sky                    # #
          Os_sky[xs_O, ys_O, "Unc", ang] <- py$O_ext_sky_U                   # #
          Os_MW_sky[xs_O, ys_O, "Data", ang] <- py$O_pnt_sky                 # #
          Os_MW_sky[xs_O, ys_O, "Unc", ang] <- py$O_pnt_sky_U                # #
                                                                             # #
          py_run_string("del(O_ext_sky, O_ext_sky_U,O_pnt_sky,O_pnt_sky_U)") # #
        }                                                                    # #
      }                                                                      # #
                                                                             # #
      print(paste0("* Saving all sources masks to fits format for set #", n, # #
                   " of band ", bs, " *"))                                   # #
                                                                             # #
      # Save source mask used to determine background                        # #
      writeFITSim(Es_ext_src_mask, file = e_ext_mask_out_fits)               # #
      writeFITSim(Os_ext_src_mask, file = o_ext_mask_out_fits)               # #
      writeFITSim(Es_pnt_src_mask, file = e_pnt_mask_out_fits)               # #
      writeFITSim(Os_pnt_src_mask, file = o_pnt_mask_out_fits)               # #
      rm(Es_ext_src_mask, Os_ext_src_mask, Es_pnt_src_mask, Os_pnt_src_mask) # #
                                                                             # #
      print(paste0("* Saving background maps inferred from all sources mask",# #
                   "s, to fits format for set #", n, " of band ", bs, " *")) # #
                                                                             # #
      # Save backgrounds for future checking                                 # #
      writeFITSim(Es_MW_sky, file = e_mw_sky_out_fits, crvaln = crval2,      # #
                  crpixn = crpix2, ctypen = ctype2, header = skyE_head)      # #
      writeFITSim(Os_MW_sky, file = o_mw_sky_out_fits, crvaln = crval2,      # #
                  crpixn = crpix2, ctypen = ctype2, header = skyO_head)      # #
      writeFITSim(Es_sky, file = e_sky_out_fits, crvaln = crval2,            # #
                  crpixn = crpix2, ctypen = ctype2, header = skyE_head)      # #
      writeFITSim(Os_sky, file = o_sky_out_fits, crvaln = crval2,            # #
                  crpixn = crpix2, ctypen = ctype2, header = skyO_head)      # #
      ##########################################################################
    }else{                                                                     #
      ########################### Loading Sky Beams ########################## #
      ##########################################################################
      print(paste0("* Loading background maps for set #", n, " of band ", bs,# #
                   " *"))                                                    # #
                                                                             # #
      Es_MW_sky[,,,] <- readFITS(e_mw_sky_out_fits)$imDat                    # #
      Os_MW_sky[,,,] <- readFITS(o_mw_sky_out_fits)$imDat                    # #
      Es_sky[,,,] <- readFITS(e_sky_out_fits)$imDat                          # #
      Os_sky[,,,] <- readFITS(o_sky_out_fits)$imDat                          # #
      ##########################################################################
    }                                                                          #
    rm(e_mw_sky_out_fits, o_mw_sky_out_fits, e_sky_out_fits, o_sky_out_fits,   #
       e_ext_mask_out_fits, o_ext_mask_out_fits, e_pnt_mask_out_fits,          #
       o_pnt_mask_out_fits, exist_flag)                                        #
                                                                               #
    I_sky_out_fits <- paste0(sky_Iflux_folder, fits_fname, "-I_Sky.fits")      #
                                                                               #
    Isky <- array(NA, dim = singDim, dimnames = singNames)                     #
                                                                               #
    if(!file.exists(I_sky_out_fits)){                                          #
      ########################## Estimating Sky Flux ######################### #
      ##########################################################################
      print(paste0("* Estimating total background flux for set #", n,        # #
                   " of band ", bs, " *"))                                   # #
                                                                             # #
      Is_sky <- Es_sky[,,"Data",] + Os_sky[,,"Data",]                        # #
      Is_sky_unc <- unc_add(Es_sky[,,"Unc",], Os_sky[,,"Unc",])              # #
      Isky[,,"Data"] <- apply(Is_sky, 1:2, median, na.rm = T)                # #
      Isky[,,"Unc"] <- unc_median(Is_sky, Is_sky_unc, 1:2)                   # #
      rm(Is_sky, Is_sky_unc)                                                 # #
                                                                             # #
      writeFITSim(Isky, file = I_sky_out_fits, crvaln = crval,               # #
                  crpixn = crpix, ctypen = ctype, header = I_head)           # #
      ##########################################################################
    }else{                                                                     #
      ############################ Loading Sky Flux ########################## #
      ##########################################################################
      print(paste0("* Loading background flux map for set #", n, " of band ",# #
                   bs, " *"))                                                # #
                                                                             # #
      Isky[,,] <- readFITS(I_sky_out_fits)$imDat                             # #
      ##########################################################################
    }                                                                          #
    rm(I_sky_out_fits)                                                         #
                                                                               #
    mwqui_out_csv <- paste0(mwPolCSV_path, fits_fname, "-MW_QU_instC.csv")     #
    mwqupx_out_csv <- paste0(MWstars_folder, fits_fname, "-MW_STARS_pol.csv")  #
                                                                               #
    exist_flag <- prod(file.exists(c(mwqui_out_csv, mwqupx_out_csv)))          #
                                                                               #
    munc_par <- c("Data", "Unc")                                               #
    beam_par <- c('E','O')                                                     #
    beam_ang_par <- c("E0", "E22", "E45", "E68", "O0", "O22", "O45", "O68")    #
    pos_par <- c("x", "y")                                                     #
    srcs_par <- c('x','y', 'Q','uQ', 'U', 'uU', 'P','uP', 'X','uX')            #
    stokes_par <-  c('Q', 'U')                                                 #
    polV_par <- c('P', 'X')                                                    #
                                                                               #
    if(!exist_flag){                                                           #
      rej_csv_4 <- paste0(mwSrcREJ_path, fits_name, "-rej4_src_lst.csv")       #
      acp_csv_4 <- paste0(mwSrcACP_path, fits_name, "-acp4_src_lst.csv")       #
      rej_csv_f <- paste0(mwSrcREJ_path, fits_name, "-rejFin_src_lst.csv")     #
      acp_csv_f <- paste0(mwSrcACP_path, fits_name, "-acpFin_src_lst.csv")     #
                                                                               #
      mwflux_out_csv <- paste0(mwSrcFLX_path, fits_fname, "-Fluxes.csv")       #
                                                                               #
      mwqur_out_csv <- paste0(mwPolCSV_path, fits_fname, "-MW_QU_raw.csv")     #
      mwpxr_out_csv <- paste0(mwPolCSV_path, fits_fname, "-MW_PX_raw.csv")     #
      mwpxi_out_csv <- paste0(mwPolCSV_path, fits_fname, "-MW_PX_instC.csv")   #
      mwqucc_out_csv <- paste0(mwPolCSV_path, fits_fname, "-MW_QU_CC.csv")     #
      mwpxcc_out_csv <- paste0(mwPolCSV_path, fits_fname, "-MW_PX_CC.csv")     #
      mwpxdeb_out_csv <- paste0(mwPolCSV_path, fits_fname, "-MW_PX_allC.csv")  #
                                                                               #
      src_polr_out_csv <- paste0(mwSrcPOL_path, fits_fname,"-srcs_raw_pol.csv")#
      src_poli_out_csv <-paste0(mwSrcPOL_path,fits_fname,"-srcs_instC_pol.csv")#
      src_polc_out_csv <-paste0(mwSrcPOL_path, fits_fname,"-srcs_allC_pol.csv")#
                                                                               #
      mwqur_out_pdf <- paste0(mwPolPDF_path, fits_fname, "-QvsU_mw_raw.pdf")   #
      mwqui_out_pdf <- paste0(mwPolPDF_path, fits_fname, "-QvsU_mw_instC.pdf") #
      mwqucc_out_pdf <- paste0(mwPolPDF_path, fits_fname, "-QvsU_mw_allC.pdf") #
                                                                               #
      #################### Performing Photometry on Beams #################### #
      ##########################################################################
                                                                             # #
      Emw_out_fits <- paste0(mw_noSky_folder, fits_fname, "-E_noMWSky.fits") # #
      Omw_out_fits <- paste0(mw_noSky_folder, fits_fname, "-O_noMWSky.fits") # #
                                                                             # #
      mw_E <- array(NA, dim = multDim, dimnames = multNames)                 # #
      mw_O <-array(NA, dim = multDim, dimnames = multNames)                  # #
                                                                             # #
      if(!prod(file.exists(c(Emw_out_fits, Omw_out_fits)))){                 # #
        ########## Subtracting Foregrounds Background from Beams ########### # #
        ########################################################################
        print(paste0("* Removing point source background from beams to mo",# # #
                     "re accurately extract Milky Way sources for set #",  # # #
                     n, " of band ", bs, " *"))                            # # #
                                                                           # # #
        mw_E[,,"Data",] <- replace_qtts(Es[,,"Data",]-Es_MW_sky[,,"Data",],# # #
                                        0, '<', 0)                         # # #
        mw_O[,,"Data",] <- replace_qtts(Os[,,"Data",]-Os_MW_sky[,,"Data",],# # #
                                        0, '<', 0)                         # # #
                                                                           # # #
        temp <- unc_add(Es[,,"Unc",], Es_MW_sky[,,"Unc",])                 # # #
        mw_E[,,"Unc",] <- bij_match_val(mw_E[,,"Data",], temp, 0, NA)      # # #
        temp <- unc_add(Os[,,"Unc",], Os_MW_sky[,,"Unc",])                 # # #
        mw_O[,,"Unc",] <- bij_match_val(mw_O[,,"Data",], temp, 0, NA)      # # #
        rm(Es_MW_sky, Os_MW_sky, temp)                                     # # #
                                                                           # # #
        writeFITSim(mw_O, file = Omw_out_fits, crvaln = crval2,            # # #
                    crpixn = crpix2, ctypen = ctype2)                      # # #
        writeFITSim(mw_E, file = Emw_out_fits, crvaln = crval2,            # # #
                    crpixn = crpix2, ctypen = ctype2)                      # # #
        ########################################################################
      }else{                                                                 # #
        ######### Loading Foregrounds Background Subtracted Beams ########## # #
        ########################################################################
        print(paste0("* Loading point source background subtracted beams ",# # #
                     "to more accurately extract Milky Way sources for se",# # #
                     "t #", n, " of band ", bs, " *"))                     # # #
                                                                           # # #
        mw_O[,,,] <- readFITS(Omw_out_fits)$imDat                          # # #
        mw_E[,,,] <- readFITS(Emw_out_fits)$imDat                          # # #
        ########################################################################
      }                                                                      # #
      rm(Emw_out_fits, Omw_out_fits)                                         # #
                                                                             # #
      print(paste0("* Listing candidate Milky Way sources for each HWP angl",# #
                   "e for set #", n, " of band ", bs, " *"))                 # #
                                                                             # #
      cdd_srcs_0 <- NULL                                                     # #
      cdd_srcs_0$E <- NULL                                                   # #
      cdd_srcs_0$O <- NULL                                                   # #
      cdd_srcs_22 <- NULL                                                    # #
      cdd_srcs_22$E <- NULL                                                  # #
      cdd_srcs_22$O <- NULL                                                  # #
      cdd_srcs_45 <- NULL                                                    # #
      cdd_srcs_45$E <- NULL                                                  # #
      cdd_srcs_45$O <- NULL                                                  # #
      cdd_srcs_68 <- NULL                                                    # #
      cdd_srcs_68$E <- NULL                                                  # #
      cdd_srcs_68$O <- NULL                                                  # #
                                                                             # #
      src_map_ang <- array(0, dim = c(dim(mw_E[,,1,])))                      # #
                                                                             # #
      for(ang in angs){                                                      # #
        as <- ang_str[ang]                                                   # #
                                                                             # #
        acp_csv_1 <- paste0(mwSrcACP_path, fits_name, "-acp1_src_lst-HWP=",  # #
                            as, ".csv")                                      # #
        rej_csv_1 <- paste0(mwSrcREJ_path, fits_name, "-rej1_src_lst-HWP=",  # #
                            as, ".csv")                                      # #
        acp_csv_1c <- paste0(mwSrcACP_path, fits_name, "-acp1_src_lst-HWP=", # #
                            as, "_comp-Beam.csv")                            # #
        rej_csv_1c <- paste0(mwSrcREJ_path, fits_name, "-rej1_src_lst-HWP=", # #
                            as, "_comp-Beam.csv")                            # #
        rej_csv_2 <- paste0(mwSrcREJ_path, fits_name, "-rej2_src_lst-HWP=",  # #
                            as, ".csv")                                      # #
        acp_csv_2 <- paste0(mwSrcACP_path, fits_name, "-acp2_src_lst-HWP=",  # #
                            as, ".csv")                                      # #
        rej_csv_2c <- paste0(mwSrcREJ_path, fits_name, "-rej2_src_lst-HWP=", # #
                             as, "_comp-Beam.csv")                           # #
        acp_csv_2c <- paste0(mwSrcACP_path, fits_name, "-acp2_src_lst-HWP=", # #
                             as, "_comp-Beam.csv")                           # #
        rej_csv_3 <- paste0(mwSrcREJ_path, fits_name, "-rej3_src_lst-HWP=",  # #
                            as, ".csv")                                      # #
        acp_csv_3 <- paste0(mwSrcACP_path, fits_name, "-acp3_src_lst-HWP=",  # #
                            as, ".csv")                                      # #
                                                                             # #
        print(paste0("--> Loading, to Python Env., beams data, with extende",# #
                     "d source background level removed, and uncertainties ",# #
                     "of HWP angle = ", as, "º..."))                         # #
                                                                             # #
        temp <- mw_E[,, "Data", ang]                                         # #
        py_run_string("E_temp = r.temp.copy(order='C')")                     # #
        temp <- mw_E[,, "Unc", ang]                                          # #
        py_run_string("E_temp_err = r.temp.copy(order='C')")                 # #
                                                                             # #
        temp <- mw_O[,, "Data", ang]                                         # #
        py_run_string("O_temp = r.temp.copy(order='C')")                     # #
        temp <- mw_O[,, "Unc", ang]                                          # #
        py_run_string("O_temp_err = r.temp.copy(order='C')")                 # #
                                                                             # #
        print(paste0("--> Loading, to Python Env., beam NA pixel mask of HW",# #
                     "P angle = ", as, "º..."))                              # #
                                                                             # #
        temp <- Es_pix_mask[,,ang]                                           # #
        temp[which(mw_E <= 0, arr.ind = T)] <- TRUE                          # #
        py_run_string("E_mask_temp = r.temp.copy(order='C')")                # #
                                                                             # #
        temp <- Os_pix_mask[,,ang]                                           # #
        temp[which(mw_O <= 0, arr.ind = T)] <- TRUE                          # #
        py_run_string("O_mask_temp = r.temp.copy(order='C')")                # #
        rm(temp)                                                             # #
                                                                             # #
        print(paste0("--> Loading O beam to get initial list of candidate M",# #
                     "ilky Way sources for HWP angle ", as, "º..."))         # #
                                                                             # #
        temp_data <- mw_O[,, "Data", ang]                                    # #
                                                                             # #
        # List sources                                                       # #
        print(paste0("--> Creating list and segmentation map of potential M",# #
                     "W sources in O beam for HWP angle ", as, "º..."))      # #
                                                                             # #
        ini_sep_cmd <- paste0("sources_py,seg_map = sep.extract(O_temp, ",   # #
                              src_sig, ", err=O_temp_err, mask=O_mask_temp,",# #
                              " segmentation_map=True)")                     # #
        py_run_string(ini_sep_cmd)                                           # #
        rm(ini_sep_cmd)                                                      # #
                                                                             # #
        print(paste0("--> Saving segmentation map of potential MW sources i",# #
                     "n O beam for HWP angle = ", as, "º to R Env. object..",# #
                     "."))                                                   # #
                                                                             # #
        if(ang == 1){                                                        # #
          tmp_seg_map <- py$seg_map                                          # #
          seg_name <- paste0(seg_path, "Ini_Seg_Map-band=", bs, "_set#=", n, # #
                             run_date, ".fits")                              # #
          writeFITSim(tmp_seg_map, seg_name)                                 # #
          rm(seg_name)                                                       # #
        }                                                                    # #
                                                                             # #
        # Save preliminary source list to a csv file to be read later        # #
        print(paste0("--> Saving list of potential MW sources in O beam for",# #
                     " HWP angle = ", as,  "º to CSV file..."))              # #
                                                                             # #
        ### Reject sources based on area and eccentricity values in O beam   # #
        lst_src_cmd <- paste0("rej_obj = sexD.list_sources_info(O_temp, sou",# #
                              "rces_py,'", acp_csv_1, "', '", rej_csv_1,"',",# #
                              max_ecc, ",", min_A, ",", max_A, ", O_temp_er",# #
                              "r)")                                          # #
                                                                             # #
        py_run_string(lst_src_cmd)                                           # #
        rm(lst_src_cmd)                                                      # #
                                                                             # #
        # Read csv file and filtering source list                            # #
        # List psf parameters of remaining sources                           # #
        print(paste0("--> Reading CSV file with list of potential MW source",# #
                     "s in O beam for HWP angle = ", as, "º to R Env. objec",# #
                     "t..."))                                                # #
                                                                             # #
        srcs_MW <- read.table(acp_csv_1, header = TRUE, sep = ";", dec = ".")# #
        rm(acp_csv_1)                                                        # #
                                                                             # #
        # Exclude rejected objects from segmentation map                     # #
        if(ang == 1){                                                        # #
          print(paste0("--> Excluding rejected objects from segmentation ma",# #
                       "p for O beam and HWP angle = ", as, "º..."))         # #
                                                                             # #
          rej_arr <- array(0, dim = length(py$rej_obj))                      # #
          for(i in 1:dim(rej_arr)){                                          # #
            rej_arr[i] <- py$rej_obj[[i]][1]                                 # #
          }                                                                  # #
          for(i in rej_arr){                                                 # #
            tmp_seg_map[which(tmp_seg_map == i, arr.ind = TRUE)] <- 0        # #
          }                                                                  # #
          rm(rej_arr)                                                        # #
        }                                                                    # #
                                                                             # #
        if(ang == 1){                                                        # #
          seg_name <- paste0(seg_path, "acp1_Seg_Map-band=", bs, "_set#=", n,# #
                             run_date, ".fits")                              # #
          writeFITSim(tmp_seg_map, seg_name)                                 # #
          rm(seg_name)                                                       # #
        }                                                                    # #
                                                                             # #
        print(paste0("--> Estimating some extra parameters for listed poten",# #
                     "tial MW sources in O beam for HWP angle = ",as,"º..."))# #
                                                                             # #
        srcs_MW$r <- (srcs_MW$a + srcs_MW$b) / 2                             # #
        srcs_MW$r_max <- ceiling(sqrt(srcs_MW$Area / pi))                    # #
                                                                             # #
        # Some variables to be used                                          # #
        ini_N <- dim(srcs_MW)[1]                                             # #
        ind_to_rm <- NULL                                                    # #
        obj_to_rm <- NULL                                                    # #
        psf_coefs_ini <- NULL                                                # #
                                                                             # #
        # Check pixels within listed source's area for bad pixels            # #
        print(paste0("--> Checking each potential source in O beam for HWP ",# #
                     "angle = ", as, "º, for bad pixels, flagging those whi",# #
                     "ch have them and saving psf parameters of those which",# #
                     " don't..."))                                           # #
                                                                             # #
        for(s in 1:ini_N){                                                   # #
          cent_pix <- round(c(srcs_MW$x_barycenter[s],                       # #
                              srcs_MW$y_barycenter[s]))                      # #
          dist <- ceiling(srcs_MW$r_max[s])                                  # #
          xs <- (cent_pix[1] - dist):(cent_pix[1] + dist)                    # #
          ys <- (cent_pix[2] - dist):(cent_pix[2] + dist)                    # #
                                                                             # #
          na_pix_flag <- length(which(is.na(temp_data[xs, ys]), arr.ind = T))# #
                                                                             # #
          if(na_pix_flag){                                                   # #
            ind_to_rm <-append(ind_to_rm, s)                                 # #
            obj_to_rm <- append(obj_to_rm, srcs_MW$Seg_n[s])                 # #
            print(paste0("----> Source #", s, " has been flagged..."))       # #
          }                                                                  # #
          else{                                                              # #
            psf_coefs_ini <- rbind(psf_coefs_ini,                            # #
                                   srcs_MW[s, c("Obj_n", "Seg_n")])          # #
          }                                                                  # #
        }                                                                    # #
        rm(ini_N, xs, ys, dist, cent_pix, na_pix_flag)                       # #
                                                                             # #
        ### Remove sources flagged for including NA pixels in O beam         # #
        print(paste0("--> Removing flagged sources in O beam for HWP angle ",# #
                     "= ", as, "º, from list of potential MW sources and fr",# #
                     "om the segmentation map..."))                          # #
                                                                             # #
        rejs_MW <- read.table(rej_csv_1, header = TRUE, sep = ";", dec = ".")# #
        rm(rej_csv_1)                                                        # #
                                                                             # #
        if(!is.null(ind_to_rm)){                                             # #
          rejs_MW <- rbind(rejs_MW, srcs_MW[ind_to_rm, 1:dim(rejs_MW)[2]])   # #
          srcs_MW <- srcs_MW[-ind_to_rm,]                                    # #
                                                                             # #
          if(ang == 1){                                                      # #
            for(o in obj_to_rm){                                             # #
              tmp_seg_map[which(tmp_seg_map == o, arr.ind = TRUE)] <- 0      # #
            }                                                                # #
          }                                                                  # #
        }                                                                    # #
                                                                             # #
        # Saving register of rejected and accepted sources                   # #
        write.table(rejs_MW, rej_csv_2, row.names = F, sep = ";", dec = ".") # #
        write.table(srcs_MW, acp_csv_2, row.names = F, sep = ";", dec = ".") # #
        rm(obj_to_rm, ind_to_rm)                                             # #
                                                                             # #
        if(ang == 1){                                                        # #
          seg_name <- paste0(seg_path, "acp2_Seg_Map-band=", bs, "_set#=", n,# #
                             run_date, ".fits")                              # #
          writeFITSim(tmp_seg_map, seg_name)                                 # #
          rm(seg_name)                                                       # #
        }                                                                    # #
                                                                             # #
        print(paste0("--> Estimating PSF parameters for list of potential M",# #
                     "W sources in O beam for HWP angle = ", as, "º..."))    # #
        srcs_MW$sig <- srcs_MW$r                                             # #
        srcs_MW$fwhm <- srcs_MW$sig * 2 * sqrt(2 * log(2))                   # #
                                                                             # #
        # This will go into python environment which indexes fits images     # #
        # differently than R hence the switch between x and y                # #
        psf_coefs_ini$y <- srcs_MW$x_barycenter - 1                          # #
        psf_coefs_ini$x <- srcs_MW$y_barycenter - 1                          # #
        psf_coefs_ini$r_fwhm <- srcs_MW$fwhm / 2                             # #
        psf_coefs_ini <- psf_coefs_ini[, -1]                                 # #
        rm(srcs_MW)                                                          # #
                                                                             # #
        # Determine radii that yield optimal SNR                             # #
        if(ang == 1){                                                        # #
          print(paste0("--> Loading, to Python Env., sources PSF parameters",# #
                       " in O beam for HWP angle = ", as, "º..."))           # #
          py_run_string("stars = r.psf_coefs_ini")                           # #
                                                                             # #
          print(paste0("--> Determining optimal SNR aperture radii in O bea",# #
                       "m for HWP angle = ", as, "º..."))                    # #
                                                                             # #
          hist_csv <- paste0(mwSrcAPT_path, fits_name, "-Optimal_aperture_r",# #
                             "adii_log.csv")                                 # #
          opt_snr_cmd <- paste0("err_arr, flux_arr, r_arr, snr_arr = sexD.i",# #
                                "ndiv_max_snr_circ_aperture(data_m = O_temp",# #
                                ", sources = stars, file = '", hist_csv, "'",# #
                                ", data_err = O_temp_err, gain = gain_m)")   # #
          py_run_string(opt_snr_cmd)                                         # #
          rm(opt_snr_cmd, hist_csv)                                          # #
                                                                             # #
          print(paste0("--> Loading optimal SNR aperture radii in O beam fo",# #
                       "r HWP angle = ", as, "º to R Env. object..."))       # #
                                                                             # #
          optO_r_arr <- array(0, dim = length(py$r_arr))                     # #
          optO_snr_arr <- array(0, dim = length(py$snr_arr))                 # #
                                                                             # #
          for(i in 1:dim(optO_r_arr)){                                       # #
            optO_r_arr[i] <- py$r_arr[[i]][1]                                # #
            optO_snr_arr[i] <- py$snr_arr[[i]][1]                            # #
          }                                                                  # #
                                                                             # #
          apt_r <- median(optO_r_arr)                                        # #
          apt_snr <- median(optO_snr_arr)                                    # #
          rm(optO_r_arr, optO_snr_arr, max_r)                                # #
                                                                             # #
          fwhm_r <- median(psf_coefs_ini$r_fwhm, na.rm = T)                  # #
          a2f_ratio <- apt_r / fwhm_r                                        # #
                                                                             # #
          opt_info <- paste0(mwSrcAPT_path, fits_name, "-Optimal_aperture_i",# #
                             "nfo.txt")                                      # #
          cat(paste0("The median values here reported were extracted from a",# #
                     " sample of N = ", length(psf_coefs_ini$r_fwhm), " MW ",# #
                     "stars."), file = opt_info, append = F)                 # #
          cat(paste0("Median optimal SNR = ", apt_snr, "\n"), file =opt_info,# #
              append = T)                                                    # #
          cat(paste0("Median optimal aperture R for HWP angle 0º = ", apt_r, # #
                     "\n"), file = opt_info, append = T)                     # #
          cat(paste0("Median FWHM R for HWP angle 0º = ", fwhm_r, "\n"),     # #
              file = opt_info, append = T)                                   # #
          cat(paste0("Ratio between median optimal aperture R and FWHM (as ",# #
                     "determined for O beam at HWP Angle 0º) = ", a2f_ratio, # #
                     "\n"), file = opt_info, append = T)                     # #
          cat("\n", opt_info, append = T)                                    # #
        }else{                                                               # #
          fwhm_r <- median(psf_coefs_ini$r_fwhm, na.rm = T)                  # #
          apt_r <- a2f_ratio * fwhm_r                                        # #
          cat(paste0("Median optimal aperture R for HWP angle ", as, "º = ", # #
                     apt_r, "\n"), file = opt_info, append = T)              # #
          cat(paste0("Median FWHM R for HWP angle 0º = ", fwhm_r, "\n"),     # #
              file = opt_info, append = T)                                   # #
        }                                                                    # #
        rm(fwhm_r)                                                           # #
                                                                             # #
        print(paste0("----> Optimal aperture radius for HWP angle ", as,     # #
                     "º is ", apt_r, " ..."))                                # #
                                                                             # #
        # Set up opposite beam to cross-check preliminary list of sources    # #
        print(paste0("--> Loading E beam to match list of candidate MW sour",# #
                     "ces for HWP angle ", as, "º..."))                      # #
        temp_data <- mw_E[,, "Data", ang]                                    # #
                                                                             # #
        # List sources in complementary beam                                 # #
        print(paste0("--> Creating list and segmentation map of potential M",# #
                     "W sources in E beam for HWP angle ", as, "º..."))      # #
                                                                             # #
        fin_sep_cmd <- paste0("sources_py, seg_map = sep.extract(E_temp,",   # #
                              src_sig, ", err = E_temp_err, mask = E_mask_t",# #
                              "emp, segmentation_map=True)")                 # #
        py_run_string(fin_sep_cmd)                                           # #
        rm(fin_sep_cmd)                                                      # #
                                                                             # #
        # Save source list to a csv file to be read later                    # #
        print(paste0("--> Saving list of potential MW sources in E beam for",# #
                     " HWP angle = ", as,  "º to CSV file..."))              # #
                                                                             # #
        ### Reject sources based on area and eccentricity values in E beam   # #
        lst_src_cmd <- paste0("rej_obj = sexD.list_sources_info(E_temp, sou",# #
                              "rces_py,'",acp_csv_1c,"', '", rej_csv_1c,"',",# #
                              max_ecc, ",", min_A, ",", max_A, ", E_temp_er",# #
                              "r)")                                          # #
                                                                             # #
        py_run_string(lst_src_cmd)                                           # #
        rm(lst_src_cmd)                                                      # #
                                                                             # #
        # Read csv file and cross-checking sources with previous selection   # #
        print(paste0("--> Reading CSV file with list of potential MW source",# #
                     "s in E beam for HWP angle = ", as, "º to R Env. objec",# #
                     "t..."))                                                # #
                                                                             # #
        srcs_MW <- read.table(acp_csv_1c, header = T, sep = ";", dec = ".")  # #
        rm(acp_csv_1c)                                                       # #
                                                                             # #
        print(paste0("--> Estimating some extra parameters for listed poten",# #
                     "tial MW sources in E beam for HWP angle = ", as, "º..",# #
                     "."))                                                   # #
                                                                             # #
        srcs_MW$r <- (srcs_MW$a + srcs_MW$b) / 2                             # #
        srcs_MW$r_max <- ceiling(sqrt(srcs_MW$Area / pi))                    # #
                                                                             # #
        # Some variables to be used                                          # #
        fin_N <- dim(srcs_MW)[1]                                             # #
        ind_to_rm <- NULL                                                    # #
        psf_coefs_fin <- NULL                                                # #
                                                                             # #
        # Check pixels within listed source's area for bad pixels            # #
        print(paste0("--> Checking each potential source in E beam for HWP ",# #
                     "angle = ", as, "º, for bad pixels, flagging those whi",# #
                     "ch have them and saving psf parameters of those which",# #
                     " don't..."))                                           # #
                                                                             # #
        for(s in 1:fin_N){                                                   # #
          cent_pix <- round(c(srcs_MW$x_barycenter[s],                       # #
                              srcs_MW$y_barycenter[s]))                      # #
          dist <- ceiling(srcs_MW$r_max[s])                                  # #
          xs <- (cent_pix[1] - dist):(cent_pix[1] + dist)                    # #
          ys <- (cent_pix[2] - dist):(cent_pix[2] + dist)                    # #
                                                                             # #
          na_pix_flag <- length(which(is.na(temp_data[xs, ys]), arr.ind = T))# #
                                                                             # #
          if(na_pix_flag){                                                   # #
            ind_to_rm <-append(ind_to_rm, s)                                 # #
            print(paste0("----> Source #", s, " has been flagged..."))       # #
          }                                                                  # #
          else{                                                              # #
            psf_coefs_fin <- rbind(psf_coefs_fin,                            # #
                                   srcs_MW[s,c("Obj_n","Seg_n")])            # #
          }                                                                  # #
        }                                                                    # #
        rm(temp_data, cent_pix, na_pix_flag, dist, fin_N, s)                 # #
                                                                             # #
        ### Remove sources flagged for including NA pixels in E beam         # #
        print(paste0("--> Removing flagged sources in E beam for HWP angle ",# #
                     "= ", as, "º, from list of potential MW sources and fr",# #
                     "om the segmentation map..."))                          # #
                                                                             # #
        rejs_MW <- read.table(rej_csv_1c, header = T, sep = ";", dec = ".")  # #
        rm(rej_csv_1c)                                                       # #
                                                                             # #
        if(!is.null(ind_to_rm)){                                             # #
          rejs_MW <- rbind(rejs_MW, srcs_MW[ind_to_rm, 1:dim(rejs_MW)[2]])   # #
          srcs_MW <- srcs_MW[-ind_to_rm,]                                    # #
        }                                                                    # #
                                                                             # #
        # Saving register of rejected and accepted sources                   # #
        write.table(rejs_MW, rej_csv_2c, row.names = F, sep = ";", dec = ".")# #
        write.table(srcs_MW, acp_csv_2c, row.names = F, sep = ";", dec = ".")# #
        rm(rej_csv_2c, acp_csv_2c, rejs_MW, ind_to_rm)                       # #
                                                                             # #
        print(paste0("--> Estimating PSF parameters for list of potential M",# #
                     "W sources in E beam for HWP angle = ", as, "º..."))    # #
                                                                             # #
        srcs_MW$sig <- srcs_MW$r                                             # #
        srcs_MW$fwhm <- srcs_MW$sig * 2 * sqrt(2 * log(2))                   # #
                                                                             # #
        # This will go into python environment which indexes for fits        # #
        # images differently from R hence the switch between x and y         # #
        psf_coefs_fin$y <- srcs_MW$x_barycenter                              # #
        psf_coefs_fin$x <- srcs_MW$y_barycenter                              # #
        psf_coefs_fin$r_fwhm <- srcs_MW$fwhm / 2                             # #
        psf_coefs_fin <- psf_coefs_fin[, -1]                                 # #
                                                                             # #
        # Matching sources in O beam to sources in E beam                    # #
        print(paste0("--> Matching candidate sources in E beam with those i",# #
                     "n O beam for HWP angle = ", as,  "º, saving those tha",# #
                     "t match and discarding those that don't..."))          # #
                                                                             # #
        ini_srcs <- 1:dim(psf_coefs_ini)[1]                                  # #
        fin_srcs <- 1:dim(psf_coefs_fin)[1]                                  # #
        ind_to_rm <- NULL                                                    # #
                                                                             # #
        srcs_MWo <- read.table(acp_csv_2, header = T, sep = ";", dec = ".")  # #
        rejs_MWo <- read.table(rej_csv_2, header = T, sep = ";", dec = ".")  # #
        rm(rej_csv_2, acp_csv_2)                                             # #
                                                                             # #
        ### Reject sources not matched across beams at the same HWP angle    # #
        for(si in ini_srcs){                                                 # #
          # Tolerance box                                                    # #
          xi <- floor(psf_coefs_ini$x[si] - 1.5)                             # #
          xf <- ceiling(psf_coefs_ini$x[si] + 1.5)                           # #
          yi <- floor(psf_coefs_ini$y[si] - 1.5)                             # #
          yf <- ceiling(psf_coefs_ini$y[si] + 1.5)                           # #
                                                                             # #
          for(sf in fin_srcs){                                               # #
            if(xi <= psf_coefs_fin$x[sf] && xf >= psf_coefs_fin$x[sf] &&     # #
               yi <= psf_coefs_fin$y[sf] && yf >= psf_coefs_fin$y[sf]){      # #
                                                                             # #
              # If a source is found to match tolerance box, the parameters  # #
              # of that source in each beam are saved to a new data frame,   # #
              # that source is removed from search list in complementary beam# #
              # and the search resumes in for a new source of first beam     # #
              temp_ini <- psf_coefs_ini[si, c('x','y','Seg_n')]              # #
              temp_fin <- psf_coefs_fin[sf, c('x', 'y')]                     # #
                                                                             # #
              switch(ang,                                                    # #
                     c(cdd_srcs_0$O <- rbind(cdd_srcs_0$O, temp_ini),        # #
                       cdd_srcs_0$E <- rbind(cdd_srcs_0$E, temp_fin)),       # #
                     c(cdd_srcs_22$O <- rbind(cdd_srcs_22$O, temp_ini),      # #
                       cdd_srcs_22$E <- rbind(cdd_srcs_22$E, temp_fin)),     # #
                     c(cdd_srcs_45$O <- rbind(cdd_srcs_45$O, temp_ini),      # #
                       cdd_srcs_45$E <- rbind(cdd_srcs_45$E, temp_fin)),     # #
                     c(cdd_srcs_68$O <- rbind(cdd_srcs_68$O, temp_ini),      # #
                       cdd_srcs_68$E <- rbind(cdd_srcs_68$E, temp_fin)))     # #
                                                                             # #
              fin_srcs <- fin_srcs[-which(fin_srcs == sf, arr.ind = TRUE)]   # #
              break                                                          # #
            }else{                                                           # #
              if(sf == fin_srcs[length(fin_srcs)]){                          # #
                ind_to_rm <- append(ind_to_rm, si)                           # #
                                                                             # #
                if(ang == 1){                                                # #
                  tmp_seg_map[which(tmp_seg_map == psf_coefs_ini[si,'Seg_n'],# #
                                    arr.ind = TRUE)] <- 0                    # #
                }                                                            # #
              }                                                              # #
            }                                                                # #
          }                                                                  # #
        }                                                                    # #
                                                                             # #
        if(!is.null(ind_to_rm)){                                             # #
          rejs_MWo <- rbind(rejs_MWo,                                        # #
                            srcs_MWo[ind_to_rm, 1:dim(rejs_MWo)[2]])         # #
          srcs_MWo <- srcs_MWo[-ind_to_rm,]                                  # #
        }                                                                    # #
        rm(ind_to_rm)                                                        # #
                                                                             # #
        write.table(rejs_MWo, rej_csv_3, row.names = F, sep = ";", dec = ".")# #
        write.table(srcs_MWo, acp_csv_3, row.names = F, sep = ";", dec = ".")# #
        rm(rej_csv_3, acp_csv_3)                                             # #
                                                                             # #
        if(ang == 1){                                                        # #
          seg_name <- paste0(seg_path, "acp3_Seg_Map-band=", bs, "_set#=", n,# #
                             run_date, ".fits")                              # #
          writeFITSim(tmp_seg_map, seg_name)                                 # #
          rm(seg_name)                                                       # #
                                                                             # #
          srcs_MW0 <- srcs_MWo                                               # #
          rejs_MW0 <- rejs_MWo                                               # #
        }                                                                    # #
                                                                             # #
        src_map_ang[,,ang] <- tmp_seg_map                                    # #
                                                                             # #
        rm(psf_coefs_ini, psf_coefs_fin, fin_srcs, temp_ini, temp_fin, si,   # #
           ini_srcs, srcs_MWo, rejs_MWo)                                     # #
      }                                                                      # #
      py_run_string("del(E_temp, O_temp, E_temp_err, O_temp_err, rej_obj)")  # #
      py_run_string("del(E_mask_temp, O_mask_temp, sources_py, seg_map)")    # #
      rm(Es_pix_mask, Os_pix_mask)                                           # #
                                                                             # #
      # Matching sources in HWP = 0º to sources in other HWP angles          # #
      print(paste0("* Matching sources across all different HWP angles for ",# #
                   "set #", n, " of band ", bs, ", saving those that match ",# #
                   "and discarding those that don't *"))                     # #
                                                                             # #
      srcs_0 <- NULL                                                         # #
      srcs_0$E <- NULL                                                       # #
      srcs_0$O <- NULL                                                       # #
      srcs_22 <- NULL                                                        # #
      srcs_22$E <- NULL                                                      # #
      srcs_22$O <- NULL                                                      # #
      srcs_45 <- NULL                                                        # #
      srcs_45$E <- NULL                                                      # #
      srcs_45$O <- NULL                                                      # #
      srcs_68 <- NULL                                                        # #
      srcs_68$E <- NULL                                                      # #
      srcs_68$O <- NULL                                                      # #
                                                                             # #
      ind_to_rm <- NULL                                                      # #
                                                                             # #
      # Check in only one beam since the complementary                       # #
      # beams were checked in the previous sub-routine                       # #
      ls_0 <- 1:dim(cdd_srcs_0$O)[1]                                         # #
      ls_22 <- 1:dim(cdd_srcs_22$O)[1]                                       # #
      ls_45 <- 1:dim(cdd_srcs_45$O)[1]                                       # #
      ls_68 <- 1:dim(cdd_srcs_68$O)[1]                                       # #
      t_0 <- cdd_srcs_0$O                                                    # #
      t_22 <- cdd_srcs_22$O                                                  # #
      t_45 <- cdd_srcs_45$O                                                  # #
      t_68 <- cdd_srcs_68$O                                                  # #
                                                                             # #
      for(s0 in ls_0){                                                       # #
        print(paste0("--> Checking source #", s0, " of ", max(ls_0), "...")) # #
                                                                             # #
        # s0 will be compared against sources in all other 3 HWP angles      # #
        # flg_aa will register in which other angles s0 is matched,          # #
        # if s0 is matched in all angles the sources is saved to final list  # #
        flg_22 <- FALSE                                                      # #
        flg_45 <- FALSE                                                      # #
        flg_68 <- FALSE                                                      # #
                                                                             # #
        pt0 <- c(t_0[s0, 'x'], t_0[s0, 'y'])                                 # #
                                                                             # #
        # In each of the next 3 loops                                        # #
        # If a source is found to match the tolerance box it is flagged,     # #
        # removed from search list and the loop breaks                       # #
        for(s22 in ls_22){                                                   # #
          pt22 <- c(t_22[s22, 'x'], t_22[s22, 'y'])                          # #
                                                                             # #
          if(dist(rbind(pt0, pt22)) <= 2){                                   # #
            flg_22 <- TRUE                                                   # #
            ls_22 <- ls_22[-which(ls_22 == s22, arr.ind = TRUE)]             # #
            break                                                            # #
          }                                                                  # #
        }                                                                    # #
        for(s45 in ls_45){                                                   # #
          pt45 <- c(t_45[s45, 'x'], t_45[s45, 'y'])                          # #
                                                                             # #
          if(dist(rbind(pt0, pt45)) <= 2){                                   # #
            flg_45 <- TRUE                                                   # #
            ls_45 <- ls_45[-which(ls_45 == s45, arr.ind = TRUE)]             # #
            break                                                            # #
          }                                                                  # #
        }                                                                    # #
        for(s68 in ls_68){                                                   # #
          pt68 <- c(t_68[s68, 'x'], t_68[s68, 'y'])                          # #
                                                                             # #
          if(dist(rbind(pt0, pt68)) <= 2){                                   # #
            flg_68 <- TRUE                                                   # #
            ls_68 <- ls_68[-which(ls_68 == s68, arr.ind = TRUE)]             # #
            break                                                            # #
          }                                                                  # #
        }                                                                    # #
        rm(pt22, pt45, pt68)                                                 # #
                                                                             # #
        # Adding candidate to final list if it was successfully matched      # #
        flg_0 <- flg_22 & flg_45 & flg_68                                    # #
        rm(flg_22, flg_45, flg_68)                                           # #
                                                                             # #
        if(flg_0){                                                           # #
          srcs_0$O <- rbind(srcs_0$O, cdd_srcs_0$O[s0, c('x','y')])          # #
          srcs_0$E <- rbind(srcs_0$E, cdd_srcs_0$E[s0, c('x','y')])          # #
          srcs_22$O <- rbind(srcs_22$O, cdd_srcs_22$O[s22, c('x','y')])      # #
          srcs_22$E <- rbind(srcs_22$E, cdd_srcs_22$E[s22, c('x','y')])      # #
          srcs_45$O <- rbind(srcs_45$O, cdd_srcs_45$O[s45, c('x','y')])      # #
          srcs_45$E <- rbind(srcs_45$E, cdd_srcs_45$E[s45, c('x','y')])      # #
          srcs_68$O <- rbind(srcs_68$O, cdd_srcs_68$O[s68, c('x','y')])      # #
          srcs_68$E <- rbind(srcs_68$E, cdd_srcs_68$E[s68, c('x','y')])      # #
                                                                             # #
          print(paste0("----> Source ", s0, " of ", max(ls_0), " has been a",# #
                       "ddded to the list of sources common to all HWP angl",# #
                       "es and beams..."))                                   # #
        }else{                                                               # #
          print(paste0("----> Source ", s0, " of ", max(ls_0), " has been r",# #
                       "ejected..."))                                        # #
                                                                             # #
          ind_to_rm <- append(ind_to_rm, s0)                                 # #
                                                                             # #
          tmp_seg_map[which(tmp_seg_map == t_0[s0,'Seg_n'],                  # #
                            arr.ind = T)] <- 0                               # #
        }                                                                    # #
        rm(flg_0)                                                            # #
      }                                                                      # #
      rm(cdd_srcs_0, cdd_srcs_22, cdd_srcs_45, cdd_srcs_68, xf, xi, yf, yi,  # #
         ls_0, ls_22, ls_45, ls_68, t_0, t_22, t_45, t_68, s0, s22, s45, s68)# #
                                                                             # #
      ### Reject sources not matched across different HWP angles             # #
      if(!is.null(ind_to_rm)){                                               # #
        rejs_MW <- rbind(rejs_MW0, srcs_MW0[ind_to_rm, 1:dim(rejs_MW0)[2]])  # #
        srcs_MW <- srcs_MW0[-ind_to_rm,]                                     # #
      }else{                                                                 # #
        rejs_MW <- rejs_MW0                                                  # #
        srcs_MW <- srcs_MW0                                                  # #
      }                                                                      # #
      rm(ind_to_rm, rejs_MW0, srcs_MW0)                                      # #
                                                                             # #
      # Saving register of rejected and accepted sources                     # #
      write.table(rejs_MW, rej_csv_4, row.names = F, sep = ";", dec = ".")   # #
      write.table(srcs_MW, acp_csv_4, row.names = F, sep = ";", dec = ".")   # #
      rm(rej_csv_4, acp_csv_4)                                               # #
                                                                             # #
      # Save segmentation map to FITS file                                   # #
      print(paste0("* Saving segmentation map, showing Milky Way sources co",# #
                   "mmon to E and O beams at all HWP angles for set #", n,   # #
                   " of band ", bs, ", to FITS file. This file can be found",# #
                   " in ", seg_path, " *"))                                  # #
                                                                             # #
      seg_name <- paste0(seg_path, "acp4_Seg_Map-band=", bs, "_set#=", n,    # #
                         run_date, ".fits")                                  # #
      writeFITSim(tmp_seg_map, seg_name)                                     # #
                                                                             # #
      src_map_sum <- apply(src_map_ang, 1:2, sum)                            # #
      rm(src_map_ang, seg_name)                                              # #
                                                                             # #
      print(paste0("* Calculating and storing E and O beams fluxes for sele",# #
                   "cted sources at all different HWP angles for set #", n,  # #
                   " of band ", bs, ", saving those that match and discardi",# #
                   "ng those that don't *"))                                 # #
                                                                             # #
      n_srcs <- dim(srcs_0$E)[1]                                             # #
                                                                             # #
      src_flux_dim <- c(n_srcs, length(beam_par), length(munc_par),          # #
                        length(ang_str))                                     # #
      src_flux_dname <- list(NULL, beam_par, munc_par, ang_str)              # #
                                                                             # #
      mw_srcs_flx <- array(NA, dim = src_flux_dim, dimnames = src_flux_dname)# #
                                                                             # #
      for(ang in angs){                                                      # #
        as <- ang_str[ang]                                                   # #
                                                                             # #
        apt_rs <- NULL                                                       # #
                                                                             # #
        # Load E beam and O beam maps to python environment                  # #
        print(paste0("--> Loading, to Python Env., E and O beam data, with ",# #
                     "extended source background level removed, and uncerta",# #
                     "inties of HWP angle = ", as, "º..."))                  # #
                                                                             # #
        temp <- mw_E[,, "Data", ang]                                         # #
        py_run_string("E_temp = r.temp.copy(order='C')")                     # #
        temp <- mw_O[,, "Data", ang]                                         # #
        py_run_string("O_temp = r.temp.copy(order='C')")                     # #
        temp <- mw_E[,, "Unc", ang]                                          # #
        py_run_string("E_temp_err = r.temp.copy(order='C')")                 # #
        temp <- mw_O[,, "Unc", ang]                                          # #
        py_run_string("O_temp_err = r.temp.copy(order='C')")                 # #
        rm(temp)                                                             # #
                                                                             # #
        switch(ang,                                                          # #
               c(Ex <- srcs_0$E[,'x'], Ox <- srcs_0$O[,'x'],                 # #
                 Ey <- srcs_0$E[,'y'], Oy <- srcs_0$O[,'y']),                # #
               c(Ex <- srcs_22$E[,'x'], Ox <- srcs_22$O[,'x'],               # #
                 Ey <- srcs_22$E[,'y'], Oy <- srcs_22$O[,'y']),              # #
               c(Ex <- srcs_45$E[,'x'], Ox <- srcs_45$O[,'x'],               # #
                 Ey <- srcs_45$E[,'y'], Oy <- srcs_45$O[,'y']),              # #
               c(Ex <- srcs_68$E[,'x'], Ox <- srcs_68$O[,'x'],               # #
                 Ey <- srcs_68$E[,'y'], Oy <- srcs_68$O[,'y'])               # #
        )                                                                    # #
                                                                             # #
        print(paste0("--> Loading, to Python Env., optimal aperture radius ",# #
                     "for HWP angle = ", as, "º..."))                        # #
                                                                             # #
        apt_rs <- c(apt_rs, apt_r)                                           # #
                                                                             # #
        temp <- apt_r                                                        # #
        py_run_string("r_opt = np.array(r.temp, ndmin = 1)")                 # #
                                                                             # #
        print(paste0("---> The optimal aperture radius for HWP angle = ", as,# #
                     "º is ", temp, " pix..."))                              # #
        rm(temp)                                                             # #
                                                                             # #
        # Calculating fluxes of sources accepted until this point            # #
        for(ns in 1:n_srcs){                                                 # #
          # Load coordinates for present source to python env                # #
          print(paste0("----> Loading, to Python Env., coordinates of sourc",# #
                       "e #", ns, " of ", n_srcs, " in E and O beam maps fo",# #
                       "r HWP angle = ", as, "º..."))                        # #
                                                                             # #
          # python treats arrays different than R, this file will be read by # #
          # R, as such x and y (rows and columns) must be switched around    # #
          temp <- Ex[ns]                                                     # #
          py_run_string("xn_E = np.array(r.temp, ndmin = 1)")                # #
          temp <- Ey[ns]                                                     # #
          py_run_string("yn_E = np.array(r.temp, ndmin = 1)")                # #
          temp <- Ox[ns]                                                     # #
          py_run_string("xn_O = np.array(r.temp, ndmin = 1)")                # #
          temp <- Oy[ns]                                                     # #
          py_run_string("yn_O = np.array(r.temp, ndmin = 1)")                # #
          rm(temp)                                                           # #
                                                                             # #
          print(paste0("----> Calculating fluxes of source #", ns, " of ",   # #
                       n_srcs, " in E and O beam maps for HWP angle = ", as, # #
                       "º..."))                                              # #
                                                                             # #
          calc_flx_E_cmd <- paste0("E_flx, E_flx_U, E_flg = sep.sum_circle(",# #
                                   "data = E_temp, x = xn_E, y = yn_E, r = ",# #
                                   "r_opt, err = E_temp_err, subpix = 0, ga",# #
                                   "in = gain_m)")                           # #
          calc_flx_O_cmd <- paste0("O_flx, O_flx_U, O_flg = sep.sum_circle(",# #
                                   "data = O_temp, x = xn_O, y = yn_O, r = ",# #
                                   "r_opt, err = O_temp_err, subpix = 0, ga",# #
                                   "in = gain_m)")                           # #
          py_run_string(calc_flx_E_cmd)                                      # #
          py_run_string(calc_flx_O_cmd)                                      # #
          rm(calc_flx_E_cmd, calc_flx_O_cmd)                                 # #
                                                                             # #
          print(paste0("----> Saving fluxes and uncertainties of source #",  # #
                       ns, " of ", n_srcs, " in E and O beam for HWP angle ",# #
                       "= ", as, "º to R Env. object..."))                   # #
                                                                             # #
          mw_srcs_flx[ns, 'E', "Data", as] <- py$E_flx                       # #
          mw_srcs_flx[ns, 'O', "Data", as] <- py$O_flx                       # #
                                                                             # #
          if(is.na(py$E_flx_U)){                                             # #
            tnoise2 <- median(mw_E[,,"Unc",ang], na.rm = T)^2 * pi * apt_r^2 # #
            mw_srcs_flx[ns,'E',"Unc",as] <- sqrt(tnoise2 + py$E_flx / gain)  # #
            rm(tnoise2)                                                      # #
          }else{                                                             # #
            mw_srcs_flx[ns,'E',"Unc",as] <- py$E_flx_U                       # #
          }                                                                  # #
          if(is.na(py$O_flx_U)){                                             # #
            tnoise2 <- median(mw_O[,,"Unc",ang], na.rm = T)^2 * pi * apt_r^2 # #
            mw_srcs_flx[ns,'O',"Unc",as] <- sqrt(tnoise2 + py$O_flx / gain)  # #
            rm(tnoise2)                                                      # #
          }else{                                                             # #
            mw_srcs_flx[ns,'O',"Unc",as] <- py$O_flx_U                       # #
          }                                                                  # #
                                                                             # #
          print(paste0("----> Flux for source #", ns, " of ", n_srcs, " in ",# #
                       "E beam for HWP angle = ", as,  "º is ", py$E_flx,    # #
                       " +- ", mw_srcs_flx[ns,'E',"Unc",as], "..."))         # #
          print(paste0("----> Flux for source #", ns, " of ", n_srcs, " in ",# #
                       "O beam for HWP angle = ", as,  "º is ", py$O_flx,    # #
                       " +- ", mw_srcs_flx[ns,'O',"Unc",as], "..."))         # #
        }                                                                    # #
      }                                                                      # #
      rm(mw_E, mw_O)                                                         # #
                                                                             # #
      apt_rad_arr[b, n] <- max(apt_rs, na.rm = T)                            # #
      rm(apt_rs)                                                             # #
                                                                             # #
      ### Reject sources with NA flux or uncertainty in 1 or more HWP angles # #
      print(paste0("* Rejecting sources with NaN as the value for either fl",# #
                   "ux or flux uncertainty at any of the beams and HWP angl",# #
                   "es for set#", n, " of band ", bs, ". *"))                # #
                                                                             # #
      ind_to_rm <- NULL                                                      # #
      for(ns in 1:n_srcs){                                                   # #
        if(length(which(is.na(mw_srcs_flx[ns,,,]), arr.ind = T))){           # #
          ind_to_rm <- append(ind_to_rm, ns)                                 # #
        }                                                                    # #
      }                                                                      # #
                                                                             # #
      if(!is.null(ind_to_rm)){                                               # #
        rejs_MW <- rbind(rejs_MW, srcs_MW[ind_to_rm, 1:dim(rejs_MW)[2]])     # #
        srcs_MW <- srcs_MW[-ind_to_rm,]                                      # #
        mw_srcs_flx <- mw_srcs_flx[-ind_to_rm,,,]                            # #
        srcs_0$O <- srcs_0$O[-ind_to_rm,]                                    # #
        srcs_22$E <- srcs_22$E[-ind_to_rm,]                                  # #
        srcs_22$O <- srcs_22$O[-ind_to_rm,]                                  # #
        srcs_45$E <- srcs_45$E[-ind_to_rm,]                                  # #
        srcs_45$O <- srcs_45$O[-ind_to_rm,]                                  # #
        srcs_68$E <- srcs_68$E[-ind_to_rm,]                                  # #
        srcs_68$O <- srcs_68$O[-ind_to_rm,]                                  # #
                                                                             # #
        for(i in ind_to_rm){                                                 # #
          src_numb <- tmp_seg_map[srcs_0$E[i,"y"], srcs_0$E[i,"x"]]          # #
          tmp_seg_map[which(tmp_seg_map == src_numb, arr.ind = T)] <- 0      # #
        }                                                                    # #
                                                                             # #
        srcs_0$E <- srcs_0$E[-ind_to_rm,]                                    # #
      }                                                                      # #
      rm(ind_to_rm)                                                          # #
                                                                             # #
      # Saving register of rejected and accepted sources                     # #
      write.table(rejs_MW, rej_csv_f, row.names = F, sep = ";", dec = ".")   # #
      write.table(srcs_MW, acp_csv_f, row.names = F, sep = ";", dec = ".")   # #
      rm(rej_csv_f, acp_csv_f, rejs_MW, srcs_MW)                             # #
                                                                             # #
      # Save segmentation map to FITS file                                   # #
      print(paste0("* Saving segmentation map, showing Milky Way sources wi",# #
                   "th valid fluxes at all beams and HWP angles for set #",  # #
                   n, " of band ", bs, ", to FITS file. This file can be fo",# #
                   "und in ", seg_path, " *"))                               # #
                                                                             # #
      seg_name <- paste0(seg_path, "Fin_Seg_Map-band=", bs, "_set#=", n,     # #
                         run_date, ".fits")                                  # #
                                                                             # #
      writeFITSim(tmp_seg_map, seg_name)                                     # #
      rm(seg_name, tmp_seg_map)                                              # #
                                                                             # #
      # Saving flux estimates of selected sources                            # #
      write.table(mw_srcs_flx, file = mwflux_out_csv, row.names = F,         # #
                  sep = ";", dec = ".")                                      # #
      rm(mwflux_out_csv)                                                     # #
                                                                             # #
      py_run_string("del(E_temp, O_temp, E_temp_err, O_temp_err,r_opt,xn_E)")# #
      py_run_string("del(yn_E, xn_O, yn_O, E_flx, E_flx_U, E_flg, O_flx)")   # #
      py_run_string("del(O_flx_U, O_flg)")                                   # #
      rm(Ex, Ox, Ey, Oy)                                                     # #
                                                                             # #
      # Creating an array to hold all sources coordinates                    # #
      # across beams and HWP angles                                          # #
      n_srcs <- dim(srcs_0$E)[1]                                             # #
                                                                             # #
      src_cd_dim <- c(n_srcs, length(beam_ang_par), length(pos_par))         # #
      src_cd_dname <- list(NULL, beam_ang_par, pos_par)                      # #
                                                                             # #
      srcs_coords <- array(NA, dim = src_cd_dim, dimnames = src_cd_dname)    # #
      rm(src_cd_dim, src_cd_dname)                                           # #
                                                                             # #
      # Loading coordinates to new array                                     # #
      srcs_coords[,"E0","y"] <- srcs_0$E[,'x']                               # #
      srcs_coords[,"O0","y"] <- srcs_0$O[,'x']                               # #
      srcs_coords[,"E0","x"] <- srcs_0$E[,'y']                               # #
      srcs_coords[,"O0","x"] <- srcs_0$O[,'y']                               # #
      srcs_coords[,"E22","y"] <- srcs_22$E[,'x']                             # #
      srcs_coords[,"O22","y"] <- srcs_22$O[,'x']                             # #
      srcs_coords[,"E22","x"] <- srcs_22$E[,'y']                             # #
      srcs_coords[,"O22","x"] <- srcs_22$O[,'y']                             # #
      srcs_coords[,"E45","y"] <- srcs_45$E[,'x']                             # #
      srcs_coords[,"O45","y"] <- srcs_45$O[,'x']                             # #
      srcs_coords[,"E45","x"] <- srcs_45$E[,'y']                             # #
      srcs_coords[,"O45","x"] <- srcs_45$O[,'y']                             # #
      srcs_coords[,"E68","y"] <- srcs_68$E[,'x']                             # #
      srcs_coords[,"O68","y"] <- srcs_68$O[,'x']                             # #
      srcs_coords[,"E68","x"] <- srcs_68$E[,'y']                             # #
      srcs_coords[,"O68","x"] <- srcs_68$O[,'y']                             # #
      rm(srcs_0, srcs_22, srcs_45, srcs_68)                                  # #
                                                                             # #
      # Creating new array to hold average position of sources               # #
      src_xy_dim <- c(n_srcs, length(pos_par))                               # #
      src_xy_dname <- list(NULL, pos_par)                                    # #
                                                                             # #
      srcs_xy <- array(NA, dim = src_xy_dim, dimnames = src_xy_dname)        # #
      rm(src_xy_dim, src_xy_dname)                                           # #
                                                                             # #
      # Calculating average position of sources across beams and HWP angles  # #
      if(n_srcs == 1){                                                       # #
        srcs_xy[,"x"] <- round(median(srcs_coords[1,,"x"]))                  # #
        srcs_xy[,"y"] <- round(median(srcs_coords[1,,"y"]))                  # #
      }else{                                                                 # #
        srcs_xy[,"x"] <- round(apply(srcs_coords[,,"x"], 1, median))         # #
        srcs_xy[,"y"] <- round(apply(srcs_coords[,,"y"], 1, median))         # #
      }                                                                      # #
      rm(srcs_coords)                                                        # #
                                                                             # #
      # Creating a new array that will hold CCD positions and polarimetry    # #
      # for the selected sources                                             # #
      srcs_dim <- c(n_srcs, length(srcs_par))                                # #
      srcs_dname <- list(NULL, srcs_par)                                     # #
                                                                             # #
      mw_srcs_QUPX <- array(NA, dim = srcs_dim, dimnames = srcs_dname)       # #
      rm(srcs_dim, srcs_dname)                                               # #
                                                                             # #
      mw_srcs_QUPX[,'x'] <- srcs_xy[,'x']                                    # #
      mw_srcs_QUPX[,'y'] <- srcs_xy[,'y']                                    # #
      ##########################################################################
                                                                               #
      #----------------------------------------------------------------------# #
                                                                               #
      ################## Foreground Polarization Estimation ################## #
      ##########################################################################
      src_pol_dim <- c(n_srcs, length(stokes_par), length(munc_par))         # #
      mw_pol_dim <- c(length(stokes_par), length(munc_par))                  # #
      mw_polV_dim <- c(length(polV_par), length(munc_par))                   # #
      src_pol_dname <- list(NULL, stokes_par, munc_par)                      # #
      mw_pol_dname<- list(stokes_par, munc_par)                              # #
      mw_polV_dname<- list(polV_par, munc_par)                               # #
                                                                             # #
      mw_srcs_pol <- array(NA, dim = src_pol_dim, dimnames = src_pol_dname)  # #
      mw_srcs_pol_cc <- array(NA, dim = src_pol_dim, dimnames =src_pol_dname)# #
      mw_srcs_pol_i <- array(NA, dim = src_pol_dim, dimnames = src_pol_dname)# #
      mw_QUr <- array(NA, dim = mw_pol_dim, dimnames = mw_pol_dname)         # #
      mw_QUcc <- array(NA, dim = mw_pol_dim, dimnames = mw_pol_dname)        # #
      mw_QUi <- array(NA, dim = mw_pol_dim, dimnames = mw_pol_dname)         # #
      mw_PXr <- array(NA, dim =  mw_polV_dim, dimnames = mw_polV_dname)      # #
      mw_PXcc <- array(NA, dim = mw_polV_dim, dimnames = mw_polV_dname)      # #
      mw_PXi <- array(NA, dim = mw_polV_dim, dimnames = mw_polV_dname)       # #
      mw_PXdeb <- array(NA, dim = mw_polV_dim, dimnames = mw_polV_dname)     # #
                                                                             # #
      # Calculation of Q and U parameters for selected sources               # #
      print(paste0("* Calculating Stokes Q and U of selected Milky Way sour",# #
                   "ces, using the ratio method [Bagnulo et al. 2009], of s",# #
                   "et #", n, " of band ", bs, " *"))                        # #
                                                                             # #
      mw_srcs_pol[,'Q',] <- stokesPar_from_dualBeams_Bagnulo(                # #
        mw_srcs_flx[,'O', "Data", c("0", "45")],                             # #
        mw_srcs_flx[,'E', "Data", c("0", "45")],                             # #
        mw_srcs_flx[,'O', "Unc", c("0", "45")],                              # #
        mw_srcs_flx[,'E', "Unc", c("0", "45")])                              # #
                                                                             # #
      mw_srcs_pol[,'U',] <- stokesPar_from_dualBeams_Bagnulo(                # #
        mw_srcs_flx[,'O', "Data", c("22.5", "67.5")],                        # #
        mw_srcs_flx[,'E', "Data", c("22.5", "67.5")],                        # #
        mw_srcs_flx[,'O', "Unc", c("22.5", "67.5")],                         # #
        mw_srcs_flx[,'E', "Unc", c("22.5", "67.5")])                         # #
      rm(mw_srcs_flx)                                                        # #
                                                                             # #
      # Output CSV with MW sources pol with sky correction                   # #
      write.table(mw_srcs_pol, src_polr_out_csv, sep = ";", dec = ".")       # #
      rm(src_polr_out_csv)                                                   # #
                                                                             # #
      # Plot Q vs U with sky correction of MW sources                        # #
      print(paste0("* Creating PDF file with raw Q vs U plot of selected Mi",# #
                   "lky Way sources for set #", n, " of band ", bs, " *"))   # #
                                                                             # #
      if(is.null(dim(mw_srcs_pol[,,"Data"]))){                               # #
        mw_QUr[,"Data"] <- mw_srcs_pol[,,"Data"]                             # #
        mw_QUr[,"Unc"] <- mw_srcs_pol[,,"Unc"]                               # #
      }else{                                                                 # #
        mw_QUr[,"Data"] <- apply(mw_srcs_pol[,,"Data"], 2, median, na.rm = T)# #
        mw_QUr[,"Unc"] <- unc_median(mw_srcs_pol[,,"Data"],                  # #
                                     mw_srcs_pol[,,"Unc"], 2)                # #
      }                                                                      # #
                                                                             # #
      temp <- P_from_QU(mw_QUr['Q', "Data"], mw_QUr['U', "Data"],            # #
                        mw_QUr['Q', "Unc"], mw_QUr['U', "Unc"])              # #
      mw_PXr['P', "Data"] <- temp$Data                                       # #
      mw_PXr['P', "Unc"] <- temp$Unc                                         # #
      temp <- X_from_QU(mw_QUr['Q', "Data"], mw_QUr['U', "Data"],            # #
                        mw_QUr['Q', "Unc"], mw_QUr['U', "Unc"])              # #
      mw_PXr['X', "Data"] <- temp$Data                                       # #
      mw_PXr['X', "Unc"] <- temp$Unc                                         # #
      rm(temp)                                                               # #
                                                                             # #
      # Output CSV with MW pol parameters with sky correction                # #
      write.table(mw_QUr, mwqur_out_csv, sep = ";", dec = ".")               # #
      write.table(mw_PXr, mwpxr_out_csv, sep = ";", dec = ".")               # #
      rm(mw_PXr, mwqur_out_csv, mwpxr_out_csv)                               # #
                                                                             # #
      # Plot Q vs U with sky correction for all MW sources                   # #
      qlim <- c(min(c(mw_srcs_pol[,'Q',"Data"] - mw_srcs_pol[,'Q',"Unc"],    # #
                      mw_QUr['Q',"Data"] - mw_QUr['Q',"Unc"]), na.rm = T),   # #
                max(c(mw_srcs_pol[,'Q',"Data"] + mw_srcs_pol[,'Q',"Unc"],    # #
                      mw_QUr['Q',"Data"] + mw_QUr['Q',"Unc"]), na.rm =T))*100# #
      ulim <- c(min(c(mw_srcs_pol[,'U',"Data"] - mw_srcs_pol[,'U',"Unc"],    # #
                      mw_QUr['U',"Data"] - mw_QUr['U',"Unc"]), na.rm = T),   # #
                max(c(mw_srcs_pol[,'U',"Data"] + mw_srcs_pol[,'U',"Unc"],    # #
                      mw_QUr['U',"Data"] + mw_QUr['U',"Unc"]), na.rm =T))*100# #
                                                                             # #
      pdf(mwqur_out_pdf, width = 15, height = 15)                            # #
      rm(mwqur_out_pdf)                                                      # #
                                                                             # #
      plot(mw_srcs_pol[,'Q', "Data"] * 100, mw_srcs_pol[,'U', "Data"] * 100, # #
           xlab = 'q (%)', ylab = 'u (%)', xlim = qlim, ylim = ulim,         # #
           pch = 20, cex = 2, col = 'black',                                 # #
           main = "Raw Q vs U of selected Milky Way sources")                # #
      points(mw_QUr['Q', "Data"] * 100, mw_QUr['U', "Data"] * 100, pch = 14, # #
             cex = 3, col = "red")                                           # #
      arrows(mw_srcs_pol[,'Q', "Data"] * 100,                                # #
             y0 = (mw_srcs_pol[,'U', "Data"] - mw_srcs_pol[,'U', "Unc"])*100,# #
             y1 = (mw_srcs_pol[,'U', "Data"] + mw_srcs_pol[,'U', "Unc"])*100,# #
             code = 3, length = 0.02, angle = 90)                            # #
      arrows(x0 = (mw_srcs_pol[,'Q', "Data"] - mw_srcs_pol[,'Q', "Unc"])*100,# #
             mw_srcs_pol[,'U', "Data"] * 100,                                # #
             x1 = (mw_srcs_pol[,'Q', "Data"] + mw_srcs_pol[,'Q', "Unc"])*100,# #
             code = 3, length = 0.02, angle = 90)                            # #
      arrows(mw_QUr['Q',"Data"] * 100,                                       # #
             y0 = (mw_QUr['U',"Data"] - mw_QUr['U',"Unc"]) * 100,            # #
             y1 = (mw_QUr['U',"Data"] + mw_QUr['U',"Unc"]) * 100,            # #
             code = 3, length = 0.03, angle = 90, col = "red")               # #
      arrows(x0 = (mw_QUr['Q',"Data"] - mw_QUr['Q',"Unc"]) * 100,            # #
             mw_QUr['U',"Data"] * 100,                                       # #
             x1 = (mw_QUr['Q',"Data"] + mw_QUr['Q',"Unc"]) * 100,            # #
             code = 3, length = 0.03, angle = 90, col = "red")               # #
      dev.off()                                                              # #
      rm(qlim, ulim)                                                         # #
                                                                             # #
      # Applying instrumental spatial correction to MW sources Stokes Q & U  # #
      print(paste0("* Applying instrumental correction to Stokes Q and U of",# #
                   " selected Milky Way sources of set #", n, " of band ",bs,# #
                   " *"))                                                    # #
                                                                             # #
      rid <- round(apt_rad_arr[b, n] * sqrt(2))                              # #
                                                                             # #
      for(ns in 1:n_srcs){                                                   # #
        xs <- (srcs_xy[ns, "x"] - rid):(srcs_xy[ns, "x"] + rid)              # #
        ys <- (srcs_xy[ns, "y"] - rid):(srcs_xy[ns, "y"] + rid)              # #
        nsQ <- mw_srcs_pol[ns,'Q',]                                          # #
        nsU <- mw_srcs_pol[ns,'U',]                                          # #
        instQns <- median(Qinst[xs, ys, "Data"], na.rm = T)                  # #
        instUns <- median(Uinst[xs, ys, "Data"], na.rm = T)                  # #
        instQns_unc <- unc_median(Qinst[xs,ys,"Data"], Qinst[xs,ys,"Unc"], 0)# #
        instUns_unc <- unc_median(Uinst[xs,ys,"Data"], Uinst[xs,ys,"Unc"], 0)# #
                                                                             # #
        mw_srcs_pol_i[ns,'Q',"Data"] <- nsQ["Data"] - instQns                # #
        mw_srcs_pol_i[ns,'U',"Data"] <- nsU["Data"] - instUns                # #
        mw_srcs_pol_i[ns,'Q',"Unc"] <- unc_add(nsQ["Unc"], instQns_unc)      # #
        mw_srcs_pol_i[ns,'U',"Unc"] <- unc_add(nsU["Unc"], instUns_unc)      # #
      }                                                                      # #
      rm(xs,ys,ns, nsQ, nsU, instQns, instUns, instQns_unc, instUns_unc, rid)# #
                                                                             # #
      # Output CSV with MW sources pol parameters with sky                   # 3
      # and instrumental spatial corrections                                 # #
      write.table(mw_srcs_pol_i, src_poli_out_csv, sep = ";", dec = ".")     # #
      rm(src_poli_out_csv)                                                   # #
                                                                             # #
      # Estimating Stokes Q and U from MW with sky                           # #
      # and spatial instrumental correction                                  # #
      print(paste0("* Calculating median Stokes Q and U of selected Milky W",# #
                   "ay sources with instrumental correction, for an estimat",# #
                   "e of the foreground polarization level of set #", n, " ",# #
                   "of band ", bs, " *"))                                    # #
                                                                             # #
      if(is.null(dim(mw_srcs_pol_i[,,"Data"]))){                             # #
        mw_QUi[,"Data"] <- mw_srcs_pol_i[,,"Data"]                           # #
        mw_QUi[,"Unc"] <- mw_srcs_pol_i[,,"Unc"]                             # #
      }else{                                                                 # #
        mw_QUi[,"Data"] <- apply(mw_srcs_pol_i[,,"Data"], 2, median, na.rm=T)# #
        mw_QUi[,"Unc"] <- unc_median(mw_srcs_pol_i[,,"Data"],                # #
                                     mw_srcs_pol_i[,,"Unc"], 2)              # #
      }                                                                      # #
                                                                             # #
      temp <- P_from_QU(mw_QUi['Q', "Data"], mw_QUi['U', "Data"],            # #
                        mw_QUi['Q', "Unc"], mw_QUi['U', "Unc"])              # #
      mw_PXi['P', "Data"] <- temp$Data                                       # #
      mw_PXi['P', "Unc"] <- temp$Unc                                         # #
      temp <- X_from_QU(mw_QUi['Q', "Data"], mw_QUi['U', "Data"],            # #
                        mw_QUi['Q', "Unc"], mw_QUi['U', "Unc"])              # #
      mw_PXi['X', "Data"] <- temp$Data                                       # #
      mw_PXi['X', "Unc"] <- temp$Unc                                         # #
      rm(temp)                                                               # #
                                                                             # #
      # Output CSV with MW pol parameters with sky                           # #
      # and instrumental spatial corrections                                 # #
      write.table(mw_QUi, mwqui_out_csv, sep = ";", dec = ".")               # #
      write.table(mw_PXi, mwpxi_out_csv, sep = ";", dec = ".")               # #
      rm(mw_PXi, mwpxi_out_csv)                                              # #
                                                                             # #
      # Plot Q vs U with sky and instrumental spacial corrections            # #
      # of MW sources                                                        # #
      print(paste0("* Creating PDF file with instrumentaly corrected Q vs U",# #
                   " plot of selected Milky Way sources for set #", n, " of",# #
                   " band ", bs, " *"))                                      # #
                                                                             # #
      qlim <- c(min(c(mw_srcs_pol_i[,'Q',"Data"] - mw_srcs_pol_i[,'Q',"Unc"],# #
                      mw_QUi['Q',"Data"] - mw_QUi['Q',"Unc"]),               # #
                    na.rm = T) * 100,                                        # #
                max(c(mw_srcs_pol_i[,'Q',"Data"] + mw_srcs_pol_i[,'Q',"Unc"],# #
                      mw_QUi['Q',"Data"] + mw_QUi['Q',"Unc"]),               # #
                    na.rm = T) * 100)                                        # #
      ulim <- c(min(c(mw_srcs_pol_i[,'U',"Data"] - mw_srcs_pol_i[,'U',"Unc"],# #
                      mw_QUi['U',"Data"] - mw_QUi['U',"Unc"]),               # #
                    na.rm = T) * 100,                                        # #
                max(c(mw_srcs_pol_i[,'U',"Data"] + mw_srcs_pol_i[,'U',"Unc"],# #
                      mw_QUi['U',"Data"] + mw_QUi['U',"Unc"]),               # #
                    na.rm = T) * 100)                                        # #
                                                                             # #
      pdf(mwqui_out_pdf, width = 15, height = 15)                            # #
      rm(mwqui_out_pdf)                                                      # #
                                                                             # #
      plot(mw_srcs_pol_i[,'Q',"Data"]*100, mw_srcs_pol_i[,'U',"Data"]*100,   # #
           xlab = 'q (%)', ylab = 'u (%)', xlim = qlim, ylim = ulim, pch =20,# #
           cex = 2, col = 'black',                                           # #
           main = "instC Q vs U of selected Milky Way sources")              # #
      points(mw_QUi['Q', "Data"] * 100, mw_QUi['U', "Data"] * 100, pch = 14, # #
             cex = 3, col = "red")                                           # #
      arrows(mw_srcs_pol_i[,'Q', "Data"] * 100,                              # #
             y0 = (mw_srcs_pol_i[,'U',"Data"] -                              # #
                     mw_srcs_pol_i[,'U',"Unc"]) * 100,                       # #
             y1 = (mw_srcs_pol_i[,'U',"Data"] +                              # #
                     mw_srcs_pol_i[,'U',"Unc"]) * 100,                       # #
             code = 3, length = 0.02, angle = 90)                            # #
      arrows(x0 = (mw_srcs_pol_i[,'Q',"Data"] -                              # #
                     mw_srcs_pol_i[,'Q',"Unc"]) * 100,                       # #
             mw_srcs_pol_i[,'U', "Data"] * 100,                              # #
             x1 = (mw_srcs_pol_i[,'Q',"Data"] +                              # #
                     mw_srcs_pol_i[,'Q',"Unc"]) * 100,                       # #
             code = 3, length = 0.02, angle = 90)                            # #
      arrows(mw_QUi['Q',"Data"] * 100,                                       # #
             y0 = (mw_QUi['U',"Data"] - mw_QUi['U',"Unc"]) * 100,            # #
             y1 = (mw_QUi['U',"Data"] + mw_QUi['U',"Unc"]) * 100,            # #
             code = 3, length = 0.03, angle = 90, col = "red")               # #
      arrows(x0 = (mw_QUi['Q',"Data"] - mw_QUi['Q',"Unc"]) * 100,            # #
             mw_QUi['U',"Data"] * 100,                                       # #
             x1 = (mw_QUi['Q',"Data"] + mw_QUi['Q',"Unc"]) * 100,            # #
             code = 3, length = 0.03, angle = 90, col = "red")               # #
      dev.off()                                                              # #
      rm(qlim, ulim)                                                         # #
                                                                             # #
      # Applying instrumental chromatic correction to MW sources Stokes Q & U# #
      print(paste0("* Applying chromatic correction to Stokes Q and U of se",# #
                   "lected Milky Way sources of set #", n, " of band ", bs,  # #
                   " *"))                                                    # #
                                                                             # #
      for(ns in 1:n_srcs){                                                   # #
        nsQ <- mw_srcs_pol_i[ns,'Q',]                                        # #
        nsU <- mw_srcs_pol_i[ns,'U',]                                        # #
        nsQc2X <- nsQ["Data"] * c2X                                          # #
        nsQs2X <- nsQ["Data"] * s2X                                          # #
        nsUc2X <- nsU["Data"] * c2X                                          # #
        nsUs2X <- nsU["Data"] * s2X                                          # #
        u_nsQc2X <- unc_mult(nsQ["Data"], c2X, nsQ["Unc"], u_c2X)            # #
        u_nsQs2X <- unc_mult(nsQ["Data"], s2X, nsQ["Unc"], u_s2X)            # #
        u_nsUc2X <- unc_mult(nsU["Data"], c2X, nsU["Unc"], u_c2X)            # #
        u_nsUs2X <- unc_mult(nsU["Data"], s2X, nsU["Unc"], u_s2X)            # #
                                                                             # #
        mw_srcs_pol_cc[ns,'Q',"Data"] <- nsQc2X + nsUs2X                     # #
        mw_srcs_pol_cc[ns,'U',"Data"] <- nsUc2X - nsQs2X                     # #
        mw_srcs_pol_cc[ns,'Q',"Unc"] <- unc_add(u_nsQc2X, u_nsUs2X)          # #
        mw_srcs_pol_cc[ns,'U',"Unc"] <- unc_add(u_nsUc2X, u_nsQs2X)          # #
                                                                             # #
        mw_srcs_QUPX[ns, c("Q","uQ")] <- c(mw_srcs_pol_cc[ns,'Q',"Data"],    # #
                                           mw_srcs_pol_cc[ns,'Q',"Unc"])     # #
        mw_srcs_QUPX[ns, c("U","uU")] <- c(mw_srcs_pol_cc[ns,'U',"Data"],    # #
                                           mw_srcs_pol_cc[ns,'U',"Unc"])     # #
                                                                             # #
        temp <- debiased_P_from_QU(mw_srcs_pol_cc[ns,'Q',"Data"],            # #
                                   mw_srcs_pol_cc[ns,'U',"Data"],            # #
                                   mw_srcs_pol_cc[ns,'Q',"Unc"],             # #
                                   mw_srcs_pol_cc[ns,'U',"Unc"])             # #
        mw_srcs_QUPX[ns, c("P","uP")] <- c(temp$Data, temp$Unc)              # #
                                                                             # #
        temp <- X_from_QU(mw_srcs_pol_cc[ns,'Q',"Data"],                     # #
                          mw_srcs_pol_cc[ns,'U',"Data"],                     # #
                          mw_srcs_pol_cc[ns,'Q',"Unc"],                      # #
                          mw_srcs_pol_cc[ns,'U',"Unc"])                      # #
        mw_srcs_QUPX[ns, c("X","uX")] <- c(temp$Data, temp$Unc)              # #
                                                                             # #
      }                                                                      # #
      rm(ns, mw_srcs_pol_i, nsQ, nsU, nsQc2X, nsQs2X, nsUc2X, nsUs2X,        # #
         u_nsQc2X, u_nsQs2X, u_nsUc2X, u_nsUs2X, n_srcs)                     # #
                                                                             # #
      # Output CSV with MW sources pol parameters with all corrections       # #
      write.table(mw_srcs_pol_cc, src_polc_out_csv, sep = ";", dec =".")     # #
      write.table(mw_srcs_QUPX, mwqupx_out_csv, sep = ";", dec = ".")        # #
      rm(src_polc_out_csv, mw_srcs_QUPX)                                     # #
                                                                             # #
      # Estimating Stokes Q and U from MW with sky                           # #
      # and all instrumental corrections                                     # #
      print(paste0("* Calculating median Stokes Q and U of selected Milky W",# #
                   "ay sources with chromatic correction, for an estimate o",# #
                   "f the foreground polarization level of set #", n, " of ",# #
                   "band ", bs, " *"))                                       # #
                                                                             # #
      if(is.null(dim(mw_srcs_pol_cc[,,"Data"]))){                            # #
        mw_QUcc[,"Data"] <- mw_srcs_pol_cc[,,"Data"]                         # #
        mw_QUcc[,"Unc"] <- mw_srcs_pol_cc[,,"Unc"]                           # #
      }else{                                                                 # #
        mw_QUcc[,"Data"] <- apply(mw_srcs_pol_cc[,,"Data"], 2,median,na.rm=T)# #
        mw_QUcc[,"Unc"] <- unc_median(mw_srcs_pol_cc[,,"Data"],              # #
                                      mw_srcs_pol_cc[,,"Unc"], 2)            # #
      }                                                                      # #
                                                                             # #
      temp <- P_from_QU(mw_QUcc['Q', "Data"], mw_QUcc['U', "Data"],          # #
                        mw_QUcc['Q', "Unc"], mw_QUcc['U', "Unc"])            # #
      mw_PXcc['P', "Data"] <- temp$Data                                      # #
      mw_PXcc['P', "Unc"] <- temp$Unc                                        # #
      temp <- X_from_QU(mw_QUcc['Q', "Data"], mw_QUcc['U', "Data"],          # #
                        mw_QUcc['Q', "Unc"], mw_QUcc['U', "Unc"])            # #
      mw_PXcc['X', "Data"] <- temp$Data                                      # #
      mw_PXcc['X', "Unc"] <- temp$Unc                                        # #
      rm(temp)                                                               # #
                                                                             # #
      # Output CSV with MW pol parameters with all corrections               # #
      write.table(mw_QUcc, mwqucc_out_csv, sep = ";", dec = ".")             # #
      write.table(mw_PXcc, mwpxcc_out_csv, sep = ";", dec = ".")             # #
      rm(mwqucc_out_csv,  mwpxcc_out_csv)                                    # #
                                                                             # #
      # Plot Q vs U with all corrections for all MW sources                  # #
      print(paste0("* Creating PDF file with chromaticaly corrected Q vs U ",# #
                   "plot of selected Milky Way sources for set #", n, " of ",# #
                   "band ", bs, " *"))                                       # #
                                                                             # #
                                                                             # #
      pdf(mwqucc_out_pdf, width = 15, height = 15)                           # #
      rm(mwqucc_out_pdf)                                                     # #
                                                                             # #
      qlim <- c(min(c(mw_srcs_pol_cc[,'Q', "Data"] -                         # #
                        mw_srcs_pol_cc[,'Q', "Unc"],                         # #
                      mw_QUcc['Q', "Data"] - mw_QUcc['Q', "Unc"]),           # #
                    na.rm = T) * 100,                                        # #
                max(c(mw_srcs_pol_cc[,'Q', "Data"] +                         # #
                        mw_srcs_pol_cc[,'Q', "Unc"],                         # #
                      mw_QUcc['Q', "Data"] + mw_QUcc['Q', "Unc"]),           # #
                    na.rm = T) * 100)                                        # #
      ulim <- c(min(c(mw_srcs_pol_cc[,'U', "Data"] -                         # #
                        mw_srcs_pol_cc[,'U', "Unc"],                         # #
                      mw_QUcc['U', "Data"] - mw_QUcc['U', "Unc"]),           # #
                    na.rm = T) * 100,                                        # #
                max(c(mw_srcs_pol_cc[,'U', "Data"] +                         # #
                        mw_srcs_pol_cc[,'U', "Unc"],                         # #
                      mw_QUcc['U', "Data"] + mw_QUcc['U', "Unc"]),           # #
                    na.rm = T) * 100)                                        # #
                                                                             # #
      plot(mw_srcs_pol_cc[,'Q',"Data"]*100, mw_srcs_pol_cc[,'U',"Data"]*100, # #
           xlab = 'q (%)', ylab = 'u (%)', xlim = qlim, ylim = ulim,         # #
           pch = 20, cex = 2, col = 'black',                                 # #
           main = "CC Q vs U of selected Milky Way sources")                 # #
      points(mw_QUcc['Q', "Data"] * 100, mw_QUcc['U', "Data"] * 100,         # #
             pch = 14, cex = 3, col = "red")                                 # #
      arrows(mw_srcs_pol_cc[,'Q', "Data"] * 100,                             # #
             y0 = (mw_srcs_pol_cc[,'U', "Data"] - mw_srcs_pol_cc[,'U',"Unc"])# #
             * 100,                                                          # #
             y1 = (mw_srcs_pol_cc[,'U', "Data"] + mw_srcs_pol_cc[,'U',"Unc"])# #
             * 100, code = 3, length = 0.02, angle = 90)                     # #
      arrows(x0 = (mw_srcs_pol_cc[,'Q', "Data"] - mw_srcs_pol_cc[,'Q',"Unc"])# #
             * 100,                                                          # #
             mw_srcs_pol_cc[,'U', "Data"] * 100,                             # #
             x1 = (mw_srcs_pol_cc[,'Q', "Data"] + mw_srcs_pol_cc[,'Q',"Unc"])# #
             * 100, code = 3, length = 0.02, angle = 90)                     # #
      arrows(mw_QUcc['Q',"Data"] * 100,                                      # #
             y0 = (mw_QUcc['U',"Data"] - mw_QUcc['U',"Unc"]) * 100,          # #
             y1 = (mw_QUcc['U',"Data"] + mw_QUcc['U',"Unc"]) * 100,          # #
             code = 3, length = 0.03, angle = 90, col = "red")               # #
      arrows(x0 = (mw_QUcc['Q',"Data"] - mw_QUcc['Q',"Unc"]) * 100,          # #
             mw_QUcc['U',"Data"] * 100,                                      # #
             x1 = (mw_QUcc['Q',"Data"] + mw_QUcc['Q',"Unc"]) * 100,          # #
             code = 3, length = 0.03, angle = 90, col = "red")               # #
      dev.off()                                                              # #
      rm(qlim, ulim)                                                         # #
                                                                             # #
      # Debiasing Milk Way's P estimation                                    # #
      print(paste0("* Applying debiasing of MW estimation of P for set #", n,# #
                   " of band ", bs, " *"))                                   # #
                                                                             # #
      temp <- debiased_P_from_QU(mw_QUcc['Q', "Data"], mw_QUcc['U', "Data"], # #
                                 mw_QUcc['Q', "Unc"], mw_QUcc['U', "Unc"])   # #
      mw_PXdeb['P', "Data"] <- temp$Data                                     # #
      mw_PXdeb['P', "Unc"] <- temp$Unc                                       # #
      mw_PXdeb['X', "Data"] <- mw_PXcc['X', "Data"]                          # #
      mw_PXdeb['X', "Unc"] <- mw_PXcc['X', "Unc"]                            # #
      rm(temp, mw_QUcc, mw_PXcc)                                             # #
                                                                             # #
      # Output CSV with MW P and X all corrections, and P debias             # #
      write.table(mw_PXdeb, mwpxdeb_out_csv, sep = ";", dec = ".")           # #
      rm(mw_PXdeb)                                                           # #
      ##########################################################################
    }else{                                                                     #
      ################ Loading Estimate of MW Stokes Parameters ############## #
      ############# with sky and instrumental spatial corrections ############ #
      ##########################################################################
      mw_QUi <- read.table(mwqui_out_csv, header = T, sep = ";", dec = ".")  # #
      ##########################################################################
    }                                                                          #
    rm(mwqui_out_csv, mwqupx_out_csv, exist_flag)                              #
                                                                               #
    E_nosky_out_fits <- paste0(trgt_noSky_folder, fits_fname,                  #
                               "-Es_noSky_noCut.fits")                         #
    O_nosky_out_fits <- paste0(trgt_noSky_folder, fits_fname,                  #
                               "-Os_noSky_noCut.fits")                         #
    I_gal_out_fits <- paste0(noSky_Iflux_folder, fits_fname,                   #
                             "-I_noSky_noCut.fits")                            #
                                                                               #
    exist_flag <- prod(file.exists(c(E_nosky_out_fits, O_nosky_out_fits,       #
                                     I_gal_out_fits)))                         #
                                                                               #
    Es_noSky <- Es                                                             #
    Os_noSky <- Os                                                             #
    Igal <- array(NA, dim = singDim, dimnames = singNames)                     #
                                                                               #
    if(!exist_flag){                                                           #
      ######################### Removing Sky from Beams ###################### #
      ##########################################################################
      print(paste0("* Removing Sky from beams for set #", n, " of band ", bs,# #
                   " *"))                                                    # #
                                                                             # #
      Es_noSky[,,"Data",] <- Es[,,"Data",] - Es_sky[,,"Data",]               # #
      Es_noSky[,,"Unc",] <- unc_add(Es[,,"Unc",], Es_sky[,,"Unc",])          # #
      Os_noSky[,,"Data",] <- Os[,,"Data",] - Os_sky[,,"Data",]               # #
      Os_noSky[,,"Unc",] <- unc_add(Os[,,"Unc",], Os_sky[,,"Unc",])          # #
                                                                             # #
      print(paste0("* Saving sky subtracted beams for set #", n, " of band ",# #
                   bs, " *"))                                                # #
                                                                             # #
      # Saving No Cut Sky corrected beams                                    # #
      writeFITSim(Es_noSky, file = E_nosky_out_fits,                         # #
                  crvaln = crval2, crpixn = crpix2, ctypen = ctype2)         # #
      writeFITSim(Os_noSky, file = O_nosky_out_fits,                         # #
                  crvaln = crval2, crpixn = crpix2, ctypen = ctype2)         # #
      ##########################################################################
                                                                               #
      #----------------------------------------------------------------------# #
                                                                               #
      ################### Estimating Sky Subtracted Flux Map ################# #
      ##########################################################################
      Is_noSky <- Es_noSky[,,"Data",] + Os_noSky[,,"Data",]                  # #
      Is_noSky_unc <- unc_add(Es_noSky[,,"Unc",], Os_noSky[,,"Unc",])        # #
                                                                             # #
      print(paste0("* Creating median sky subtracted intensity map for set ",# #
                   "#", n, " of band ", bs, " *"))                           # #
                                                                             # #
      temp <- apply(Is_noSky, 1:2, median, na.rm = T)                        # #
      temp[NAs_ind] <- NA                                                    # #
                                                                             # #
      Igal[,,"Data"] <- temp                                                 # #
                                                                             # #
      temp <- unc_median(Is_noSky, Is_noSky_unc, 1:2)                        # #
      temp[NAs_ind] <- NA                                                    # #
                                                                             # #
      Igal[,,"Unc"] <- temp                                                  # #
      rm(Is_noSky, Is_noSky_unc, temp, t_nas, lt_nas, m_temp, s_temp)        # #
                                                                             # #
      writeFITSim(Igal, file = I_gal_out_fits, crvaln = crval,               # #
                  crpixn = crpix, ctypen = ctype, header = I_head)           # #
      ##########################################################################
    }else{                                                                     #
      ###################### Loading Sky Subtracted Beams #################### #
      ##########################################################################
      print(paste0("* Loading sky subtracted beams for set #",n, " of band ",# #
                   bs, " *"))                                                # #
                                                                             # #
      Es_noSky[,,,] <- readFITS(E_nosky_out_fits)$imDat                      # #
      Os_noSky[,,,] <- readFITS(O_nosky_out_fits)$imDat                      # #
      ##########################################################################
                                                                               #
      #----------------------------------------------------------------------# #
                                                                               #
      #################### Loading Sky Subtracted Flux Map ################### #
      ##########################################################################
      print(paste0("* Loading median sky subtracted intensity map for set #",# #
                   n, " of band ", bs, " *"))                                # #
                                                                             # #
      Igal[,,] <- readFITS(I_gal_out_fits)$imDat                             # #
      ##########################################################################
    }                                                                          #
    rm(E_nosky_out_fits, O_nosky_out_fits, I_gal_out_fits, exist_flag)         #
                                                                               #
    Ehc_out_fits <- paste0(trgt_noSky_folder, fits_fname, "-Es_noSky_Cut.fits")#
    Ohc_out_fits <- paste0(trgt_noSky_folder, fits_fname, "-Os_noSky_Cut.fits")#
    Ihc_out_fits <- paste0(noSky_Iflux_folder, fits_fname, "-I_noSky_Cut.fits")#
    pdf_hist <- paste0(trgt_noSky_folder, fits_fname,"-beams_Histograms.pdf")  #
                                                                               #
    exist_flag <- prod(file.exists(c(Ehc_out_fits, Ohc_out_fits, Ihc_out_fits))#
                       )                                                       #
                                                                               #
    if(!exist_flag){                                                           #
      #### Creating Sky Corrected Beams and Flux Maps with Different Cuts #### #
      ##########################################################################
      print(paste0("* Creating sky subtracted beams and intensity maps with",# #
                   " different statistical cuts for set #", n, " of band ",  # #
                   bs, " *"))                                                # #
                                                                             # #
      pdf(pdf_hist, width = 80, height = 55)                                 # #
      rm(pdf_hist)                                                           # #
      par(mfrow = c(1, 2))                                                   # #
      par(mar = c(15, 15, 10, 5))                                            # #
      par(cex.main = 7)                                                      # #
      par(cex.axis = 6)                                                      # #
      par(cex.lab = 6)                                                       # #
                                                                             # #
      # Deriving adequate statistics to determine a reasonable flux level    # #
      # consistent with sky, create a plot of those measures, then remove it # #
      # from beams                                                           # #
      beam_cuts <- array(NA, dim = c(2, length(angs)),                       # #
                         dimnames = list(c("O", "E"), ang_str))              # #
      for(ang in angs){                                                      # #
        as <- ang_str[ang]                                                   # #
                                                                             # #
        Odat <- Os_noSky[,,"Data", ang]                                      # #
        Edat <- Es_noSky[,,"Data", ang]                                      # #
        Odat <- Odat[which(!is.na(Odat))]                                    # #
        Edat <- Edat[which(!is.na(Edat))]                                    # #
                                                                             # #
        Osky <- Os_sky[,,"Data", ang]                                        # #
        Esky <- Es_sky[,,"Data", ang]                                        # #
        Osky <- Osky[which(!is.na(Osky))]                                    # #
        Esky <- Esky[which(!is.na(Esky))]                                    # #
                                                                             # #
        Odat <- Odat[which(Odat < (5 * sd(Osky)))]                           # #
        Edat <- Edat[which(Edat < (5 * sd(Esky)))]                           # #
        ON <- length(Odat)                                                   # #
        EN <- length(Edat)                                                   # #
                                                                             # #
        # Defining sigma-clip like intervals                                 # #
        temp_Odat <- Odat[which(Odat <= 0, arr.ind = T)]                     # #
        mirr_Odat <- c(-temp_Odat, temp_Odat)                                # #
        mO <- median(mirr_Odat, na.rm = T)                                   # #
        sdO <- mad(mirr_Odat, na.rm = T)                                     # #
        sclO <- length(mirr_Odat)                                            # #
                                                                             # #
        # Defining sigma-clip like intervals                                 # #
        temp_Edat <- Edat[which(Edat <= 0, arr.ind = T)]                     # #
        mirr_Edat <- c(-temp_Edat, temp_Edat)                                # #
        mE <- median(mirr_Edat, na.rm = T)                                   # #
        sdE <- mad(mirr_Edat, na.rm = T)                                     # #
        sclE <- length(mirr_Edat)                                            # #
        rm(mirr_Odat, mirr_Edat, temp_Odat, temp_Edat)                       # #
                                                                             # #
        Zhist_len <- 200                                                     # #
                                                                             # #
        # Setting up breaks for zoom plot                                    # #
        sd5O <- 5 * sdO                                                      # #
        sd5E <- 5 * sdE                                                      # #
        Obrks_zoom <- seq(mO - sd5O, mO + sd5O, length.out = Zhist_len)      # #
        Ebrks_zoom <- seq(mE - sd5E, mE + sd5E, length.out = Zhist_len)      # #
                                                                             # #
        cutO <- mO                                                           # #
        tol_flag <- FALSE                                                    # #
                                                                             # #
        while(!tol_flag){                                                    # #
          # Finding interval where deviation sdO_cut*sd from the mean is     # #
          brks_s <- Obrks_zoom[which(Obrks_zoom < cutO, arr.ind = T)]        # #
          brks_b <- Obrks_zoom[which(Obrks_zoom > cutO, arr.ind = T)]        # #
                                                                             # #
          # If the search has gone beyond the available range                # #
          if(length(brks_b) == 0){                                           # #
            tol_flag <- TRUE                                                 # #
            break                                                            # #
          }                                                                  # #
                                                                             # #
          x1 <- brks_s[length(brks_s)]                                       # #
          x2 <- brks_b[1]                                                    # #
          del_x <- x2 - x1                                                   # #
                                                                             # #
          # Estimating the amount of pixels within the interval              # #
          p1 <- pnorm(x1, mO, sdO)                                           # #
          p2 <- pnorm(x2, mO, sdO)                                           # #
          del_p <- p2 - p1                                                   # #
          est_N <- del_p * sclO                                              # #
                                                                             # #
          # Counting amount of pixels in that interval                       # #
          brks_N <- length(Odat[which((Odat >= x1) & (Odat <= x2))])         # #
                                                                             # #
          # Tolerance probability, if predicted number is too small          # #
          # cutO is not update and the search cycle breaks                   # #
          tol_N <- q_prob * brks_N                                           # #
                                                                             # #
          if(est_N < tol_N){                                                 # #
            tol_flag <- TRUE                                                 # #
            break                                                            # #
          }                                                                  # #
          cutO <- cutO + del_x / 2                                           # #
        }                                                                    # #
                                                                             # #
        cutE <- mE                                                           # #
        tol_flag <- FALSE                                                    # #
                                                                             # #
        while(!tol_flag){                                                    # #
          brks_s <- Ebrks_zoom[which(Ebrks_zoom < cutE, arr.ind = T)]        # #
          brks_b <- Ebrks_zoom[which(Ebrks_zoom > cutE, arr.ind = T)]        # #
                                                                             # #
          if(length(brks_b) == 0){                                           # #
            tol_flag <- TRUE                                                 # #
            break                                                            # #
          }                                                                  # #
                                                                             # #
          x1 <- brks_s[length(brks_s)]                                       # #
          x2 <- brks_b[1]                                                    # #
          del_x <- x2 - x1                                                   # #
                                                                             # #
          # Estimating the amount of pixels within the interval              # #
          p1 <- pnorm(x1, mE, sdE)                                           # #
          p2 <- pnorm(x2, mE, sdE)                                           # #
          del_p <- p2 - p1                                                   # #
          est_N <- del_p * sclE                                              # #
                                                                             # #
          # Counting amount of pixels in that interval                       # #
          brks_N <- length(Edat[which((Edat >= x1) & (Edat <= x2))])         # #
                                                                             # #
          # Tolerance probability, if predicted number is too small          # #
          # cutO is not update and the search cycle breaks                   # #
          tol_N <- q_prob * brks_N                                           # #
                                                                             # #
          if(est_N < tol_N){                                                 # #
            tol_flag <- TRUE                                                 # #
            break                                                            # #
          }                                                                  # #
          cutE <- cutE + del_x / 2                                           # #
        }                                                                    # #
        rm(x1, x2, p1, p2, del_x, del_p, est_N, tol_N, tol_flag, ON, EN,     # #
           brks_s, brks_b, brks_N)                                           # #
                                                                             # #
        # Trimming distribution tales for zoom plot                          # #
        temp_Odat <- Odat                                                    # #
        temp_Edat <- Edat                                                    # #
        temp_Odat[which(Odat < Obrks_zoom[1], arr.ind = T)] <- NA            # #
        temp_Odat[which(Odat > Obrks_zoom[Zhist_len], arr.ind = T)] <- NA    # #
        Odat <- temp_Odat                                                    # #
        temp_Edat[which(Edat < Ebrks_zoom[1], arr.ind = T)] <- NA            # #
        temp_Edat[which(Edat > Ebrks_zoom[Zhist_len], arr.ind = T)] <- NA    # #
        Edat <- temp_Edat                                                    # #
        rm(temp_Odat, temp_Edat)                                             # #
                                                                             # #
        Ohist_title <- paste0("Histogram of O beam w/ HWP at ", as, "º")     # #
        Ehist_title <- paste0("Histogram of E beam w/ HWP at ", as, "º")     # #
                                                                             # #
        Oh <- hist(Odat, breaks = Obrks_zoom, freq = T, main = Ohist_title,  # #
                   xlab = "", ylab = "", axes = F)                           # #
                                                                             # #
        ################# Determining limits of plot axis ################## # #
        ########################################################################
        max_counts <- Oh$counts[which.max(Oh$counts)]                      # # #
        maxPow <- floor(log10(max_counts))                                 # # #
        ymax <- round(max_counts / (10^maxPow)) * 10^maxPow                # # #
        rm(max_counts)                                                     # # #
                                                                           # # #
        yticks <- seq(0, ymax, length.out = 5)                             # # #
        ylabels <- c("0")                                                  # # #
        for(i in 2:length(yticks)){                                        # # #
          powi <- floor(log10(yticks[i]))                                  # # #
          yi <- round(yticks[i] / (10^powi), 1)                            # # #
          ylabels <- c(ylabels, paste0(yi, " E+", powi))                   # # #
        }                                                                  # # #
                                                                           # # #
        max_brk <- (Oh$breaks[Zhist_len] + Oh$breaks[Zhist_len - 1]) / 2   # # #
        maxPow <- floor(log10(max_brk))                                    # # #
        xmax <- round(max_brk / (10^maxPow)) * 10^maxPow                   # # #
        rm(max_brk)                                                        # # #
                                                                           # # #
        min_brk <- (Oh$breaks[1] + Oh$breaks[2]) / 2                       # # #
        if(min_brk == 0){                                                  # # #
          xmin <- 0                                                        # # #
        }else{                                                             # # #
          maxPow <- floor(log10(abs(min_brk)))                             # # #
          xmin <- round(min_brk / (10^maxPow)) * 10^maxPow                 # # #
        }                                                                  # # #
        rm(min_brk, maxPow)                                                # # #
                                                                           # # #
        if(xmin == 0){                                                     # # #
          xticks <- seq(xmin, xmax, length.out = 5)                        # # #
        }else{                                                             # # #
          delx <- (xmax - xmin) / 5                                        # # #
          xticks <- seq(0, xmin, by = -delx)                               # # #
          xticks <- sort(c(xticks[-1], seq(0, xmax, by = delx)))           # # #
          rm(delx)                                                         # # #
        }                                                                  # # #
        xlabels <- paste0("", round(xticks))                               # # #
        ########################################################################
                                                                             # #
        # Plotting Axes                                                      # #
        axis(1, at = xticks, labels = xlabels, line = 3, lwd = 0)            # #
        axis(1, at = xticks, labels = rep("", length(xticks)), lwd.ticks = 2)# #
        axis(2, at = yticks, labels = ylabels, line = 1, lwd = 0)            # #
        axis(2, at = yticks, labels = rep("", length(yticks)), lwd.ticks = 2)# #
                                                                             # #
        curve(dnorm(x, mO, sdO) * sclO * diff(Oh$breaks)[1], add = T,        # #
              col = "red")                                                   # #
        segments(x0 = cutO, y0 = 0, x1 = cutO, y1 = ymax,                    # #
                 col = "blue", lty = 3, lwd = 2)                             # #
        segments(x0 = mO, y0 = 0, x1 = mO, y1 = ymax,                        # #
                 col = "black", lty = 1, lwd = 2)                            # #
        mtext("Flux (arbitrary units)", 1, line = 12, cex = 7)               # #
        mtext("Frequency", 2, line = 10, cex = 7)                            # #
        rm(mO)                                                               # #
                                                                             # #
        Eh <- hist(Edat, breaks = Ebrks_zoom, freq = T, main = Ehist_title,  # #
                   xlab = "", ylab = "", axes = FALSE)                       # #
                                                                             # #
        ################# Determining limits of plot axis ################## # #
        ########################################################################
        max_counts <- Eh$counts[which.max(Eh$counts)]                      # # #
        maxPow <- floor(log10(max_counts))                                 # # #
        ymax <- round(max_counts / (10^maxPow)) * 10^maxPow                # # #
        rm(max_counts)                                                     # # #
                                                                           # # #
        yticks <- seq(0, ymax, length.out = 5)                             # # #
        ylabels <- c("0")                                                  # # #
        for(i in 2:length(yticks)){                                        # # #
          powi <- floor(log10(yticks[i]))                                  # # #
          yi <- round(yticks[i] / (10^powi), 1)                            # # #
          ylabels <- c(ylabels, paste0(yi, " E+", powi))                   # # #
        }                                                                  # # #
                                                                           # # #
        max_brk <- (Eh$breaks[Zhist_len] + Eh$breaks[Zhist_len - 1]) / 2   # # #
        maxPow <- floor(log10(max_brk))                                    # # #
        xmax <- round(max_brk / (10^maxPow)) * 10^maxPow                   # # #
        rm(max_brk)                                                        # # #
                                                                           # # #
        min_brk <- (Eh$breaks[1] + Eh$breaks[2]) / 2                       # # #
        if(min_brk == 0){                                                  # # #
          xmin <- 0                                                        # # #
        }else{                                                             # # #
          maxPow <- floor(log10(abs(min_brk)))                             # # #
          xmin <- round(min_brk / (10^maxPow)) * 10^maxPow                 # # #
        }                                                                  # # #
        rm(min_brk, maxPow)                                                # # #
                                                                           # # #
        if(xmin == 0){                                                     # # #
          xticks <- seq(xmin, xmax, length.out = 5)                        # # #
        }else{                                                             # # #
          delx <- (xmax - xmin) / 5                                        # # #
          xticks <- seq(0, xmin, by = -delx)                               # # #
          xticks <- sort(c(xticks[-1], seq(0, xmax, by = delx)))           # # #
          rm(delx)                                                         # # #
        }                                                                  # # #
        xlabels <- paste0("", round(xticks))                               # # #
        ########################################################################
                                                                             # #
        # Plotting Axes                                                      # #
        axis(1, at = xticks, labels = xlabels, line = 3, lwd = 0)            # #
        axis(1, at = xticks, labels = rep("", length(xticks)), lwd.ticks = 2)# #
        axis(2, at = yticks, labels = ylabels, line = 1, lwd = 0)            # #
        axis(2, at = yticks, labels = rep("", length(yticks)), lwd.ticks = 2)# #
                                                                             # #
        curve(dnorm(x, mE, sdE) * sclE * diff(Eh$breaks)[1], add = T,        # #
              col = "red")                                                   # #
        segments(x0 = cutE, y0 = 0, x1 = cutE, y1 = ymax,                    # #
                 col = "blue", lty = 3, lwd = 2)                             # #
        segments(x0 = mE, y0 = 0, x1 = mE, y1 = ymax,                        # #
                 col = "black", lty = 1, lwd = 2)                            # #
        mtext("Flux (arbitrary units)", 1, line = 12, cex = 7)               # #
        mtext("Frequency", 2, line = 10, cex = 7)                            # #
        rm(mE, Ohist_title, Ehist_title, Odat, Edat, Obrks_zoom, Ebrks_zoom, # #
           Zhist_len)                                                        # #
                                                                             # #
        beam_cuts["O", ang] <- cutO                                          # #
        beam_cuts["E", ang] <- cutE                                          # #
      }                                                                      # #
      dev.off()                                                              # #
      rm(Oh, Eh)                                                             # #
                                                                             # #
      # SKY FLUX CUTS                                                        # #
      for(ang in angs){                                                      # #
        # Sky cut                                                            # #
        Es_noSky[,,"Data",ang] <- replace_qtts(Es_noSky[,,"Data", ang],      # #
                                               beam_cuts["E", ang], '<=', NA)# #
        Os_noSky[,,"Data",ang] <- replace_qtts(Os_noSky[,,"Data", ang],      # #
                                               beam_cuts["O", ang], '<=', NA)# #
                                                                             # #
        Es_noSky[,, "Unc",ang] <- bij_match_val(Es_noSky[,,"Data", ang],     # #
                                                Es_noSky[,,"Unc", ang],      # #
                                                NA, NA)                      # #
        Os_noSky[,, "Unc",ang] <- bij_match_val(Os_noSky[,,"Data", ang],     # #
                                                Os_noSky[,,"Unc", ang],      # #
                                                NA, NA)                      # #
      }                                                                      # #
      rm(beam_cuts)                                                          # #
                                                                             # #
      # Sum of beams across different angles, sky cut                        # #
      Is_noSky <- Es_noSky[,,"Data",] + Os_noSky[,,"Data",]                  # #
      Is_noSky_unc <- unc_add(Es_noSky[,,"Unc",], Os_noSky[,,"Unc",])        # #
                                                                             # #
      # Median of sum of beams across different angles, sky cut              # #
      temp <- apply(Is_noSky, 1:2, median, na.rm = T)                        # #
      temp[NAs_ind] <- NA                                                    # #
      Igal[,,"Data"] <- temp                                                 # #
                                                                             # #
      # Uncertainty of the median of sum of beams across all angles, sky cut # #
      temp <- unc_median(Is_noSky, Is_noSky_unc, 1:2)                        # #
      temp[NAs_ind] <- NA                                                    # #
      Igal[,,"Unc"] <- temp                                                  # #
      rm(Is_noSky, Is_noSky_unc, temp)                                       # #
                                                                             # #
      # Saving Cut Sky corrected beams                                       # #
      writeFITSim(Es_noSky, file = Ehc_out_fits,                             # #
                  crvaln = crval2, crpixn = crpix2, ctypen = ctype2)         # #
      writeFITSim(Os_noSky, file = Ohc_out_fits,                             # #
                  crvaln = crval2, crpixn = crpix2, ctypen = ctype2)         # #
      writeFITSim(Igal, file = Ihc_out_fits, crvaln = crval, crpixn = crpix, # #
                  ctypen = ctype, header = I_head)                           # #
      ##########################################################################
    }else{                                                                     #
      ##### Loading Sky Corrected Beams and Flux Maps with Different Cuts #### #
      ##########################################################################
      print(paste0("* Loading sky subtracted beams and intensity maps with ",# #
                   "different statistical cuts for set #", n, " of band ",   # #
                   bs, " *"))                                                # #
                                                                             # #
      Igal[,,] <- readFITS(Ihc_out_fits)$imDat                               # #
      Es_noSky[,,,] <- readFITS(Ehc_out_fits)$imDat                          # #
      Os_noSky[,,,] <- readFITS(Ohc_out_fits)$imDat                          # #
      ##########################################################################
    }                                                                          #
    rm(Ehc_out_fits, Ohc_out_fits, Ihc_out_fits, exist_flag, I_head,           #
       skyE_head, skyO_head)                                                   #
                                                                               #
    ########################## Binning Spatial Data ########################## #
    ############################################################################
    if(binS > 1){                                                            # #
      temp_head <- modVal('BINNING', toString(binS), "Bin length", temp_head)# #
      sky_head <- modVal('BINNING', toString(binS), "Bin length", sky_head)  # #
                                                                             # #
      Ebhc_out_fits <- paste0(trgt_noSky_bin_folder, fits_fname,             # #
                              "-Bin_Ebeams_noSky.fits")                      # #
      Obhc_out_fits <- paste0(trgt_noSky_bin_folder, fits_fname,             # #
                              "-Bin_Obeams_noSky.fits")                      # #
      Ebsky_out_fits <- paste0(trgt_sky_bin_folder, fits_fname,              # #
                               "-Bin_Sky_Ebeams.fits")                       # #
      Obsky_out_fits <- paste0(trgt_sky_bin_folder, fits_fname,              # #
                               "-Bin_Sky_Obeams.fits")                       # #
      Ebobs_out_fits <- paste0(obs_bin_folder, fits_fname,                   # #
                               "-Bin_ObsEbeams.fits")                        # #
      Obobs_out_fits <- paste0(obs_bin_folder, fits_fname,                   # #
                               "-Bin_ObsObeams.fits")                        # #
                                                                             # #
      exist_flag <- prod(file.exists(c(Ebhc_out_fits, Obhc_out_fits,         # #
                                       Ebsky_out_fits, Obsky_out_fits,       # #
                                       Ebobs_out_fits, Obobs_out_fits)))     # #
                                                                             # #
      if(!exist_flag){                                                       # #
        ########################## Binning Beams ########################### # #
        ########################################################################
        binEs <- array(NA, dim = bin_dims, dimnames = multNames)           # # #
        binOs <- array(NA, dim = bin_dims, dimnames = multNames)           # # #
        binEs_sky <- array(NA, dim = bin_dims, dimnames = multNames)       # # #
        binOs_sky <- array(NA, dim = bin_dims, dimnames = multNames)       # # #
        binEs_noSky <- array(NA, dim = bin_dims, dimnames = multNames)     # # #
        binOs_noSky <- array(NA, dim = bin_dims, dimnames = multNames)     # # #
                                                                           # # #
        print(paste0("* Binning E and O beam maps pixels by their median ",# # #
                     "within ", binS, "x", binS, " pixel windows, for set",# # #
                     " #", n, " of band ", bs, " *"))                      # # #
                                                                           # # #
        for(ang in angs){                                                  # # #
          binEs[,,,ang] <- bin_code_map(Es[,,,ang], binS, ref_crpix)       # # #
          binOs[,,,ang] <- bin_code_map(Os[,,,ang], binS, ref_crpix)       # # #
          binEs_noSky[,,,ang] <- bin_code_map(Es_noSky[,,,ang], binS,      # # #
                                              ref_crpix)                   # # #
          binOs_noSky[,,,ang] <- bin_code_map(Os_noSky[,,,ang], binS,      # # #
                                              ref_crpix)                   # # #
          binEs_sky[,,,ang] <- bin_code_map(Es_sky[,,,ang], binS,          # # #
                                            ref_crpix)                     # # #
          binOs_sky[,,,ang] <- bin_code_map(Os_sky[,,,ang], binS,          # # #
                                            ref_crpix)                     # # #
        }                                                                  # # #
        rm(Es, Os, Es_noSky, Os_noSky, ang, Es_sky, Os_sky)                # # #
                                                                           # # #
        Es <- binEs                                                        # # #
        Os <- binOs                                                        # # #
        Es_noSky <- binEs_noSky                                            # # #
        Os_noSky <- binOs_noSky                                            # # #
        Es_sky <- binEs_sky                                                # # #
        Os_sky <- binOs_sky                                                # # #
        rm(binEs, binOs, binEs_noSky, binOs_noSky, binEs_sky, binOs_sky)   # # #
                                                                           # # #
        writeFITSim(Es_noSky, file = Ebhc_out_fits, crvaln = crval2,       # # #
                    crpixn = crpix2, ctypen = ctype2)                      # # #
        writeFITSim(Os_noSky, file = Obhc_out_fits, crvaln = crval2,       # # #
                    crpixn = crpix2, ctypen = ctype2)                      # # #
        writeFITSim(Es_sky, file = Ebsky_out_fits, crvaln = crval2,        # # #
                    crpixn = crpix2, ctypen = ctype2)                      # # #
        writeFITSim(Os_sky, file = Obsky_out_fits, crvaln = crval2,        # # #
                    crpixn = crpix2, ctypen = ctype2)                      # # #
        writeFITSim(Es, file = Ebobs_out_fits, crvaln = crval2,            # # #
                    crpixn = crpix2, ctypen = ctype2)                      # # #
        writeFITSim(Os, file = Obobs_out_fits, crvaln = crval2,            # # #
                    crpixn = crpix2, ctypen = ctype2)                      # # #
        ########################################################################
      }else{                                                                 # #
        ####################### Loading Binned Beams ####################### # #
        ########################################################################
        Es <- array(NA, dim = bin_dims, dimnames = multNames)              # # #
        Os <- array(NA, dim = bin_dims, dimnames = multNames)              # # #
        Es_sky <- array(NA, dim = bin_dims, dimnames = multNames)          # # #
        Os_sky <- array(NA, dim = bin_dims, dimnames = multNames)          # # #
        Es_noSky <- array(NA, dim = bin_dims, dimnames = multNames)        # # #
        Os_noSky <- array(NA, dim = bin_dims, dimnames = multNames)        # # #
                                                                           # # #
        Es[,,,] <- readFITS(Ebobs_out_fits)$imDat                          # # #
        Os[,,,] <- readFITS(Obobs_out_fits)$imDat                          # # #
        Es_sky[,,,] <- readFITS(Ebsky_out_fits)$imDat                      # # #
        Os_sky[,,,] <- readFITS(Obsky_out_fits)$imDat                      # # #
        Es_noSky[,,,] <- readFITS(Ebhc_out_fits)$imDat                     # # #
        Os_noSky[,,,] <- readFITS(Obhc_out_fits)$imDat                     # # #
        ########################################################################
      }                                                                      # #
      rm(Ebobs_out_fits, Obobs_out_fits, Ebsky_out_fits, Obsky_out_fits,     # #
         Ebhc_out_fits, Obhc_out_fits, exist_flag)                           # #
                                                                             # #
      print(paste0("* Binning Isky and Igal maps pixels by their median wit",# #
                   "hin ", binS, "x", binS, " pixel windows, for set #", n,  # #
                   " of band ", bs, " *"))                                   # #
                                                                             # #
      bIsky <- bin_code_map(Isky, binS, ref_crpix)                           # #
      bIgal <- bin_code_map(Igal, binS, ref_crpix)                           # #
    }                                                                        # #
    ############################################################################
                                                                               #
                                                                               #
    Ibhc_stats_out_csv <- paste0(trgt_noSky_stat_folder, fits_fname,           #
                                 "-Bin_I_Gal_stats.csv")                       #
    Isky_stats_out_csv <- paste0(trgt_sky_stat_folder, fits_fname,             #
                                 "-Bin_I_Sky_stats.csv")                       #
    Bbhc_stats_out_csv <- paste0(trgt_noSky_stat_folder, fits_fname,           #
                                 "-Bin_Beams_stats.csv")                       #
    Bbobs_stats_out_csv <- paste0(obs_stat_folder, fits_fname,                 #
                                  "-Bin_Beams_stats.csv")                      #
    Bbsky_stats_out_csv <- paste0(trgt_sky_stat_folder, fits_fname,            #
                                  "-Bin_Sky_beam_stats.csv")                   #
                                                                               #
    exist_flag <- prod(file.exists(c(Ibhc_stats_out_csv, Isky_stats_out_csv,   #
                                     Bbhc_stats_out_csv, Bbobs_stats_out_csv,  #
                                     Bbsky_stats_out_csv)))                    #
                                                                               #
    stat_par <- c("min", "max", "median", "mad")                               #
    pos_par <- c("val", "x", "y")                                              #
                                                                               #
    if(!exist_flag){                                                           #
      ############# Producing Statistics for Beams and Flux Maps ############# #
      ##########################################################################
      print(paste0("* Producing statistics for beams and intensity maps wit",# #
                   "h different statistical cuts for set #", n, " of band ", # #
                   bs, " *"))                                                # #
                                                                             # #
      beam_stat_dim <- c(length(beam_par), length(ang_str), length(stat_par),# #
                         length(pos_par))                                    # #
      I_stat_dim <- c(length(stat_par), length(pos_par))                     # #
                                                                             # #
      beam_stat_dname <- list(beam_par, ang_str, stat_par, pos_par)          # #
      I_stat_dname <- list(stat_par, pos_par)                                # #
      rm(stat_par, pos_par)                                                  # #
                                                                             # #
      beam_stats <- array(NA, dim = beam_stat_dim,                           # #
                          dimnames = beam_stat_dname)                        # #
      beam_sky_stats<- array(NA, dim = beam_stat_dim,                        # #
                             dimnames = beam_stat_dname)                     # #
      beam_noSky_stats<- array(NA, dim = beam_stat_dim,                      # #
                               dimnames = beam_stat_dname)                   # #
                                                                             # #
      I_sky_stats<- array(NA, dim = I_stat_dim, dimnames = I_stat_dname)     # #
      I_gal_stats<- array(NA, dim = I_stat_dim, dimnames = I_stat_dname)     # #
      rm(beam_stat_dim, beam_stat_dname, I_stat_dim, I_stat_dname)           # #
                                                                             # #
      for(ang in angs){                                                      # #
        beam_stats["E", ang,, "val"] <- c(min(Es[,,"Data", ang], na.rm = T), # #
                                          max(Es[,,"Data", ang], na.rm = T), # #
                                          median(Es[,,"Data", ang], na.rm=T),# #
                                          mad(Es[,,"Data", ang], na.rm = T)) # #
        beam_stats["E", ang, "min", c("x","y")] <-                           # #
          which(Es[,,"Data", ang] ==                                         # #
                  beam_stats["E", ang, "min", "val"], arr.ind = T)[1,]       # #
        beam_stats["E", ang, "max", c("x","y")] <-                           # #
          which(Es[,,"Data", ang] ==                                         # #
                  beam_stats["E", ang, "max", "val"], arr.ind = T)[1,]       # #
                                                                             # #
        beam_stats["O", ang,, "val"] <- c(min(Os[,,"Data", ang], na.rm = T), # #
                                          max(Os[,,"Data", ang], na.rm = T), # #
                                          median(Os[,,"Data", ang], na.rm=T),# #
                                          mad(Os[,,"Data", ang], na.rm = T)) # #
        beam_stats["O", ang, "min", c("x","y")] <-                           # #
          which(Os[,,"Data", ang] ==                                         # #
                  beam_stats["O", ang, "min", "val"], arr.ind = T)[1,]       # #
        beam_stats["O", ang, "max", c("x","y")] <-                           # #
          which(Os[,,"Data", ang] ==                                         # #
                  beam_stats["O", ang, "max", "val"], arr.ind = T)[1,]       # #
                                                                             # #
        beam_sky_stats["E", ang,, "val"] <-                                  # #
          c(min(Es_sky[,,"Data", ang], na.rm = T),                           # #
            max(Es_sky[,,"Data", ang], na.rm = T),                           # #
            median(Es_sky[,,"Data", ang], na.rm = T),                        # #
            mad(Es_sky[,,"Data", ang], na.rm = T))                           # #
        beam_sky_stats["E", ang, "min", c("x","y")] <-                       # #
          which(Es_sky[,,"Data", ang] ==                                     # #
                  beam_sky_stats["E", ang, "min", "val"], arr.ind = T)[1,]   # #
        beam_sky_stats["E", ang, "max", c("x","y")] <-                       # #
          which(Es_sky[,,"Data", ang] ==                                     # #
                  beam_sky_stats["E", ang, "max", "val"], arr.ind = T)[1,]   # #
                                                                             # #
        beam_sky_stats["O", ang,, "val"] <-                                  # #
          c(min(Os_sky[,,"Data", ang], na.rm = T),                           # #
            max(Os_sky[,,"Data", ang], na.rm = T),                           # #
            median(Os_sky[,,"Data", ang], na.rm = T),                        # #
            mad(Os_sky[,,"Data",ang], na.rm = T))                            # #
        beam_sky_stats["O", ang, "min", c("x","y")] <-                       # #
          which(Os_sky[,,"Data", ang] ==                                     # #
                  beam_sky_stats["O", ang, "min", "val"], arr.ind = T)[1,]   # #
        beam_sky_stats["O", ang, "max", c("x","y")] <-                       # #
          which(Os_sky[,,"Data", ang] ==                                     # #
                  beam_sky_stats["O", ang, "max", "val"], arr.ind = T)[1,]   # #
                                                                             # #
        beam_noSky_stats["E", ang,, "val"] <-                                # #
          c(min(Es_noSky[,,"Data", ang], na.rm = T),                         # #
            max(Es_noSky[,,"Data", ang], na.rm = T),                         # #
            median(Es_noSky[,,"Data", ang], na.rm = T),                      # #
            mad(Es_noSky[,,"Data", ang], na.rm = T))                         # #
        beam_noSky_stats["E", ang, "min", c("x","y")] <-                     # #
          which(Es_noSky[,,"Data", ang] ==                                   # #
                  beam_noSky_stats["E", ang, "min", "val"], arr.ind = T)[1,] # #
        beam_noSky_stats["E", ang, "max", c("x","y")] <-                     # #
          which(Es_noSky[,,"Data", ang] ==                                   # #
                  beam_noSky_stats["E", ang, "max", "val"], arr.ind = T)[1,] # #
                                                                             # #
        beam_noSky_stats["O", ang,, "val"] <-                                # #
          c(min(Os_noSky[,,"Data", ang], na.rm = T),                         # #
            max(Os_noSky[,,"Data", ang], na.rm = T),                         # #
            median(Os_noSky[,,"Data", ang], na.rm = T),                      # #
            mad(Os_noSky[,,"Data", ang], na.rm = T))                         # #
        beam_noSky_stats["O", ang, "min", c("x","y")] <-                     # #
          which(Os_noSky[,,"Data", ang] ==                                   # #
                  beam_noSky_stats["O", ang, "min", "val"], arr.ind = T)[1,] # #
        beam_noSky_stats["O", ang, "max", c("x","y")] <-                     # #
          which(Os_noSky[,,"Data", ang] ==                                   # #
                  beam_noSky_stats["O", ang, "max", "val"], arr.ind = T)[1,] # #
      }                                                                      # #
                                                                             # #
      write.table(beam_stats, Bbobs_stats_out_csv, sep = ";", dec = ".")     # #
      write.table(beam_sky_stats, Bbsky_stats_out_csv, sep = ";", dec = ".") # #
      write.table(beam_noSky_stats, Bbhc_stats_out_csv, sep = ";", dec = ".")# #
      rm(beam_noSky_stats, beam_sky_stats, beam_stats)                       # #
                                                                             # #
      I_sky_stats[,"val"] <- c(min(bIsky[,,"Data"], na.rm=T),                # #
                               max(bIsky[,,"Data"], na.rm=T),                # #
                               median(bIsky[,,"Data"], na.rm = T),           # #
                               mad(bIsky[,,"Data"], na.rm=T))                # #
      I_sky_stats["min", c("x","y")] <-                                      # #
        which(bIsky[,,"Data"] == I_sky_stats["min", "val"], arr.ind = T)[1,] # #
      I_sky_stats["max", c("x","y")] <-                                      # #
        which(bIsky[,,"Data"] == I_sky_stats["max", "val"], arr.ind = T)[1,] # #
                                                                             # #
      I_gal_stats[,"val"] <- c(min(bIgal[,,"Data"], na.rm=T),                # #
                               max(bIgal[,,"Data"], na.rm=T),                # #
                               median(bIgal[,,"Data"], na.rm = T),           # #
                               mad(bIgal[,,"Data"], na.rm=T))                # #
      I_gal_stats["min", c("x","y")] <-                                      # #
        which(bIgal[,,"Data"] == I_gal_stats["min", "val"], arr.ind = T)[1,] # #
      I_gal_stats["max", c("x","y")] <-                                      # #
        which(bIgal[,,"Data"] == I_gal_stats["max", "val"], arr.ind = T)[1,] # #
                                                                             # #
      write.table(I_gal_stats, Ibhc_stats_out_csv, sep = ";", dec = ".")     # #
      write.table(I_sky_stats, Isky_stats_out_csv, sep = ";", dec = ".")     # #
      rm(I_gal_stats, I_sky_stats)                                           # #
      ##########################################################################
    }                                                                          #
    rm(Ibhc_stats_out_csv, exist_flag, Isky_stats_out_csv, Bbobs_stats_out_csv,#
       Bbsky_stats_out_csv, Bbhc_stats_out_csv)                                #
                                                                               #
    ############# Estimating Polarization of Observational Target ############ #
    ############################################################################
    print(paste0("* Estimating polarimetric parameters Q, U, P and X maps, ",# #
                 "using the ratio method [Bagnulo et al. 2009], for set #",  # #
                 n, " of band ", bs, " *"))                                  # #
                                                                             # #
    #Raw polarization                                                        # #
    Q_obs <- stokesPar_from_dualBeams_Bagnulo(Os[,,"Data", c(1, 3)],         # #
                                              Es[,,"Data", c(1, 3)],         # #
                                              Os[,,"Unc", c(1, 3)],          # #
                                              Es[,,"Unc", c(1, 3)])          # #
    U_obs <- stokesPar_from_dualBeams_Bagnulo(Os[,,"Data", c(2, 4)],         # #
                                              Es[,,"Data", c(2, 4)],         # #
                                              Os[,,"Unc", c(2, 4)],          # #
                                              Es[,,"Unc", c(2, 4)])          # #
    rm(Os, Es)                                                               # #
                                                                             # #
    #Raw Sky polarization                                                    # #
    Q_sky_raw <- stokesPar_from_dualBeams_Bagnulo(Os_sky[,,"Data", c(1, 3)], # #
                                                  Es_sky[,,"Data", c(1, 3)], # #
                                                  Os_sky[,,"Unc", c(1, 3)],  # #
                                                  Es_sky[,,"Unc", c(1, 3)])  # #
    U_sky_raw <- stokesPar_from_dualBeams_Bagnulo(Os_sky[,,"Data", c(2, 4)], # #
                                                  Es_sky[,,"Data", c(2, 4)], # #
                                                  Os_sky[,,"Unc", c(2, 4)],  # #
                                                  Es_sky[,,"Unc", c(2, 4)])  # #
    rm(Os_sky, Es_sky)                                                       # #
                                                                             # #
    #Sky Statistics                                                          # #
    sky_stats <- array(NA, dim = 6,                                          # #
                       dimnames = list(c("Obs_date", "Median(Isky/Igal)",    # #
                                       "Median(Psky)", "Mad(Psky)",          # #
                                       "Median(Xsky)", "Mad(Xsky)")))        # #
    Rsg <- array(NA, dim = binDim, dimnames = singNames)                     # #
                                                                             # #
    Rsg[,,"Data"] <- bIsky[,,"Data"] / bIgal[,,"Data"]                       # #
    Rsg[,,"Unc"] <- unc_div(bIsky[,,"Data"], bIgal[,,"Data"], bIsky[,,"Unc"],# #
                            bIgal[,,"Unc"], Rsg[,,"Data"])                   # #
    rm(bIsky, bIgal)                                                         # #
                                                                             # #
    sky_stats["Obs_date"] <- date_name                                       # #
    sky_stats["Median(Isky/Igal)"] <- median(Rsg[,,"Data"], na.rm = T)       # #
    rm(Rsg, date_name)                                                       # #
                                                                             # #
    #Sky (cut) and Instrumentation subtracted polarization                   # #
    Q_skyC <- stokesPar_from_dualBeams_Bagnulo(Os_noSky[,,"Data", c(1, 3)],  # #
                                               Es_noSky[,,"Data", c(1, 3)],  # #
                                               Os_noSky[,,"Unc", c(1, 3)],   # #
                                               Es_noSky[,,"Unc", c(1, 3)])   # #
    U_skyC <- stokesPar_from_dualBeams_Bagnulo(Os_noSky[,,"Data", c(2, 4)],  # #
                                               Es_noSky[,,"Data", c(2, 4)],  # #
                                               Os_noSky[,,"Unc", c(2, 4)],   # #
                                               Es_noSky[,,"Unc", c(2, 4)])   # #
    rm(Os_noSky, Es_noSky)                                                   # #
                                                                             # #
    Q_sky_instC <- array(NA, dim = dim(Q_obs), dimnames = dimnames(Q_obs))   # #
    U_sky_instC <- array(NA, dim = dim(U_obs), dimnames = dimnames(U_obs))   # #
                                                                             # #
    Q_sky_instC[,,"Data"] <- Q_skyC[,,"Data"] - bQinst[,,"Data"]             # #
    U_sky_instC[,,"Data"] <- U_skyC[,,"Data"] - bUinst[,,"Data"]             # #
    Q_sky_instC[,,"Unc"] <- unc_add(Q_skyC[,,"Unc"], bQinst[,,"Unc"])        # #
    U_sky_instC[,,"Unc"] <- unc_add(U_skyC[,,"Unc"], bUinst[,,"Unc"])        # #
    rm(Q_skyC, U_skyC)                                                       # #
                                                                             # #
    #Instrumentation (inst), MW and sky subtracted polarization              # #
    #and chromatic correction                                                # #
    Q_gal <- array(NA, dim = dim(Q_obs), dimnames = dimnames(Q_obs))         # #
    U_gal <- array(NA, dim = dim(U_obs), dimnames = dimnames(U_obs))         # #
    Q_sky_allC <- array(NA, dim = dim(Q_obs), dimnames = dimnames(Q_obs))    # #
    U_sky_allC <- array(NA, dim = dim(Q_obs), dimnames = dimnames(Q_obs))    # #
                                                                             # #
    Q_gal[,,"Data"] <- Q_sky_instC[,,"Data"] - mw_QUi['Q', "Data"]           # #
    U_gal[,,"Data"] <- U_sky_instC[,,"Data"] - mw_QUi['U', "Data"]           # #
    Q_gal[,,"Unc"] <- unc_add(Q_sky_instC[,,"Unc"], mw_QUi['Q',"Unc"]^2)     # #
    U_gal[,,"Unc"] <- unc_add(U_sky_instC[,,"Unc"], mw_QUi['U',"Unc"]^2)     # #
                                                                             # #
    Q_sky_allC[,,"Data"] <- Q_sky_raw[,,"Data"] - bQinst[,,"Data"] -         # #
      mw_QUi['Q', "Data"]                                                    # #
    U_sky_allC[,,"Data"] <- U_sky_raw[,,"Data"] - bUinst[,,"Data"] -         # #
      mw_QUi['U', "Data"]                                                    # #
    Q_sky_allC[,,"Unc"] <- unc_add(unc_add(Q_sky_raw[,,"Unc"],               # #
                                           bQinst[,,"Unc"]),                 # #
                                   mw_QUi['Q', "Unc"])                       # #
    U_sky_allC[,,"Unc"]<- unc_add(unc_add(U_sky_raw[,,"Unc"],                # #
                                          bUinst[,,"Unc"]),                  # #
                                  mw_QUi['U', "Unc"])                        # #
    rm(mw_QUi, Q_sky_raw, U_sky_raw)                                         # #
                                                                             # #
    #Target chromatic correction terms                                       # #
    Qc2X <- Q_gal[,,"Data"] * c2X                                            # #
    Qs2X <- Q_gal[,,"Data"] * s2X                                            # #
    Uc2X <- U_gal[,,"Data"] * c2X                                            # #
    Us2X <- U_gal[,,"Data"] * s2X                                            # #
    u_Qc2X <- unc_mult(Q_gal[,,"Data"], c2X, Q_gal[,,"Unc"], u_c2X)          # #
    u_Qs2X <- unc_mult(Q_gal[,,"Data"], s2X, Q_gal[,,"Unc"], u_s2X)          # #
    u_Uc2X <- unc_mult(U_gal[,,"Data"], c2X, U_gal[,,"Unc"], u_c2X)          # #
    u_Us2X <- unc_mult(U_gal[,,"Data"], s2X, U_gal[,,"Unc"], u_s2X)          # #
    rm(Q_gal, U_gal)                                                         # #
                                                                             # #
    Q_allC <- array(NA, dim = dim(Q_obs), dimnames = dimnames(Q_obs))        # #
    U_allC <- array(NA, dim = dim(U_obs), dimnames = dimnames(U_obs))        # #
                                                                             # #
    # Determining target Q and U with chromatic correction                   # #
    Q_allC[,,"Data"] <- Qc2X + Us2X                                          # #
    U_allC[,,"Data"] <- Uc2X + Qs2X                                          # #
    Q_allC[,,"Unc"] <- unc_add(u_Qc2X, u_Us2X)                               # #
    U_allC[,,"Unc"] <- unc_add(u_Uc2X, u_Qs2X)                               # #
                                                                             # #
    #Sky chromatic correction terms                                          # #
    Qc2X <- Q_sky_allC[,,"Data"] * c2X                                       # #
    Qs2X <- Q_sky_allC[,,"Data"] * s2X                                       # #
    Uc2X <- U_sky_allC[,,"Data"] * c2X                                       # #
    Us2X <- U_sky_allC[,,"Data"] * s2X                                       # #
    u_Qc2X <- unc_mult(Q_sky_allC[,,"Data"], c2X, Q_sky_allC[,,"Unc"], u_c2X)# #
    u_Qs2X <- unc_mult(Q_sky_allC[,,"Data"], s2X, Q_sky_allC[,,"Unc"], u_s2X)# #
    u_Uc2X <- unc_mult(U_sky_allC[,,"Data"], c2X, U_sky_allC[,,"Unc"], u_c2X)# #
    u_Us2X <- unc_mult(U_sky_allC[,,"Data"], s2X, U_sky_allC[,,"Unc"], u_s2X)# #
                                                                             # #
    # Determining sky Q and U with chromatic correction                      # #
    Q_sky_allC[,,"Data"] <- Qc2X + Us2X                                      # #
    U_sky_allC[,,"Data"] <- Uc2X + Qs2X                                      # #
    Q_sky_allC[,,"Unc"] <- unc_add(u_Qc2X, u_Us2X)                           # #
    U_sky_allC[,,"Unc"] <- unc_add(u_Uc2X, u_Qs2X)                           # #
    rm(Qc2X, Qs2X, Uc2X, Us2X, u_Qc2X, u_Qs2X, u_Uc2X, u_Us2X)               # #
    ############################################################################
                                                                               #
    #------------------------------------------------------------------------# #
                                                                               #
    ######################### Generating Output Files ######################## #
    ############################################################################
    print(paste0("* Creating FITS files for the different estimations of th",# #
                 "e polarimetric parameters maps of Q, U, P and X, and PDF ",# #
                 "files with arrow maps of the different estimations of the",# #
                 " polarization field (P and X); for set #", n, " of band ", # #
                 bs, " *"))                                                  # #
                                                                             # #
    pdf_list <- paste0("-Pol_", file_label)                                  # #
    pdf_list_X <- paste0("-X_", file_label)                                  # #
    pdf_list_P <- paste0("-P_", file_label)                                  # #
                                                                             # #
    Pt <- array(NA, dim = binDim, dimnames = singNames)                      # #
    Xt <- array(NA, dim = binDim, dimnames = singNames)                      # #
                                                                             # #
    rowi <- round(mean(c(1, binS)))                                          # #
    coli <- round(mean(c(1, binS)))                                          # #
    binrows <- seq(rowi, dimsObj[1], binS)                                   # #
    bincols <- seq(coli, dimsObj[2], binS)                                   # #
    Nrows <- length(binrows)                                                 # #
    Ncols <- length(bincols)                                                 # #
    binrows_r <- rep(binrows, Ncols)                                         # #
    bincols_r <- rep(bincols, each = Nrows)                                  # #
    xs <- (binrows_r - 1) / dimsObj[1]                                       # #
    ys <- (bincols_r - 1) / dimsObj[2]                                       # #
                                                                             # #
    for(out_ind in 1:length(file_label)){                                    # #
      out_pdf <- paste0(folder_t[out_ind], "PDFs/")                          # #
      if(!dir.exists(out_pdf)){                                              # #
        mkdir(out_pdf)                                                       # #
      }                                                                      # #
      out_stat <- paste0(folder_t[out_ind], "Stats/")                        # #
      if(!dir.exists(out_stat)){                                             # #
        mkdir(out_stat)                                                      # #
      }                                                                      # #
      out_img <- paste0(folder_t[out_ind], "FITS/")                          # #
      if(!dir.exists(out_img)){                                              # #
        mkdir(out_img)                                                       # #
      }                                                                      # #
                                                                             # #
      img_name <- paste0(out_img, fits_fname)                                # #
      pdf_name <- paste0(out_pdf, fits_fname, pdf_list[out_ind], ".pdf")     # #
      pdf_name_X <- paste0(out_pdf, fits_fname, pdf_list_X[out_ind], ".pdf") # #
      pdf_name_P <- paste0(out_pdf, fits_fname, pdf_list_P[out_ind], ".pdf") # #
      stat_name <- paste0(out_stat, fits_fname, file_label[out_ind],         # #
                          "-Pol_Stats.csv")                                  # #
      file_name <- paste0(file_label[out_ind],".fits")                       # #
      file_type <- file_type_1[out_ind]                                      # #
      out_label <- pdf_labels[out_ind]                                       # #
                                                                             # #
      deb_flag <- FALSE                                                      # #
      sky_flag <- FALSE                                                      # #
                                                                             # #
      switch(out_ind,                                                        # #
             c(Qt <- Q_obs, Ut <- U_obs, It <- Iobs[,,"Data"], I <- Iobs,    # #
               rm(Q_obs, U_obs), out_head <- temp_head),                     # #
             c(Qt <- Q_sky_instC, Ut <- U_sky_instC, It <- Igal[,,"Data"],   # #
               I <- Igal, rm(Q_sky_instC, U_sky_instC),                      # #
               temp_head <- modVal('CORR-SKY', 'YES',                        # #
                                   "Was sky pol. subtracted?", temp_head),   # #
               temp_head <- modVal('CORR-INS', 'YES',                        # #
                                   "Was instr. pol. subtracted?", temp_head),# #
               out_head <- temp_head),                                       # #
             c(Qt <- Q_allC, Ut <- U_allC, It <- Igal[,,"Data"],             # #
               rm(Q_allC, U_allC, Igal),                                     # #
               temp_head <- modVal('CORR-MW', 'YES',                         # #
                                   "Was MW pol. subtracted?", temp_head),    # #
               temp_head <- modVal('CORR-CHR', 'YES',                        # #
                                   "Was pol. angle corrected?", temp_head),  # #
               out_head <- temp_head, deb_flag <- TRUE, rm(temp_head)),      # #
             c(Qt <- Q_sky_allC, Ut <- U_sky_allC, It <- Isky[,,"Data"],     # #
               I <- Isky, rm(Q_sky_allC, U_sky_allC), out_head <- sky_head,  # #
               sky_flag <- TRUE, deb_flag <- TRUE, rm(sky_head))             # #
      )                                                                      # #
                                                                             # #
      print(paste0("** Calculating P and X for instance ",                   # #
                   pdf_labels[out_ind], "... **"))                           # #
                                                                             # #
      P <- P_from_QU(Qt[,,"Data"], Ut[,,"Data"], Qt[,,"Unc"], Ut[,,"Unc"])   # #
      X <- X_from_QU(Qt[,,"Data"], Ut[,,"Data"], Qt[,,"Unc"], Ut[,,"Unc"])   # #
                                                                             # #
      pol_stats <- array(NA, dim = c(4, 4, 3),                               # #
                         dimnames = list(c("Q", "U", "P", "X"),              # #
                                         c("min", "max", "median", "mad"),   # #
                                         c("val", "x", "y")))                # #
                                                                             # #
      P_NA <- which(is.na(P$Data), arr.ind = T)                              # #
      P$Unc[P_NA] <- NA                                                      # #
      X$Data[P_NA] <- NA                                                     # #
      X$Unc[P_NA] <- NA                                                      # #
      P_flg <- length(P_NA) / 2 != length(P$Data)                            # #
      rm(P_NA)                                                               # #
                                                                             # #
      if(P_flg){                                                             # #
        Pf_med <- median(P$Data, na.rm = T)                                  # #
        pol_stats["P","mad","val"] <- mad(P$Data, na.rm = T)                 # #
        Pf_unc <- unc_median(P$Data, P$Unc, 0)                               # #
        Pf_max <- max(P$Data, na.rm = T)                                     # #
        Pf_min <- min(P$Data, na.rm = T)                                     # #
        pol_stats["P",c("min","max","median"),"val"] <-                      # #
          c(Pf_min, Pf_max, Pf_med)                                          # #
                                                                             # #
        Xf_max <- max(X$Data, na.rm = T)                                     # #
        Xf_min <- min(X$Data, na.rm = T)                                     # #
        Xf_med <- median(X$Data, na.rm = T)                                  # #
        pol_stats["X","mad","val"] <- mad(X$Data, na.rm = T)                 # #
        Xf_unc <- unc_median(X$Data, X$Unc, 0)                               # #
                                                                             # #
        pol_stats["X",c("min","max","median"),"val"] <-                      # #
          c(Xf_min,Xf_max,Xf_med)                                            # #
        pol_stats["Q",c("min","max","median","mad"),"val"] <-                # #
          c(min(Qt[,,"Data"], na.rm = T), max(Qt[,,"Data"], na.rm = T),      # #
            median(Qt[,,"Data"], na.rm = T), mad(Qt[,,"Data"], na.rm = T))   # #
        pol_stats["U",c("min","max","median","mad"),"val"] <-                # #
          c(min(Ut[,,"Data"], na.rm = T), max(Ut[,,"Data"], na.rm = T),      # #
            median(Ut[,,"Data"], na.rm = T), mad(Ut[,,"Data"], na.rm = T))   # #
      }                                                                      # #
                                                                             # #
      Pt[,,"Data"] <- P$Data                                                 # #
      Pt[,,"Unc"] <- P$Unc                                                   # #
      rm(P)                                                                  # #
                                                                             # #
      Xt[,,"Data"] <- X$Data                                                 # #
      Xt[,,"Unc"] <- X$Unc                                                   # #
      rm(X)                                                                  # #
                                                                             # #
      Q <- bin_expand_map(Qt, binS, c(dimsObj, 2), ref_crpix)                # #
                                                                             # #
      Qtemp <- Q[,,"Data"]                                                   # #
      Qtemp[NAs_ind] <- NA                                                   # #
      Q[,,"Data"] <- Qtemp                                                   # #
                                                                             # #
      if(P_flg){                                                             # #
        pol_stats["Q","min",c("x","y")] <-                                   # #
          which(Q[,,"Data"] == pol_stats["Q","min","val"], arr.ind = T)[1,]  # #
        pol_stats["Q","max",c("x","y")] <-                                   # #
          which(Q[,,"Data"] == pol_stats["Q","max","val"], arr.ind = T)[1,]  # #
      }                                                                      # #
                                                                             # #
      Qtemp <- Q[,,"Unc"]                                                    # #
      Qtemp[NAs_ind] <- NA                                                   # #
      Q[,,"Unc"] <- Qtemp                                                    # #
      rm(Qtemp)                                                              # #
                                                                             # #
      out_head <- modVal('FILETYP1', file_type, "File type", out_head)       # #
      out_head <- modVal('FILETYP2', 'STOKES Q', "File sub-type", out_head)  # #
      writeFITSim(Q, file = paste0(img_name, "-Q_", file_name),              # #
                  crvaln = crval, crpixn = crpix, ctypen = ctype,            # #
                  header = out_head)                                         # #
      rm(Q)                                                                  # #
                                                                             # #
      U <- bin_expand_map(Ut, binS, c(dimsObj, 2), ref_crpix)                # #
                                                                             # #
      Utemp <- U[,,"Data"]                                                   # #
      Utemp[NAs_ind] <- NA                                                   # #
      U[,,"Data"] <- Utemp                                                   # #
                                                                             # #
      if(P_flg){                                                             # #
        pol_stats["U","min",c("x","y")] <-                                   # #
          which(U[,,"Data"] == pol_stats["U","min","val"], arr.ind = T)[1,]  # #
        pol_stats["U","max",c("x","y")] <-                                   # #
          which(U[,,"Data"] == pol_stats["U","max","val"], arr.ind = T)[1,]  # #
      }                                                                      # #
                                                                             # #
      Utemp <- U[,,"Unc"]                                                    # #
      Utemp[NAs_ind] <- NA                                                   # #
      U[,,"Unc"] <- Utemp                                                    # #
      rm(Utemp)                                                              # #
                                                                             # #
      out_head <- modVal('FILETYP2', 'STOKES U', "File sub-type", out_head)  # #
      writeFITSim(U, file = paste0(img_name, "-U_", file_name),              # #
                  crvaln = crval, crpixn = crpix, ctypen = ctype,            # #
                  header = out_head)                                         # #
      rm(U)                                                                  # #
                                                                             # #
      P <- array(NA, dim = pxDim, dimnames = pxNames)                        # #
      X <- array(NA, dim = pxDim, dimnames = pxNames)                        # #
                                                                             # #
      tP <- bin_expand_map(Pt, binS, c(dimsObj, 2), ref_crpix)               # #
                                                                             # #
      Ptemp <- tP[,,"Data"]                                                  # #
      Ptemp[NAs_ind] <- NA                                                   # #
      P[,,"Data"] <- Ptemp                                                   # #
                                                                             # #
      if(P_flg){                                                             # #
        pol_stats["P","min",c("x","y")] <- which(P[,,"Data"] == Pf_min,      # #
                                                 arr.ind = T)[1,]            # #
        pol_stats["P","max",c("x","y")] <- which(P[,,"Data"] == Pf_max,      # #
                                                 arr.ind = T)[1,]            # #
      }                                                                      # #
                                                                             # #
      Ptemp <- tP[,,"Unc"]                                                   # #
      Ptemp[NAs_ind] <- NA                                                   # #
      P[,,"Unc"] <- Ptemp                                                    # #
      P[,,"Unc(%)"] <- P[,,"Unc"] / P[,,"Data"] * 100                        # #
      rm(Ptemp, tP)                                                          # #
                                                                             # #
      tI <- bin_code_map(I, binS, ref_crpix)                                 # #
      tIP <- array(NA, dim = bin_dims[1:3], dimnames = singNames)            # #
      tIP[,,"Data"] <- Pt[,,"Data"] * tI[,,"Data"]                           # #
      tIP[,,"Unc"] <- unc_mult(Pt[,,"Data"], tI[,,"Data"], Pt[,,"Unc"],      # #
                               tI[,,"Unc"])                                  # #
      IP <-  bin_expand_map(tIP, binS, c(dimsObj, 2), ref_crpix)             # #
      rm(tIP)                                                                # #
                                                                             # #
      temp <- IP[,,"Data"]                                                   # #
      temp[NAs_ind] <- NA                                                    # #
      IP[,,"Data"] <- temp                                                   # #
      temp <- IP[,,"Unc"]                                                    # #
      temp[NAs_ind] <- NA                                                    # #
      IP[,,"Unc"] <- temp                                                    # #
      IPg_med <- median((IP[,,"Data"] * maskGal), na.rm = T)                 # #
      IPg_max <- max((IP[,,"Data"] * maskGal), na.rm = T)                    # #
      IPg_min <- min((IP[,,"Data"] * maskGal), na.rm = T)                    # #
      IPf_med <- median(IP[,,"Data"], na.rm = T)                             # #
      IPf_max <- max(IP[,,"Data"], na.rm = T)                                # #
      IPf_min <- min(IP[,,"Data"], na.rm = T)                                # #
      rm(temp)                                                               # #
                                                                             # #
      out_head <- modVal('FILETYP2', 'POL DEG P', "File sub-type", out_head) # #
      writeFITSim(P, file = paste0(img_name, "-P_", file_name),              # #
                  crvaln = crval, crpixn = crpix, ctypen = ctype,            # #
                  header = out_head)                                         # #
                                                                             # #
      out_head <- modVal('FILETYP2', 'POL IxP', "File sub-type", out_head)   # #
      writeFITSim(IP, file = paste0(img_name, "-P_Flux-", file_name),        # #
                  crvaln = crval, crpixn = crpix, ctypen = ctype,            # #
                  header = out_head)                                         # #
                                                                             # #
      tX <- bin_expand_map(Xt, binS, c(dimsObj, 2), ref_crpix)               # #
                                                                             # #
      Xtemp <- tX[,,"Data"]                                                  # #
      Xtemp[NAs_ind] <- NA                                                   # #
      X[,,"Data"] <- Xtemp                                                   # #
                                                                             # #
      if(P_flg){                                                             # #
        pol_stats["X","min",c("x","y")] <- which(X[,,"Data"] == Xf_min,      # #
                                                 arr.ind = T)[1,]            # #
        pol_stats["X","max",c("x","y")] <- which(X[,,"Data"] == Xf_max,      # #
                                                 arr.ind = T)[1,]            # #
      }                                                                      # #
                                                                             # #
      write.table(pol_stats, stat_name, sep = ";", dec = ".")                # #
                                                                             # #
      Xtemp <- tX[,,"Unc"]                                                   # #
      Xtemp[NAs_ind] <- NA                                                   # #
      X[,,"Unc"] <- Xtemp                                                    # #
      rm(Xtemp, tX)                                                          # #
                                                                             # #
      X[,,"Unc(%)"] <- abs(X[,,"Unc"] / X[,,"Data"]) * 100                   # #
                                                                             # #
      out_head <- modVal('FILETYP2', 'POL ANG X', "File sub-type", out_head) # #
      writeFITSim(X, file = paste0(img_name, "-X_", file_name),              # #
                  crvaln = crval, crpixn = crpix, ctypen = ctype,            # #
                  header = out_head)                                         # #
                                                                             # #
      print(paste0("** The FITS files for the parameters Q, U, P and X for ",# #
                   "instance '", pdf_labels[out_ind], "' have been created.",# #
                   " **"))                                                   # #
                                                                             # #
      dif_x <- crpix[1] - 0:dimsObj[1]                                       # #
      dif_y <- 0:dimsObj[2] - crpix[2]                                       # #
      raseq <- ref_coords[1] + dif_x * scl_x / cosd                          # #
      decseq <- ref_coords[2] + dif_y * scl_y                                # #
                                                                             # #
      ra_bot <- round(raseq[round(0.1 * dimsObj[1])], 2)                     # #
      ra_top <- round(raseq[round(0.9 * dimsObj[1])], 2)                     # #
      dec_bot <- round(decseq[round(0.235 * dimsObj[2])], 2)                 # #
      dec_top <- round(decseq[round(0.935 * dimsObj[2])], 2)                 # #
                                                                             # #
      xticks <- seq(0.1, .9, length.out = 5)                                 # #
      yticks <- seq(0.235, .935, length.out = 5)                             # #
      leg_X_ticks <- seq(-90, 90, 30)                                        # #
      leg_P_ticks <- seq(0, max(gP_crp[,,1], na.rm = T), length.out = 7)*100 # #
      xlabs <- round(seq(ra_bot, ra_top, length.out = 5), 2)                 # #
      ylabs <- round(seq(dec_bot, dec_top, length.out = 5), 2)               # #
      leg_X_labs <- paste0(leg_X_ticks, "°")                                 # #
      leg_P_labs <- paste0(leg_P_ticks, "%")                                 # #
                                                                             # #
      # POLARIZATION ANGLE X, PDF                                            # #
      pdf(pdf_name_X, width = 30, height = 30)                               # #
      par(cex.main = 6)                                                      # #
      par(cex.axis = 4)                                                      # #
      par(cex.lab = 5)                                                       # #
                                                                             # #
      # Generate a circular color scale with 10 colors                       # #
      num_colors <- 181                                                      # #
      hue_range <- c(0, 360)                                                 # #
      hues <- seq(hue_range[1], hue_range[2], length.out = num_colors)       # #
      colors <- hcl(h = hues, c = 180, l = 65, alpha = 1)                    # #
      rm(num_colors, hue_range, hues)                                        # #
                                                                             # #
      # Main Plot                                                            # #
      image.plot(X[,,1], xlab = ' ', ylab = ' ', xaxt = "n", yaxt = "n",     # #
                 col = colors, bigplot = c(.09, .06 + delt, .09, .95),       # #
                 smallplot = c(.07 + delt, .08 + delt, .09, .95),            # #
                 legend.lab = "Pol. Ang. (deg)", legend.cex = 6,             # #
                 legend.line = 15, axis.args = list(at = leg_X_ticks,        # #
                                                    labels = leg_X_labs))    # #
      # Main title                                                           # #
      mtext(paste0("Polarization Angle Map, ", out_label), cex = 4, line = 2)# #
      # X-axis                                                               # #
      axis(1, at = xticks, labels = paste0(xlabs, "°"), line = 3, lwd = 0)   # #
      axis(1, at = xticks, labels = rep("", length(xticks)), lwd.ticks = 1)  # #
      mtext("Ra (deg)", side = 1, cex = 6, line = 11)                        # #
      # Y-axis                                                               # #
      axis(2, at = yticks, labels = paste0(ylabs, "°"))                      # #
      mtext("Dec (deg)", side = 2, cex = 6, line = 8)                        # #
      # N/E referencial                                                      # #
      arrows(x0 = delt+.1, y0 = 0.05, x1 = delt+.1, y1 = 0.1, length = 0.25, # #
             col = "red", lwd = 5)                                           # #
      arrows(x0 = delt+.1, y0 = 0.05, x1 = delt+.1 - 0.05, y1 = 0.05,        # #
             length = 0.25, col = "red", lwd = 5)                            # #
      text(x = delt+.1, y = 0.105, labels = "N", pos = 3, col = "red",       # #
           cex = 3)                                                          # #
      text(x = delt+.1-0.053, y = 0.05, labels = "E", pos = 2, col = "red",  # #
           cex = 3)                                                          # #
      dev.off()                                                              # #
                                                                             # #
      # POLARIZATION ANGLE P, PDF                                            # #
      pdf(pdf_name_P, width = 30, height = 30)                               # #
      par(cex.main = 6)                                                      # #
      par(cex.axis = 4)                                                      # #
      par(cex.lab = 5)                                                       # #
                                                                             # #
      # Generate a circular color scale with 10 colors                       # #
      num_colors <- 181                                                      # #
      hue_range <- c(0, 270)                                                 # #
      hues <- seq(hue_range[1], hue_range[2], length.out = num_colors)       # #
      colors <- hcl(h = hues, c = 180, l = 65, alpha = 1)                    # #
      rm(num_colors, hue_range, hues)                                        # #
                                                                             # #
      # Main Plot                                                            # #
      image.plot(P[,,1] * 100, xlab = ' ', ylab = ' ', xaxt = "n", yaxt ="n",# #
                 col = colors, bigplot = c(.09, .06 + delt, .09, .95),       # #
                 smallplot = c(.07 + delt, .08 + delt, .09, .95),            # #
                 legend.lab = "Pol. Degree (%)", legend.cex = 6,             # #
                 legend.line = 15, axis.args = list(at = leg_P_ticks,        # #
                                                    labels = leg_P_labs))    # #
      # Main title                                                           # #
      mtext(paste0("Polarization Degree Map, ", out_label), cex = 4, line =2)# #
      # X-axis                                                               # #
      axis(1, at = xticks, labels = paste0(xlabs, "°"), line = 3, lwd = 0)   # #
      axis(1, at = xticks, labels = rep("", length(xticks)), lwd.ticks = 1)  # #
      mtext("Ra (deg)", side = 1, cex = 6, line = 11)                        # #
      # Y-axis                                                               # #
      axis(2, at = yticks, labels = paste0(ylabs, "°"))                      # #
      mtext("Dec (deg)", side = 2, cex = 6, line = 8)                        # #
      # N/E referencial                                                      # #
      arrows(x0 = delt+.1, y0 = 0.05, x1 = delt+.1, y1 = 0.1, length = 0.25, # #
             col = "red", lwd = 5)                                           # #
      arrows(x0 = delt+.1, y0 = 0.05, x1 = delt+.1 - 0.05, y1 = 0.05,        # #
             length = 0.25, col = "red", lwd = 5)                            # #
      text(x = delt+.1, y = 0.105, labels = "N", pos = 3, col = "red",       # #
           cex = 3)                                                          # #
      text(x = delt+.1-0.053, y = 0.05, labels = "E", pos = 2, col = "red",  # #
           cex = 3)                                                          # #
      dev.off()                                                              # #
                                                                             # #
      legticks <- round(seq(min(It[which(!is.infinite(It))], na.rm = T),     # #
                            max(It[which(!is.infinite(It))], na.rm = T),     # #
                            length.out = 5))                                 # #
      leglabs <- paste0(legticks)                                            # #
                                                                             # #
      # POLARIZATION VECTOR FIELD, PDF                                       # #
      pdf(pdf_name, width = 30, height = 30, bg = "black")                   # #
      par(cex.main = 6)                                                      # #
      par(cex.axis = 4)                                                      # #
      par(col.axis = "white")                                                # #
      par(cex.lab = 5)                                                       # #
                                                                             # #
      temp_X <- as.vector((X[,,"Data"] * maskGal)[binrows, bincols])         # #
      temp_P <- as.vector((P[,,"Data"] * maskGal)[binrows, bincols])         # #
      temp_IP <- as.vector((IP[,,"Data"] * maskGal)[binrows, bincols])       # #
      temp_u_X <- as.vector((X[,,"Unc"] * maskGal)[binrows, bincols])        # #
      temp_u_P <- as.vector((P[,,"Unc"] * maskGal)[binrows, bincols])        # #
      temp_u_IP <- as.vector((IP[,,"Unc"] * maskGal)[binrows, bincols])      # #
      rm(P, IP)                                                              # #
                                                                             # #
      min_P <- round(min(temp_P, na.rm = T) * 100, 2)                        # #
      med_P <- round(median(temp_P, na.rm = T) * 100, 2)                     # #
      unc_med_P <- round(unc_median(temp_P, temp_u_P, 0) * 100, 2)           # #
      max_P <- round(max(temp_P, na.rm = T) * 100, 2)                        # #
      min_IP <- round(min(temp_IP, na.rm = T), 2)                            # #
      med_IP <- round(median(temp_IP, na.rm = T), 2)                         # #
      unc_med_IP <- round(unc_median(temp_IP, temp_u_IP, 0), 2)              # #
      max_IP <- round(max(temp_IP, na.rm = T), 2)                            # #
      min_X <- round(min(temp_X, na.rm = T), 2)                              # #
      med_X <- round(median(temp_X, na.rm = T), 2)                           # #
      unc_med_X <- round(unc_median(temp_X, temp_u_X, 0), 2)                 # #
      max_X <- round(max(temp_X, na.rm = T), 2)                              # #
                                                                             # #
      na_inds <- union(which(is.na(temp_P), arr.ind = TRUE),                 # #
                       which(is.na(temp_X), arr.ind = TRUE))                 # #
                                                                             # #
      if(length(na_inds) != 0){                                              # #
        temp_x <- xs[-na_inds]                                               # #
        temp_y <- ys[-na_inds]                                               # #
        temp_P <- temp_P[-na_inds]                                           # #
        temp_IP <- temp_IP[-na_inds]                                         # #
        temp_X <- temp_X[-na_inds]                                           # #
      }                                                                      # #
      rm(na_inds)                                                            # #
                                                                             # #
      bkg <- asinh(It)                                                       # #
      bkg[which(is.na(bkg), arr.ind = T)] <- 0                               # #
      legticks <- round(seq(min(bkg[which(!is.infinite(bkg))], na.rm = T),   # #
                            max(bkg[which(!is.infinite(bkg))], na.rm = T),   # #
                            length.out = 5))                                 # #
      leglabs <- paste0(legticks)                                            # #
                                                                             # #
      # Main Plot                                                            # #
      image.plot(bkg, xlab = ' ', ylab = ' ', xaxt = "n", yaxt = "n",        # #
                 col = gray.colors(64, 0, 1, 0.6),                           # #
                 bigplot = c(.1, .1 + delt, .1, .9),                         # #
                 smallplot = c(.11 + delt, .14 + delt, .1, .9),              # #
                 legend.lab = "asinh(Counts)", legend.cex = 6,               # #
                 legend.line = 15, axis.args = list(at = legticks,           # #
                                                    labels = leglabs,        # #
                                                    col = "white"))          # #
      # Arrow Field                                                          # #
      vectorField(temp_X, temp_P, temp_x, temp_y, scale = .5, headspan = 0,  # #
                  vecspec = "deg", col = "green")                            # #
      # Center Points                                                        # #
      text(tgt_xy[1] / dimsObj[1], tgt_xy[2] / dimsObj[2], "+", col = "cyan",# #
           cex = 4)                                                          # #
      text(tgt_x / dimsObj[1], tgt_y / dimsObj[2], "x", col = "yellow",      # #
           cex = 4)                                                          # #
      # Main Title                                                           # #
      mtext(paste0("Polarization Map, ", out_label), cex = 4, col = "white", # #
            line = 2)                                                        # #
      # X-axis                                                               # #
      axis(1, at = xticks, labels = paste0(xlabs, "°"), line = 3, lwd = 0)   # #
      axis(1, at = xticks, labels = rep("", length(xticks)), lwd = 2,        # #
           lwd.ticks = 1)                                                    # #
      mtext("Ra (deg)", side = 1, cex = 6, line = 11, col = "white")         # #
      # Y-axis                                                               # #
      axis(2, at = yticks, labels = paste0(ylabs, "°"), lwd = 2)             # #
      mtext("Dec (deg)", side = 2, cex = 6, line = 8, col = "white")         # #
      # N/E referencial                                                      # #
      arrows(x0 = delt+.1, y0 = 0.05, x1 = delt+.1, y1 = 0.1, length = 0.25, # #
             col = "red", lwd = 5)                                           # #
      arrows(x0 = delt+.1, y0 = 0.05, x1 = delt+.1 - 0.05, y1 = 0.05,        # #
             length = 0.25, col = "red", lwd = 5)                            # #
      text(x = delt+.1, y = 0.105, labels = "N", pos = 3, col = "red",       # #
           cex = 3)                                                          # #
      text(x = delt+.1-0.053, y = 0.05, labels = "E", pos = 2, col = "red",  # #
           cex = 3)                                                          # #
      # Statistical info                                                     # #
      mtext(paste0("min P_gal = ", min_P, "%, median P_gal = ", med_P,       # #
                   "% +- ", unc_med_P, "%, max P_gal = ", max_P, "%"),       # #
            side = 1, line = -8, cex = 2, col = "white")                     # #
      mtext(paste0("min X_gal = ", min_X, "º, median X_gal = ", med_X,       # #
                   "º +- ", unc_med_X, "º, max X_gal = ", max_X, "º"),       # #
            side = 1, line = -6, cex = 2, col = "white")                     # #
                                                                             # #
      # Main Plot                                                            # #
      image.plot(bkg, xlab = ' ', ylab = ' ', xaxt = "n", yaxt = "n",        # #
                 col = gray.colors(64, 0, 1, 0.6),                           # #
                 bigplot = c(.1, .1 + delt, .1, .9),                         # #
                 smallplot = c(.11 + delt, .14 + delt, .1, .9),              # #
                 legend.lab = "asinh(Counts)", legend.cex = 6,               # #
                 legend.line = 15, axis.args = list(at = legticks,           # #
                                                    labels = leglabs,        # #
                                                    col = "white"))          # #
      # Arrow Field                                                          # #
      vectorField(temp_X, temp_IP, temp_x, temp_y, scale = .5, headspan = 0, # #
                  vecspec = "deg", col = "green")                            # #
      # Center Points                                                        # #
      text(tgt_xy[1] / dimsObj[1], tgt_xy[2] / dimsObj[2], "+", col = "cyan",# #
           cex = 4)                                                          # #
      text(tgt_x / dimsObj[1], tgt_y / dimsObj[2], "x", col = "yellow",      # #
           cex = 4)                                                          # #
      # Main Title                                                           # #
      mtext(paste0("Polarized Flux Map, ", out_label), cex = 4, line = 2,    # #
            col = "white")                                                   # #
      # X-axis                                                               # #
      axis(1, at = xticks, labels = paste0(xlabs, "°"), line = 3, lwd = 0)   # #
      axis(1, at = xticks, labels = rep("", length(xticks)), lwd = 2,        # #
           lwd.ticks = 1)                                                    # #
      mtext("Ra (deg)", side = 1, cex = 6, line = 11, col = "white")         # #
      # Y-axis                                                               # #
      axis(2, at = yticks, labels = paste0(ylabs, "°"), lwd = 2)             # #
      mtext("Dec (deg)", side = 2, cex = 6, line = 8, col = "white")         # #
      # N/E referencial                                                      # #
      arrows(x0 = delt+.1, y0 = 0.05, x1 = delt+.1, y1 = 0.1, length = 0.25, # #
             col = "red", lwd = 5)                                           # #
      arrows(x0 = delt+.1, y0 = 0.05, x1 = delt+.1 - 0.05, y1 = 0.05,        # #
             length = 0.25, col = "red", lwd = 5)                            # #
      text(x = delt+.1, y = 0.105, labels = "N", pos = 3, col = "red",       # #
           cex = 3)                                                          # #
      text(x = delt+.1-0.053, y = 0.05, labels = "E", pos = 2, col = "red",  # #
           cex = 3)                                                          # #
      # Statistical info                                                     # #
      mtext(paste0("min IP_gal = ", min_IP, "%, median IP_gal = ", med_IP,   # #
                   "% +- ", unc_med_IP, "%, max IP_gal = ", max_IP, "%"),    # #
            side = 1, line = -8, cex = 2, col = "white")                     # #
      mtext(paste0("min X_gal = ", min_X, "º, median X_gal = ", med_X,       # #
                   "º +- ", unc_med_X, "º, max X_gal = ", max_X, "º"),       # #
            side = 1, line = -6, cex = 2, col = "white")                     # #
      dev.off()                                                              # #
                                                                             # #
      if(deb_flag){                                                          # #
        print(paste0("** Debiasing P for instance ", pdf_labels[out_ind],    # #
                     "... **"))                                              # #
        P <- debiased_P_from_QU(Qt[,,"Data"], Ut[,,"Data"], Qt[,,"Unc"],     # #
                                Ut[,,"Unc"])                                 # #
        X <- X_from_QU(Qt[,,"Data"], Ut[,,"Data"], Qt[,,"Unc"], Ut[,,"Unc"]) # #
                                                                             # #
        P_NA <- which(is.na(P$Data), arr.ind = T)                            # #
        P$Unc[P_NA] <- NA                                                    # #
        P_flg <- length(P_NA) / 2 != length(P$Data)                          # #
        rm(P_NA)                                                             # #
                                                                             # #
        if(P_flg){                                                           # #
          Pf_med <- median(P$Data, na.rm = T)                                # #
          pol_stats["P","mad","val"] <- mad(P$Data, na.rm = T)               # #
          Pf_unc <- unc_median(P$Data, P$Unc, 0)                             # #
          Pf_max <- max(P$Data, na.rm = T)                                   # #
          Pf_min <- min(P$Data, na.rm = T)                                   # #
          pol_stats["P",c("min","max","median"),"val"] <-                    # #
            c(Pf_min, Pf_max, Pf_med)                                        # #
        }                                                                    # #
                                                                             # #
        if(sky_flag){                                                        # #
          sky_stats["Median(Psky)"] <- pol_stats["P", "median", "val"]       # #
          sky_stats["Mad(Psky)"] <- pol_stats["P", "mad", "val"]             # #
          sky_stats["Median(Xsky)"] <- pol_stats["X", "median", "val"]       # #
          sky_stats["Mad(Xsky)"] <- pol_stats["X", "mad", "val"]             # #
                                                                             # #
          stat_name <- paste0(out_stat, fits_fname, "-sky_Stats.csv")        # #
                                                                             # #
          write.table(sky_stats, stat_name, sep = ";", dec = ".")            # #
          rm(sky_stats, stat_name)                                           # #
        }                                                                    # #
                                                                             # #
        Pt[,,"Data"] <- P$Data                                               # #
        Pt[,,"Unc"] <- P$Unc                                                 # #
        Xt[,,"Data"] <- X$Data                                               # #
        Xt[,,"Unc"] <- X$Unc                                                 # #
        rm(P,X)                                                              # #
                                                                             # #
        P <- array(NA, dim = pxDim, dimnames = pxNames)                      # #
        X <- array(NA, dim = pxDim, dimnames = pxNames)                      # #
                                                                             # #
        tP <- bin_expand_map(Pt, binS, c(dimsObj, 2), ref_crpix)             # #
        tX <- bin_expand_map(Xt, binS, c(dimsObj, 2), ref_crpix)             # #
                                                                             # #
        Ptemp <- tP[,,"Data"]                                                # #
        Ptemp[NAs_ind] <- NA                                                 # #
        P[,,"Data"] <- Ptemp                                                 # #
        Xtemp <- tX[,,"Data"]                                                # #
        Xtemp[NAs_ind] <- NA                                                 # #
        X[,,"Data"] <- Xtemp                                                 # #
        if(P_flg){                                                           # #
          pol_stats["P","min",c("x","y")] <- which(P[,,"Data"] == Pf_min,    # #
                                                   arr.ind = T)[1,]          # #
          pol_stats["P","max",c("x","y")] <- which(P[,,"Data"] == Pf_max,    # #
                                                   arr.ind = T)[1,]          # #
          pol_stats["X","min",c("x","y")] <- which(X[,,"Data"] == Pf_min,    # #
                                                   arr.ind = T)[1,]          # #
          pol_stats["X","max",c("x","y")] <- which(X[,,"Data"] == Pf_max,    # #
                                                   arr.ind = T)[1,]          # #
        }                                                                    # #
                                                                             # #
        Ptemp <- tP[,,"Unc"]                                                 # #
        Ptemp[NAs_ind] <- NA                                                 # #
        P[,,"Unc"] <- Ptemp                                                  # #
        P[,,"Unc(%)"] <- P[,,"Unc"] / P[,,"Data"] * 100                      # #
        Xtemp <- tX[,,"Unc"]                                                 # #
        Xtemp[NAs_ind] <- NA                                                 # #
        X[,,"Unc"] <- Xtemp                                                  # #
        X[,,"Unc(%)"] <- X[,,"Unc"] / X[,,"Data"] * 100                      # #
        rm(Ptemp, tP, Xtemp, tX)                                             # #
                                                                             # #
        tI <- bin_code_map(I, binS, ref_crpix)                               # #
        tIP <- array(NA, dim = bin_dims[1:3], dimnames = singNames)          # #
        tIP[,,"Data"] <- Pt[,,"Data"] * tI[,,"Data"]                         # #
        tIP[,,"Unc"] <- unc_mult(Pt[,,"Data"], tI[,,"Data"], Pt[,,"Unc"],    # #
                                 tI[,,"Unc"])                                # #
        IP <-  bin_expand_map(tIP, binS, c(dimsObj, 2), ref_crpix)           # #
        rm(tIP)                                                              # #
                                                                             # #
        temp <- IP[,,"Data"]                                                 # #
        temp[NAs_ind] <- NA                                                  # #
        IP[,,"Data"] <- temp                                                 # #
        temp <- IP[,,"Unc"]                                                  # #
        temp[NAs_ind] <- NA                                                  # #
        IP[,,"Unc"] <- temp                                                  # #
        rm(temp)                                                             # #
                                                                             # #
        out_head <- modVal('FILETYP2', 'POL DEG P', "File sub-type",         # #
                           out_head)                                         # #
        out_head <- modVal('CORR-PBIAS', 'YES', "Was pol. bias corrected?",  # #
                           out_head)                                         # #
        writeFITSim(P, file = paste0(img_name, "-debP_", file_name),         # #
                    crvaln = crval, crpixn = crpix, ctypen = ctype,          # #
                    header = out_head)                                       # #
                                                                             # #
        out_head <- modVal('FILETYP2', 'POL IxP', "File sub-type", out_head) # #
        writeFITSim(IP, file = paste0(img_name, "-debP_Flux-", file_name),   # #
                    crvaln = crval, crpixn = crpix, ctypen = ctype,          # #
                    header = out_head)                                       # #
                                                                             # #
        out_head <- modVal('CORR-PBIAS', 'NO', "Was pol. bias corrected?",   # #
                           out_head)                                         # #
                                                                             # #
        stat_name <- paste0(out_stat, fits_fname, file_label[out_ind],       # #
                            "-debPol_Stats.csv")                             # #
                                                                             # #
        write.table(pol_stats, stat_name, sep = ";", dec = ".")              # #
                                                                             # #
        print(paste0("** The FITS file of debiased P for instance '",        # #
                     pdf_labels[out_ind], "' has been created... **"))       # #
                                                                             # #
        pdf_name <- paste0(out_pdf, fits_fname, pdf_list[out_ind],"-deb.pdf")# #
                                                                             # #
        # POLARIZATION VECTOR FIELD, PDF                                     # #
        pdf(pdf_name, width = 30, height = 30, bg = "black")                 # #
        par(cex.main = 6)                                                    # #
        par(cex.axis = 4)                                                    # #
        par(col.axis = "white")                                              # #
        par(cex.lab = 5)                                                     # #
                                                                             # #
        temp_X <- as.vector((X[,,"Data"] * maskGal)[binrows, bincols])       # #
        temp_P <- as.vector((P[,,"Data"] * maskGal)[binrows, bincols])       # #
        temp_IP <- as.vector((IP[,,"Data"] * maskGal)[binrows, bincols])     # #
        temp_u_X <- as.vector((X[,,"Unc"] * maskGal)[binrows, bincols])      # #
        temp_u_P <- as.vector((P[,,"Unc"] * maskGal)[binrows, bincols])      # #
        temp_u_IP <- as.vector((IP[,,"Unc"] * maskGal)[binrows, bincols])    # #
        rm(P, IP, X)                                                         # #
                                                                             # #
        min_P <- round(min(temp_P, na.rm = T) * 100, 2)                      # #
        med_P <- round(median(temp_P, na.rm = T) * 100, 2)                   # #
        unc_med_P <- round(unc_median(temp_P, temp_u_P, 0) * 100, 2)         # #
        max_P <- round(max(temp_P, na.rm = T) * 100, 2)                      # #
        min_IP <- round(min(temp_IP, na.rm = T), 2)                          # #
        med_IP <- round(median(temp_IP, na.rm = T), 2)                       # #
        unc_med_IP <- round(unc_median(temp_IP, temp_u_IP, 0), 2)            # #
        max_IP <- round(max(temp_IP, na.rm = T), 2)                          # #
        min_X <- round(min(temp_X, na.rm = T), 2)                            # #
        med_X <- round(median(temp_X, na.rm = T), 2)                         # #
        unc_med_X <- round(unc_median(temp_X, temp_u_X, 0), 2)               # #
        max_X <- round(max(temp_X, na.rm = T), 2)                            # #
                                                                             # #
        na_inds <- union(which(is.na(temp_P), arr.ind = TRUE),               # #
                         which(is.na(temp_X), arr.ind = TRUE))               # #
                                                                             # #
        if(length(na_inds) != 0){                                            # #
          temp_x <- xs[-na_inds]                                             # #
          temp_y <- ys[-na_inds]                                             # #
          temp_P <- temp_P[-na_inds]                                         # #
          temp_IP <- temp_IP[-na_inds]                                       # #
          temp_X <- temp_X[-na_inds]                                         # #
        }                                                                    # #
        rm(na_inds)                                                          # #
                                                                             # #
        bkg <- asinh(It)                                                     # #
        bkg[which(is.na(bkg), arr.ind = T)] <- 0                             # #
        legticks <- round(seq(min(bkg[which(!is.infinite(bkg))], na.rm = T), # #
                              max(bkg[which(!is.infinite(bkg))], na.rm = T), # #
                              length.out = 5))                               # #
        leglabs <- paste0(legticks)                                          # #
                                                                             # #
                                                                             # #
        # Main Plot                                                          # #
        image.plot(bkg, xlab = ' ', ylab = ' ', xaxt = "n", yaxt = "n",      # #
                   col = gray.colors(64, 0, 1, 0.6),                         # #
                   bigplot = c(.1, .1 + delt, .1, .9),                       # #
                   smallplot = c(.11 + delt, .14 + delt, .1, .9),            # #
                   legend.lab = "asinh(Counts)", legend.cex = 6,             # #
                   legend.line = 15, axis.args = list(at = legticks,         # #
                                                      labels = leglabs,      # #
                                                      col = "white"))        # #
        # Arrow Field                                                        # #
        vectorField(temp_X, temp_P, temp_x, temp_y, scale = .5, headspan = 0,# #
                    vecspec = "deg", col = "green")                          # #
        # Center Points                                                      # #
        text(tgt_xy[1] / dimsObj[1], tgt_xy[2] / dimsObj[2], "+", cex = 4,   # #
             col = "cyan")                                                   # #
        text(tgt_x / dimsObj[1], tgt_y / dimsObj[2], "x", col = "yellow",    # #
             cex = 4)                                                        # #
        # Main Title                                                         # #
        mtext(paste0("Polarization Map, ", out_label), cex = 4, line = 2,    # #
              col = "white")                                                 # #
        # X-axis                                                             # #
        axis(1, at = xticks, labels = paste0(xlabs, "°"), line = 3, lwd = 0) # #
        axis(1, at = xticks, labels = rep("", length(xticks)), lwd = 2,      # #
             lwd.ticks = 1)                                                  # #
        mtext("Ra (deg)", side = 1, cex = 6, line = 11, col = "white")       # #
        # Y-axis                                                             # #
        axis(2, at = yticks, labels = paste0(ylabs, "°"), lwd = 2)           # #
        mtext("Dec (deg)", side = 2, cex = 6, line = 8, col = "white")       # #
        # N/E referencial                                                    # #
        arrows(x0 = delt+.1, y0 = 0.05, x1 = delt+.1, y1 = 0.1, length =0.25,# #
               col = "red", lwd = 5)                                         # #
        arrows(x0 = delt+.1, y0 = 0.05, x1 = delt+.1 - 0.05, y1 = 0.05,      # #
               length = 0.25, col = "red", lwd = 5)                          # #
        text(x = delt+.1, y = 0.105, labels = "N", pos = 3, col = "red",     # #
             cex = 3)                                                        # #
        text(x = delt+.1-0.053, y = 0.05, labels = "E", pos = 2, col = "red",# #
             cex = 3)                                                        # #
        # Statistical info                                                   # #
        mtext(paste0("min P_gal = ", min_P, "%, median P_gal = ", med_P,     # #
                     "% +- ", unc_med_P, "%, max P_gal = ", max_P, "%"),     # #
              side = 1, line = -8, cex = 2, col = "white")                   # #
        mtext(paste0("min X_gal = ", min_X, "º, median X_gal = ", med_X,     # #
                     "º +- ", unc_med_X, "º, max X_gal = ", max_X, "º"),     # #
              side = 1, line = -6, cex = 2, col = "white")                   # #
                                                                             # #
        # Main Plot                                                          # #
        image.plot(bkg, xlab = ' ', ylab = ' ', xaxt = "n", yaxt = "n",      # #
                   col = gray.colors(64, 0, 1, 0.6),                         # #
                   bigplot = c(.1, .1 + delt, .1, .9),                       # #
                   smallplot = c(.11 + delt, .14 + delt, .1, .9),            # #
                   legend.lab = "asinh(Counts)", legend.cex = 6,             # #
                   legend.line = 15, axis.args = list(at = legticks,         # #
                                                      labels = leglabs,      # #
                                                      col = "white"))        # #
        # Arrow Field                                                        # #
        vectorField(temp_X, temp_IP, temp_x, temp_y, scale = .5,             # #
                    headspan = 0, vecspec = "deg", col = "green")            # #
        # Center Points                                                      # #
        points(tgt_xy[1] / dimsObj[1], tgt_xy[2] / dimsObj[2], pch = 3,      # #
               cex = 4, col = "cyan")                                        # #
        points(tgt_x / dimsObj[1], tgt_y / dimsObj[2], pch = 4, cex = 4,     # #
               col = "yellow")                                               # #
        # Main Title                                                         # #
        mtext(paste0("Polarized Flux Map, ", out_label), cex = 4, line = 2,  # #
              col = "white")                                                 # #
        # X-axis                                                             # #
        axis(1, at = xticks, labels = paste0(xlabs, "°"), line = 3, lwd = 0) # #
        axis(1, at = xticks, labels = rep("", length(xticks)), lwd = 2,      # #
             lwd.ticks = 1)                                                  # #
        mtext("Ra (deg)", side = 1, cex = 6, line = 11, col = "white")       # #
        # Y-axis                                                             # #
        axis(2, at = yticks, labels = paste0(ylabs, "°"), lwd = 2)           # #
        mtext("Dec (deg)", side = 2, cex = 6, line = 8, col = "white")       # #
        # N/E referencial                                                    # #
        arrows(x0 = delt+.1, y0 = 0.05, x1 = delt+.1, y1 = 0.1, length =0.25,# #
               col = "red", lwd = 5)                                         # #
        arrows(x0 = delt+.1, y0 = 0.05, x1 = delt+.1 - 0.05, y1 = 0.05,      # #
               length = 0.25, col = "red", lwd = 5)                          # #
        text(x = delt+.1, y = 0.105, labels = "N", pos = 3, col = "red",     # #
             cex = 3)                                                        # #
        text(x = delt+.1-0.053, y = 0.05, labels = "E", pos = 2, col = "red",# #
             cex = 3)                                                        # #
        # Statistical info                                                   # #
        mtext(paste0("min IP_gal = ", min_IP, ", median IP_gal = ", med_IP,  # #
                     " +- ", unc_med_IP, " max IP_gal = ", max_IP),          # #
              side = 1, line = -8, cex = 2, col = "white")                   # #
        mtext(paste0("min X_gal = ", min_X, "º, median X_gal = ", med_X,     # #
                     "º +- ", unc_med_X, "º, max X_gal = ", max_X, "º"),     # #
              side = 1, line = -6, cex = 2, col = "white")                   # #
        dev.off()                                                            # #
                                                                             # #
        deb_flag <- FALSE                                                    # #
      }                                                                      # #
      rm(temp_P, temp_x, temp_X, temp_y, It, pol_stats, Qt, Ut, X, IP)       # #
                                                                             # #
      print(paste0("** The pdfs for the polarization maps of instance ",     # #
                   pdf_labels[out_ind], " have been created... **"))         # #
    }                                                                        # #
    rm(img_name, file_name, file_type, out_label, pdf_name, out_head, xs, ys,# #
       pdf_list, rowi, coli, binrows, bincols, NAs_ind, out_ind, Pt, Xt, I)  # #
    ############################################################################
  }                                                                            #
  rm(n, Qinst, Uinst, crpix, crval, ctype, fits_name, band_ind_E, band_ind_O,  #
     numb_ind_E, numb_ind_O, crpix2, crval2, ctype2, bQinst, bUinst, bin_dims, #
     binDim)                                                                   #
  print(paste0("### Finished processing beams of band ", bs, " of ", og_obj,   #
               " to bin size ", binS, "x", binS, " ###"))                      #
}                                                                              #
                                                                               #
if(!apt_R_flag){                                                               #
  write.table(apt_rad_arr, paste0(mwSrcPOL_path, "optimal_apt_radius_tab.csv"),#
              sep = ";", dec = ".")                                            #
  rm(apt_rad_arr)                                                              #
}                                                                              #
                                                                               #
print(paste0("##### All beam files of ", og_obj, " have been processed to bin",#
             " size ", binS, "x", binS, " #####"))                             #
                                                                               #
rm(DATE, TARGET, band_str, pdf_labels, file_type_1, file_label, folder_t,      #
   pnt_src_t, ext_src_t, src_sig, fileList, pathList, instList, b, inst_ind,   #
   dimsObj, E_ind, O_ind, E_unc_ind, O_unc_ind, ang_str, numbs_arr, sat, max_A,#
   max_ecc, unc_numbs_arr, run_date, singDim, singNames, multDim, multNames,   #
   bands, skyPol_path, angs, ind_0_E, ind_22_E, ind_45_E, ind_68_E, ind_0_O,   #
   ind_22_O, ind_45_O, ind_68_O, max_N, mwPol_path, REFDEC, REFRA, REFX, REFY, #
   TYPEX, TYPEY)                                                               #
################################################################################