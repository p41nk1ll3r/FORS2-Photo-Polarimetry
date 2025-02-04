rm(list=ls())
require(FITSio)
require(stringr)
require(jjb)
require(colorspace)
require(ggplot2)
require(plot3D)
require(fields)
require(plotrix)
require(dipsaus)
require(reticulate)

og_obj <- "NGC-3351"
print("Provide the name of the target: ")
og_obj <- edit(og_obj)

binS <- 5
print("Define a valid bin size (integer > 0) for the maps to be produced:")
binS <- edit(binS)

while(binS < 1 || binS %% 1 != 0){
  print("Invalid input. Please input a positive integer.")
  binS <- edit(binS)
}

print(paste0("Processing OFFSETS of target ", og_obj, " with Binning of ", binS,
             "x", binS, " pix²."))

home_folder <- "/media/joaomfras/GalaxyPol/Pol-Gal/"
libs_folder <- paste0(home_folder, "Polarimetric-Imaging-Reduction-Scripts-(FO",
                      "RS2)/Commit/")
lib_head_path <- paste0(libs_folder, "fetch_fits_header_info.R")
lib_pro_path <- paste0(libs_folder, "process_lib.R")
source(lib_head_path)
source(lib_pro_path)
rm(lib_head_path, lib_pro_path)

py_run_string("import numpy as np")
py_run_string("from astropy.io import fits")

################################################################################
################################# User Prompts #################################
################################################################################
og_folder <- paste0(home_folder, og_obj, "/Processed_OBJ/")                    #
print("Provide an input path: ")                                               #
og_folder <- edit(og_folder)                                                   #
                                                                               #
astrometry_folder <- paste0(home_folder, og_obj, "/Astrometry/")               #
if(!dir.exists(astrometry_folder)){                                            #
  mkdir(astrometry_folder)                                                     #
}                                                                              #
astrometry_output <- paste0(astrometry_folder, "Results/")                     #
if(!dir.exists(astrometry_output)){                                            #
  mkdir(astrometry_output)                                                     #
}                                                                              #
astrometry_input <- paste0(astrometry_folder, "Input/")                        #
if(!dir.exists(astrometry_input)){                                             #
  mkdir(astrometry_input)                                                      #
}                                                                              #
                                                                               #
astro_cfg_path <- "/etc/astrometry.cfg"                                        #
print("Provide a path for the astrometry .cfg file: ")                         #
astro_cfg_path <- edit(astro_cfg_path)                                         #
                                                                               #
mask_folder <- paste0(home_folder, og_obj, "/Merged_Masks/")                   #
if(!dir.exists(mask_folder)){                                                  #
  mkdir(mask_folder)                                                           #
}                                                                              #
                                                                               #
star_mask_folder <- paste0(mask_folder, "Star_Masks/")                         #
if(!dir.exists(star_mask_folder)){                                             #
  mkdir(star_mask_folder)                                                      #
}                                                                              #
gal_mask_folder <- paste0(mask_folder, "Gal_Masks/")                           #
if(!dir.exists(gal_mask_folder)){                                              #
  mkdir(gal_mask_folder)                                                       #
}                                                                              #
na_mask_folder <- paste0(mask_folder, "NA_Masks/")                             #
if(!dir.exists(na_mask_folder)){                                               #
  mkdir(na_mask_folder)                                                        #
}                                                                              #
                                                                               #
bin_folder <- paste0("bin", binS, "/")                                         #
                                                                               #
stokes_folder <- paste0(og_folder, "Processed_Offsets/Stokes/")                #
stopifnot(dir.exists(stokes_folder))                                           #
                                                                               #
input_folder <- paste0(stokes_folder, bin_folder)                              #
stopifnot(dir.exists(input_folder))                                            #
                                                                               #
obs_Iflux_folder <- paste0(og_folder, "Processed_Offsets/I_flux/Obs/")         #
stopifnot(dir.exists(obs_Iflux_folder))                                        #
                                                                               #
sky_Iflux_folder <- paste0(og_folder, "Processed_Offsets/I_flux/Sky/")         #
stopifnot(dir.exists(sky_Iflux_folder))                                        #
                                                                               #
noSky_Iflux_folder <- paste0(og_folder, "Processed_Offsets/I_flux/Sky_C/")     #
stopifnot(dir.exists(noSky_Iflux_folder))                                      #
                                                                               #
obspol_folder <- paste0(input_folder, "Obs/FITS/")                             #
stopifnot(dir.exists(obspol_folder))                                           #
                                                                               #
sky_instCpol_folder <- paste0(input_folder, "Sky&Inst_C/FITS/")                #
stopifnot(dir.exists(sky_instCpol_folder))                                     #
                                                                               #
skypol_folder <- paste0(input_folder, "Sky/FITS/")                             #
stopifnot(dir.exists(skypol_folder))                                           #
                                                                               #
mwpol_folder <- paste0(stokes_folder, "MW/Lists_of_Sources/Accepted/Pol/")     #
stopifnot(dir.exists(mwpol_folder))                                            #
                                                                               #
mult <- 1.2                                                                    #
                                                                               #
output_folder <- paste0(og_folder, "Merged_Offsets_SMAx", mult, "/")           #
if(!dir.exists(output_folder)){                                                #
  mkdir(output_folder)                                                         #
}                                                                              #
                                                                               #
outstokes_path <- paste0(output_folder, "Stokes/")                             #
if(!dir.exists(outstokes_path)){                                               #
  mkdir(outstokes_path)                                                        #
}                                                                              #
outpol_path <- paste0(outstokes_path, bin_folder)                              #
if(!dir.exists(outpol_path)){                                                  #
  mkdir(outpol_path)                                                           #
}                                                                              #
flux_path_obs <- paste0(output_folder, "I_flux/Obs/")                          #
if(!dir.exists(flux_path_obs)){                                                #
  mkdir(flux_path_obs)                                                         #
}                                                                              #
flux_path_sky <- paste0(output_folder, "I_flux/Sky/")                          #
if(!dir.exists(flux_path_sky)){                                                #
  mkdir(flux_path_sky)                                                         #
}                                                                              #
flux_path_noSky <- paste0(output_folder, "I_flux/Sky_C/")                      #
if(!dir.exists(flux_path_noSky)){                                              #
  mkdir(flux_path_noSky)                                                       #
}                                                                              #
                                                                               #
# Path for files without any polarization correction                           #
obsPol_path <- paste0(outpol_path, "Obs/")                                     #
if(!dir.exists(obsPol_path)){                                                  #
  mkdir(obsPol_path)                                                           #
}                                                                              #
                                                                               #
# Path for files related to Milky Way polarization                             #
mwPol_path <- paste0(outstokes_path, "MW/")                                    #
if(!dir.exists(mwPol_path)){                                                   #
  mkdir(mwPol_path)                                                            #
}                                                                              #
                                                                               #
# Path for files related to MW stars polarization                              #
MWstars_folder <- paste0(home_folder, "Milky_Way_stars/")                      #
print("Provide an output path for Milky Way star polarimetry: ")               #
MWstars_folder <- edit(MWstars_folder)                                         #
if(!dir.exists(MWstars_folder)){                                               #
  mkdir(MWstars_folder)                                                        #
}                                                                              #
                                                                               #
# Path for file related to stars photometric aperture radius                   #
starsApt_R_path <- paste0(mwpol_folder, "optimal_apt_radius_tab.csv")          #
if(!file.exists(starsApt_R_path)){                                             #
  print(paste0("The path provided for the MW stars photometric aperture radiu",#
               "s is not valid. Returning NULL."))                             #
  return(NULL)                                                                 #
}                                                                              #
                                                                               #
# Path for files with all polarization corrections                             #
allCPol_path <- paste0(outpol_path, "All_C/")                                  #
if(!dir.exists(allCPol_path)){                                                 #
  mkdir(allCPol_path)                                                          #
}                                                                              #
                                                                               #
# Creating all Output Folder Tree for files with all polarization corrections  #
allC_versions <- c("skyCut", "skyCut&starsRM")                                 #
allC_outs <- c("Dense", "Sparse", "MxUnc")                                     #
                                                                               #
allC_Vfolders <- c("Sky_Cut/", "Sky_Cut-&-Stars_RM/")                          #
allC_Ofolders <- paste0(allC_outs, "/")                                        #
                                                                               #
allC_paths_dim <- c(length(allC_versions), length(allC_outs))                  #
allC_paths_dname <- list(allC_versions, allC_outs)                             #
                                                                               #
allCPol_paths <- array("", dim = allC_paths_dim, dimnames = allC_paths_dname)  #
rm(allC_paths_dim, allC_paths_dname)                                           #
                                                                               #
for(f in 1:length(allC_versions)){                                             #
  VF <- allC_Vfolders[f]                                                       #
                                                                               #
  for(sf in 1:length(allC_outs)){                                              #
    OF <- allC_Ofolders[sf]                                                    #
                                                                               #
    allCPol_paths[f, sf] <- paste0(allCPol_path, VF, OF)                       #
    if(!dir.exists(allCPol_paths[f, sf] )){                                    #
      mkdir(allCPol_paths[f, sf] )                                             #
    }                                                                          #
  }                                                                            #
}                                                                              #
rm(allC_versions, allC_Vfolders, allC_outs, allC_Ofolders, VF, OF, f, sf)      #
                                                                               #
                                                                               #
# Path for files related to background polarization                            #
skyPol_path <- paste0(outpol_path, "Sky/")                                     #
if(!dir.exists(skyPol_path)){                                                  #
  mkdir(skyPol_path)                                                           #
}                                                                              #
################################################################################

################################################################################
#-_----------------------------------_---_-----------------------------------_-#
#/ \--------------------------------/ \|/ \---------------------------------/ \#
#\_/--------------------------------\_/|\_/---------------------------------\_/#
#------------------------------------------------------------------------------#
################################################################################

################################################################################
########### Setting up files lists, confirming bands and HWP angles ############
################################################################################
# Keywords of interest in fits files' headers                                  #
REFX <- "CRPIX1"                                                               #
REFY <- "CRPIX2"                                                               #
REFRA <- "RA     "                                                             #
REFDEC <- "DEC    "                                                            #
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
BINS <- "BINNING ="                                                            #
PIXSCALE <- "HIERARCH ESO INS PIXSCALE"                                        #
BINX <- "HIERARCH ESO DET WIN1 BINX"                                           #
BINY <- "HIERARCH ESO DET WIN1 BINY"                                           #
EXPT <- "EXPTIME"                                                              #
                                                                               #
# List of output pdf labels                                                    #
pdf_labels <- c("obs", "all C", "sky")                                         #
                                                                               #
# List of input FITS file keywords                                             #
file_type_in <- c('OBSERV', 'SKY_INST C', 'SKY')                               #
                                                                               #
# List of output FITS file keywords                                            #
file_type_out <- c('OBSERV', 'ALL C', 'SKY')                                   #
                                                                               #
# List of output file labels                                                   #
file_label <- c('obs', 'allC', 'sky')                                          #
                                                                               #
# List of output paths                                                         #
folder_t <- c(obsPol_path, allCPol_path, skyPol_path)                          #
rm(obsPol_path, allCPol_path, skyPol_path)                                     #
                                                                               #
band_str <- c("u_HIGH", "b_HIGH", "v_HIGH", "R_SPECIAL", "I_BESS")             #
bands <- 1:length(band_str)                                                    #
nB <- length(bands)                                                            #
dimsObj <- NULL                                                                #
################################################################################

################################################################################
#-_----------------------------------_---_-----------------------------------_-#
#/ \--------------------------------/ \|/ \---------------------------------/ \#
#\_/--------------------------------\_/|\_/---------------------------------\_/#
#------------------------------------------------------------------------------#
################################################################################

################################################################################
############################ I FLUX FILE PROCESSOR #############################
################################################################################
obs_pathList <- list.files(obs_Iflux_folder, full.names = T)                   #
obs_fileList <- list.files(obs_Iflux_folder, full.names = F)                   #
sky_pathList <- list.files(sky_Iflux_folder, full.names = T)                   #
sky_fileList <- list.files(sky_Iflux_folder, full.names = F)                   #
noSky_pathList <- list.files(noSky_Iflux_folder, full.names = T)               #
noSky_fileList <- list.files(noSky_Iflux_folder, full.names = F)               #
                                                                               #
print("### Merging all available offsets of targets I flux maps ###")          #
                                                                               #
Ipaths_labels <- c("obs", "skyC", "sky")                                       #
Ipaths_dim <- c(length(Ipaths_labels), nB)                                     #
Ipaths_dname <- list(Ipaths_labels, band_str)                                  #
label_tag <- "skyC"                                                            #
                                                                               #
Stacked_Ipaths <- array("", dim = Ipaths_dim, dimnames = Ipaths_dname)         #
Merged_Ipaths <- array("", dim = Ipaths_dim, dimnames = Ipaths_dname)          #
Astro_Ipaths <- array("", dim = nB, dimnames = list(band_str))                 #
na_mask_paths <- array("", dim = nB, dimnames = list(band_str))                #
                                                                               #
ref_pix_par <- c('x', 'y', 'ra', 'dec')                                        #
ref_mrg_par <- c("x", "y", "ra", "dec", "scl_x", "scl_y", "scl_flg")           #
munc_par <- c("Data", "Unc")                                                   #
                                                                               #
I_dname <- list(NULL, NULL, munc_par)                                          #
                                                                               #
for(b in bands){                                                               #
  bs <- band_str[b]                                                            #
                                                                               #
  obs_b_ind <- grep(bs, obs_fileList)                                          #
  sky_b_ind <- grep(bs, sky_fileList)                                          #
  noSky_b_ind <- grep(bs, noSky_fileList)                                      #
                                                                               #
  ifluxN <- length(noSky_b_ind)                                                #
                                                                               #
  if(length(obs_b_ind) != 0){                                                  #
    obs_tpathList <- obs_pathList[obs_b_ind]                                   #
    obs_tfileList <- obs_fileList[obs_b_ind]                                   #
  }                                                                            #
  if(length(sky_b_ind) != 0){                                                  #
    sky_tpathList <- sky_pathList[sky_b_ind]                                   #
    sky_tfileList <- sky_fileList[sky_b_ind]                                   #
  }                                                                            #
  if(length(noSky_b_ind) != 0){                                                #
    noSky_tpathList <- noSky_pathList[noSky_b_ind]                             #
    noSky_tfileList <- noSky_fileList[noSky_b_ind]                             #
  }                                                                            #
                                                                               #
  if(length(obs_b_ind) == 0 & length(sky_b_ind) == 0 &                         #
     length(noSky_b_ind) == 0){                                                #
    next                                                                       #
  }                                                                            #
                                                                               #
  ref_pix_dim <- c(ifluxN, length(ref_pix_par))                                #
  ref_pix_dname <- list(NULL, ref_pix_par)                                     #
                                                                               #
  ref_pixs <- array(NA, dim = ref_pix_dim, dimnames = ref_pix_dname)           #
                                                                               #
  if(is.null(dimsObj)){                                                        #
    dimsObj <- dim(readFITS(noSky_tpathList[1])$imDat)[1:2]                    #
    I_dim <- c(dimsObj, length(munc_par))                                      #
  }                                                                            #
                                                                               #
  basename <- strsplit(obs_tfileList, "_#")[[1]][1]                            #
  stackname <- paste0(basename, "-Stacked_I-Flux_")                            #
  mergename <- paste0(basename, "-Merged_I-Flux_")                             #
                                                                               #
  Astro_Ipaths[b] <- paste0(astrometry_input, basename, "_Astro_I-Flux.fits")  #
                                                                               #
  skyCind <- grep("I_noSky_Cut", noSky_tfileList)                              #
  AskyCind <- grep("I_noSky_noCut", noSky_tfileList)                           #
                                                                               #
  Stacked_Ipaths["obs", b] <- paste0(flux_path_obs, stackname, "obs.fits")     #
  Stacked_Ipaths["skyC", b] <- paste0(flux_path_noSky, stackname, "skyC.fits") #
  Stacked_Ipaths["sky", b] <- paste0(flux_path_sky, stackname, "sky.fits")     #
                                                                               #
  Merged_Ipaths["obs", b] <- paste0(flux_path_obs, mergename, "obs.fits")      #
  Merged_Ipaths["skyC", b] <- paste0(flux_path_noSky, mergename, "skyC.fits")  #
  Merged_Ipaths["sky", b] <- paste0(flux_path_sky, mergename, "sky.fits")      #
                                                                               #
  for(type in 1:length(Ipaths_labels)){                                        #
    switch(type,                                                               #
           c(ttpathList <- obs_tpathList, ttfileList <- obs_tfileList),        #
           c(ttpathList <- noSky_tpathList[skyCind],                           #
             ttfileList <- noSky_tfileList[skyCind],                           #
             atpathList <- noSky_tpathList[AskyCind]),                         #
           c(ttpathList <- sky_tpathList, ttfileList <- sky_tfileList)         #
    )                                                                          #
                                                                               #
    this_type <- Ipaths_labels[type]                                           #
                                                                               #
    N <- length(ttfileList)                                                    #
                                                                               #
    # Checking for results of previous runs                                    #
    stackflag <- file.exists(Stacked_Ipaths[type, b])                          #
    mergeflag <- file.exists(Merged_Ipaths[type, b])                           #
                                                                               #
    switch(type,                                                               #
           astroflag <- 1, astroflag <- file.exists(Astro_Ipaths[b]),          #
           astroflag <- 1)                                                     #
                                                                               #
    if(N != 1){                                                                #
      print(paste0("* Merging I flux offsets in ", bs, " band for ",           #
                   this_type, " case... *"))                                   #
                                                                               #
      A_dim <- c(dimsObj, N)                                                   #
      tmp_dim <- c(dimsObj, length(munc_par), N)                               #
      tmp_dname <- list(NULL, NULL, munc_par, NULL)                            #
      ref_mrg_dim <- c(length(ref_mrg_par), N)                                 #
      ref_mrg_dname <- list(ref_mrg_par, NULL)                                 #
                                                                               #
      temp_arr <- array(NA, dim = tmp_dim, dimnames = tmp_dname)               #
                                                                               #
      if(this_type == label_tag){                                              #
        temp_ast <- array(NA, dim = A_dim)                                     #
        rm(A_dim)                                                              #
      }                                                                        #
                                                                               #
      for(n in 1:N){                                                           #
        temp_fits <- readFITS(ttpathList[n])                                   #
        temp_data <- temp_fits$imDat                                           #
        temp_head <- temp_fits$header[-c(8:10)]                                #
        rm(temp_fits)                                                          #
                                                                               #
        if(this_type == label_tag){                                            #
          ast_fits <- readFITS(atpathList[n])                                  #
          ast_data <- ast_fits$imDat                                           #
          ast_head <- ast_fits$header[-c(8:10)]                                #
          rm(ast_fits)                                                         #
        }                                                                      #
                                                                               #
        refx <- get_fits_header_num(temp_head, REFX)                           #
        refy <- get_fits_header_num(temp_head, REFY)                           #
        refz <- 1                                                              #
        refz2 <- 1                                                             #
        valx <- get_fits_header_num(temp_head, VALX)                           #
        valy <- get_fits_header_num(temp_head, VALY)                           #
        valz <- NA                                                             #
        valz2 <- NA                                                            #
        typex <- "  'PIXEL     '              / Coordinate system of x-axis "  #
        typey <- "  'PIXEL     '              / Coordinate system of y-axis "  #
        typez <- "  'DATA & UNC'              / Coordinate system of z1-axis " #
        typez2 <- "  'OFFSET    '             / Coordinate system of z2-axis " #
        crpix <- c(refx, refy, refz, refz2)                                    #
        crval <- c(valx, valy, valz, valz2)                                    #
        ctype <- c(typex, typey, typez, typez2)                                #
        binx <- get_fits_header_num(temp_head, BINX)                           #
        biny <- get_fits_header_num(temp_head, BINY)                           #
        pxscl <- get_fits_header_num(temp_head, PIXSCALE)                      #
        exptime <- get_fits_header_num(temp_head, EXPT)                        #
        rm(refx, refy, refz, refz2, valx, valy, valz, valz2,                   #
           typex, typey, typez, typez2)                                        #
                                                                               #
        if(n == 1){                                                            #
          ref_mrg <- array(NA, dim = ref_mrg_dim, dimnames = ref_mrg_dname)    #
          hold_crpix <- crpix                                                  #
          hold_crval <- crval                                                  #
          hold_ctype <- ctype                                                  #
          hold_coords <- c(get_fits_header_num(temp_head, REFRA),              #
                           get_fits_header_num(temp_head, REFDEC))             #
        }                                                                      #
                                                                               #
        ref_mrg["x", n] <- crpix[1]                                            #
        ref_mrg["y", n] <- crpix[2]                                            #
        ref_mrg["ra", n] <- get_fits_header_num(temp_head, REFRA)              #
        ref_mrg["dec", n] <- get_fits_header_num(temp_head, REFDEC)            #
        ref_mrg["scl_x", n] <- pxscl * binx                                    #
        ref_mrg["scl_y", n] <- pxscl * biny                                    #
        ref_mrg["scl_flg", n] <- 1                                             #
                                                                               #
        if(this_type == label_tag){                                            #
          ref_pixs[n, c('x', 'y')] <- crpix[1:2]                               #
          ref_pixs[n, c('ra', 'dec')] <- ref_mrg[c("ra","dec"), n]             #
        }                                                                      #
                                                                               #
        if(n == 1){                                                            #
          scl_x <- ref_mrg["scl_x", 1] / 3600                                  #
          scl_y <- ref_mrg["scl_y", 1] / 3600                                  #
          this_head <- temp_head                                               #
          rm(temp_head)                                                        #
                                                                               #
          if(this_type == label_tag){                                          #
            astro_head <- ast_head                                             #
            rm(ast_head)                                                       #
          }                                                                    #
        }                                                                      #
                                                                               #
        if(n != 1 && (ref_mrg["scl_x", n] != ref_mrg["scl_x", 1] ||            #
                      ref_mrg["scl_y", n] != ref_mrg["scl_y", 1])){            #
          ref_mrg["scl_flg", n] <- 0                                           #
                                                                               #
          print(paste0("WARNING: Offset #", n, " has a different scale than t",#
                       "hat of Offset #1, and as such it will not merged with",#
                       " the remaining offsets."))                             #
        }                                                                      #
        rm(binx, biny, pxscl)                                                  #
                                                                               #
        # Skipping calculation in case files already exist                     #
        if(stackflag & mergeflag & astroflag){                                 #
        }else{                                                                 #
          if(n == 1){                                                          #
            temp_arr[,,"Data", n] <- temp_data[,,1] / exptime                  #
            temp_arr[,,"Unc", n] <- unc_div(temp_data[,,1], exptime,           #
                                            temp_data[,,2], 0.0001,            #
                                            temp_arr[,,"Data", n])             #
                                                                               #
                                                                               #
            if(this_type == label_tag){                                        #
              temp_ast[,,n] <- ast_data[,,1] / exptime                         #
            }                                                                  #
          }else{                                                               #
            if(ref_mrg["scl_flg", n]){                                         #
              dx <- ref_mrg["x", 1] - ref_mrg["x", n]                          #
              dy <- ref_mrg["y", 1] - ref_mrg["y", n]                          #
                                                                               #
              cosd <- cos(mean(c(ref_mrg["dec", 1], ref_mrg["dec", n]))        #
                          * pi / 180)                                          #
                                                                               #
              dra_x <- (ref_mrg["ra", 1] - ref_mrg["ra", n]) * cosd / scl_x    #
              ddec_y <- (ref_mrg["dec", 1] - ref_mrg["dec", n]) / scl_y        #
              difx <- round(dx + dra_x)                                        #
              dify <- round(dy + ddec_y)                                       #
              rm(dx, dy, dra_x, ddec_y)                                        #
                                                                               #
              if(difx >= 0){                                                   #
                hold_xs <- 1:(dimsObj[1] - difx)                               #
                off_xs <- (1 + difx):dimsObj[1]                                #
              }else{                                                           #
                hold_xs <- (1 - difx):dimsObj[1]                               #
                off_xs <- 1:(dimsObj[1] + difx)                                #
              }                                                                #
              if(dify <= 0){                                                   #
                hold_ys <- (1 - dify):dimsObj[2]                               #
                off_ys <- 1:(dimsObj[2] + dify)                                #
              }else{                                                           #
                hold_ys <- 1:(dimsObj[2] - dify)                               #
                off_ys <- (1 + dify):dimsObj[2]                                #
              }                                                                #
              rm(difx, dify)                                                   #
            }                                                                  #
                                                                               #
            temp_arr[hold_xs, hold_ys, "Data", n] <-                           #
              temp_data[off_xs, off_ys, 1] / exptime                           #
            temp_arr[hold_xs, hold_ys, "Unc", n] <-                            #
              unc_div(temp_data[off_xs, off_ys, 1], exptime,                   #
                      temp_data[off_xs, off_ys, 2], 0.0001,                    #
                      temp_arr[hold_xs, hold_ys, "Data", n])                   #
                                                                               #
            if(this_type == label_tag){                                        #
              temp_ast[hold_xs, hold_ys, n] <-                                 #
                ast_data[off_xs, off_ys, 1] / exptime                          #
            }                                                                  #
          }                                                                    #
        }                                                                      #
      }                                                                        #
                                                                               #
      # Skipping file creation if they already exist                           #
      if(stackflag & mergeflag & astroflag){                                   #
      }else{                                                                   #
        thisI <- array(NA, dim = I_dim, dimnames = I_dname)                    #
                                                                               #
        thisI[,,1] <- apply(temp_arr[,,"Data",], 1:2, median, na.rm = T)       #
        thisI[,,2] <- unc_median(temp_arr[,,"Data",], temp_arr[,,"Unc",], 1:2) #
                                                                               #
        if(this_type == label_tag){                                            #
          print(paste0("* Building Astrometry input file for ", bs,            #
                       " band... *"))                                          #
          astroI <- apply(temp_ast, 1:2, median, na.rm = T)                    #
          rm(temp_ast)                                                         #
                                                                               #
          for(y in 400:1975){                                                  #
            for(x in  190:1850){                                               #
              if(is.na(astroI[x,y])){                                          #
                deltit <- 7                                                    #
                t_samp <- NULL                                                 #
                                                                               #
                while(length(t_samp) == 0){                                    #
                  xi <- (-deltit):deltit + x                                   #
                  yi <- (-deltit):deltit + y                                   #
                  t_samp <- as.vector(astroI[xi, yi])                          #
                  t_samp <- t_samp[which(!is.na(t_samp), arr.ind = T)]         #
                                                                               #
                  deltit <- deltit + 2                                         #
                }                                                              #
                astroI[x,y] <- sample(t_samp, 1)                               #
              }                                                                #
            }                                                                  #
          }                                                                    #
          rm(x, y, t_samp, xi, yi, deltit)                                     #
        }                                                                      #
                                                                               #
        switch(type,                                                           #
               c(writeFITSim(temp_arr, Stacked_Ipaths["obs", b],               #
                             crpixn = hold_crpix, crvaln = hold_crval,         #
                             ctypen = hold_ctype, header = this_head),         #
                 writeFITSim(thisI, Merged_Ipaths["obs", b],                   #
                             crpixn = hold_crpix[1:3], crvaln =hold_crval[1:3],#
                             ctypen = hold_ctype[1:3], header = this_head)),   #
               c(writeFITSim(temp_arr, Stacked_Ipaths["skyC", b],              #
                             crpixn = hold_crpix, crvaln = hold_crval,         #
                             ctypen = hold_ctype, header = this_head),         #
                 writeFITSim(thisI, Merged_Ipaths["skyC", b],                  #
                             crpixn = hold_crpix[1:3], crvaln =hold_crval[1:3],#
                             ctypen = hold_ctype[1:3], header = this_head),    #
                 writeFITSim(astroI, Astro_Ipaths[b], crpixn = hold_crpix[1:2],#
                             crvaln = hold_crval[1:2], ctypen =hold_ctype[1:2],#
                             header = astro_head)),                            #
               c(writeFITSim(temp_arr, Stacked_Ipaths["sky", b],               #
                             crpixn = hold_crpix, crvaln = hold_crval,         #
                             ctypen = hold_ctype, header = this_head),         #
                 writeFITSim(thisI, Merged_Ipaths["sky", b],                   #
                             crpixn = hold_crpix[1:3], crvaln =hold_crval[1:3],#
                             ctypen = hold_ctype[1:3], header = this_head))    #
        )                                                                      #
      }                                                                        #
    }else{                                                                     #
      print(paste0("* There is only 1 I flux offset in ", bs, " band for ",    #
                   this_type, " case. Creating a normalized copy to output fo",#
                   "lder. *"))                                                 #
                                                                               #
      temp_arr <- array(NA, dim = I_dim, dimnames = I_dname)                   #
                                                                               #
      temp_fits <- readFITS(ttpathList[1])                                     #
      temp_arr[,,"Data"] <- temp_fits$imDat[,,1]                               #
      temp_arr[,,"Unc"] <- temp_fits$imDat[,,2]                                #
      this_head <- temp_fits$header[-c(8:10)]                                  #
      rm(temp_fits)                                                            #
                                                                               #
      if(this_type == label_tag){                                              #
        ast_fits <- readFITS(atpathList[1])                                    #
        astroI <- ast_fits$imDat[,,1]                                          #
        astro_head <- ast_fits$header[-c(8:10)]                                #
        rm(ast_fits)                                                           #
      }                                                                        #
                                                                               #
      refx <- get_fits_header_num(this_head, REFX)                             #
      refy <- get_fits_header_num(this_head, REFY)                             #
      refz <- 1                                                                #
      valx <- get_fits_header_num(this_head, VALX)                             #
      valy <- get_fits_header_num(this_head, VALY)                             #
      valz <- NA                                                               #
      typex <- "  'PIXEL     '              / Coordinate system of x-axis "    #
      typey <- "  'PIXEL     '              / Coordinate system of y-axis "    #
      typez <- "  'DATA & UNC'              / Coordinate system of z-axis "    #
      hold_crpix <- c(refx, refy, refz)                                        #
      hold_crval <- c(valx, valy, valz)                                        #
      hold_ctype <- c(typex, typey, typez)                                     #
                                                                               #
      ref_pixs <- hold_crpix[1:2]                                              #
                                                                               #
      # Skipping calculations if files already exist                           #
      if(stackflag & astroflag){                                               #
      }else{                                                                   #
        if(this_type == label_tag){                                            #
          astro_NAs <- which(is.na(astroI), arr.ind = T)                       #
          aNA_l <- length(astro_NAs) / 2                                       #
                                                                               #
          m_astro <- median(astroI, na.rm = T)                                 #
          s_astro <- mad(astroI, na.rm = T)                                    #
                                                                               #
          astroI[astro_NAs] <- rnorm(aNA_l, m_astro, s_astro)                  #
          rm(aNA_l, astro_NAs, m_astro, s_astro)                               #
        }                                                                      #
                                                                               #
        switch(type,                                                           #
               writeFITSim(temp_arr, Stacked_Ipaths["obs", b],                 #
                           crpixn = hold_crpix, crvaln = hold_crval,           #
                           ctypen = hold_ctype, header = this_head),           #
               c(writeFITSim(temp_arr, Stacked_Ipaths["skyC", b],              #
                           crpixn = hold_crpix, crvaln = hold_crval,           #
                           ctypen = hold_ctype, header = this_head),           #
                 writeFITSim(astroI, Astro_Ipaths[b], crpixn = hold_crpix[1:2],#
                             crvaln =hold_crval[1:2], ctypen = hold_ctype[1:2],#
                             header = astro_head)),                            #
               writeFITSim(temp_arr, Stacked_Ipaths["sky", b],                 #
                           crpixn = hold_crpix, crvaln = hold_crval,           #
                           ctypen = hold_ctype, header = this_head)            #
        )                                                                      #
      }                                                                        #
    }                                                                          #
                                                                               #
    na_mask_paths[b] <- paste0(na_mask_folder, basename,"_NA-Mask.fits")       #
    maskflag <- file.exists(na_mask_paths[b])                                  #
                                                                               #
    # Check if NA mask fits file exists if it doesn't create it                #
    if(maskflag){                                                              #
      print(paste0("* NA mask for ", bs, " band already exists... *"))         #
    }else{                                                                     #
      print(paste0("* Creating NA mask for ", bs, " band... *"))               #
                                                                               #
      if(N != 1){                                                              #
        thisI <- readFITS(Merged_Ipaths["obs", b])$imDat[,,1]                  #
      }else{                                                                   #
        thisI <- readFITS(Stacked_Ipaths["obs", b])$imDat[,,1]                 #
      }                                                                        #
                                                                               #
      na_mask <- array(1, dim = dimsObj)                                       #
      na_mask[which(is.na(thisI), arr.ind = T)] <- NA                          #
                                                                               #
      writeFITSim(na_mask, na_mask_paths[bs])                                  #
      rm(thisI, na_mask)                                                       #
    }                                                                          #
  }                                                                            #
                                                                               #
  switch(b,                                                                    #
         ref_pixs_U <- ref_pixs, ref_pixs_B <- ref_pixs,                       #
         ref_pixs_V <- ref_pixs, ref_pixs_R <- ref_pixs,                       #
         ref_pixs_I <- ref_pixs)                                               #
  rm(ref_pixs)                                                                 #
}                                                                              #
rm(temp_arr, thisI, hold_crpix, hold_crval, hold_ctype, this_head, basename,   #
   astroI, astro_head)                                                         #
################################################################################

################################################################################
#-_----------------------------------_---_-----------------------------------_-#
#/ \--------------------------------/ \|/ \---------------------------------/ \#
#\_/--------------------------------\_/|\_/---------------------------------\_/#
#------------------------------------------------------------------------------#
################################################################################

################################################################################
################### Performing Astrometry on I Flux Images #####################
################################################################################
astroPaths <- array(NA, dim = nB, dimnames = list(band_str))                   #
indexPaths <- array(NA, dim = nB, dimnames = list(band_str))                   #
rdlsPaths <- array(NA, dim = nB, dimnames = list(band_str))                    #
                                                                               #
astroList <- list.files(astrometry_output, full.names = T)                     #
                                                                               #
for(b in bands){                                                               #
  bs <- band_str[b]                                                            #
  iN <- grep(bs, Astro_Ipaths)                                                 #
                                                                               #
  if(length(iN) == 0){                                                         #
    next                                                                       #
  }                                                                            #
                                                                               #
  temp_path <- Astro_Ipaths[iN]                                                #
                                                                               #
  iName <- strsplit(temp_path, "/")[[1]][length(strsplit(temp_path,"/")[[1]])] #
                                                                               #
  print(paste0("Looking for astrometry results for ", iName, "..."))           #
                                                                               #
  astro_name <- paste0(strsplit(iName, ".fits")[[1]], ".new")                  #
  index_name <- paste0(strsplit(iName, ".fits")[[1]], "-indx.xyls")            #
  rdls_name <- paste0(strsplit(iName, ".fits")[[1]], ".rdls")                  #
                                                                               #
  is_astro <- grep(astro_name, astroList)                                      #
  is_index <- grep(index_name, astroList)                                      #
  is_rdls <- grep(rdls_name, astroList)                                        #
                                                                               #
  match_flag <- length(is_astro) * length(is_index) * length(is_rdls)          #
                                                                               #
  if(match_flag == 0){                                                         #
    print(paste0("Complete astrometry results for ", iName, " have not been f",#
                 "ound. Proceeding with creation of astrometry products for t",#
                 "his file to be used later on."))                             #
                                                                               #
    tempI <- readFITS(temp_path)                                               #
                                                                               #
    header <- tempI$header[-c(6, 21:25)]                                       #
    rm(tempI)                                                                  #
                                                                               #
    ref_coords <- c(get_fits_header_num(header, REFRA),                        #
                    get_fits_header_num(header, REFDEC))                       #
    rm(header)                                                                 #
                                                                               #
    temp_I_path <- paste0(astrometry_output, strsplit(iName, ".fits")[[1]])    #
                                                                               #
    thisAstroPath <- paste0(temp_I_path, ".new")                               #
    thisIndexPath <- paste0(temp_I_path, "-indx.xyls")                         #
    thisRdlsPath <- paste0(temp_I_path, ".rdls")                               #
                                                                               #
    astrometry_cmd <- paste0("solve-field --config ", astro_cfg_path, " --ra ",#
                             ref_coords[1], " --dec ", ref_coords[2], " --rad",#
                             "ius .25 --dir ", astrometry_output," ",          #
                             temp_path)                                        #
    system(astrometry_cmd)                                                     #
                                                                               #
    astroPaths[bs] <- thisAstroPath                                            #
    indexPaths[bs] <- thisIndexPath                                            #
    rdlsPaths[bs] <- thisRdlsPath                                              #
    rm(thisAstroPath, thisIndexPath, thisRdlsPath , ref_coords)                #
  }else{                                                                       #
    print(paste0("Astrometry results for ", iName, " have been found. Adding ",#
                 "them to list of files to be used later on."))                #
                                                                               #
    astroPaths[bs] <- astroList[is_astro]                                      #
    indexPaths[bs] <- astroList[is_index]                                      #
    rdlsPaths[bs] <- astroList[is_rdls]                                        #
  }                                                                            #
  rm(is_astro, is_index, is_rdls, match_flag, astro_name, index_name,          #
     rdls_name, iName, temp_path)                                              #
}                                                                              #
rm(astroList, iN, bs)                                                          #
                                                                               #
# Checking if all needed bands have an astrometry file associated to them      #
# If not, associate the astrometry file of the closest band                    #
for(b in bands){                                                               #
  bs <- band_str[b]                                                            #
  iN <- grep(bs, astroPaths)                                                   #
                                                                               #
  if(length(iN) == 0){                                                         #
    next                                                                       #
  }                                                                            #
                                                                               #
  max_ts <- nB - 1                                                             #
  t_bs <- NULL                                                                 #
  ts <- 1                                                                      #
                                                                               #
  while(length(t_bs) != max_ts){                                               #
    if(b == 1){                                                                #
      poss_ts <- b + ts                                                        #
    }                                                                          #
    if(b == nB){                                                               #
      poss_ts <- b - ts                                                        #
    }                                                                          #
    if(b > 1 & b < nB){                                                        #
      poss_ts <- c(b + ts, b - ts)                                             #
    }                                                                          #
                                                                               #
    neg_ts <- which(poss_ts <= 0, arr.ind = T)                                 #
                                                                               #
    if(length(neg_ts) != 0){                                                   #
      poss_ts <- poss_ts[-neg_ts]                                              #
    }                                                                          #
                                                                               #
    out_ts <- which(poss_ts > nB, arr.ind = T)                                 #
                                                                               #
    if(length(out_ts) != 0){                                                   #
      poss_ts <- poss_ts[-out_ts]                                              #
    }                                                                          #
                                                                               #
    rm_ts <- which(!is.na(match(poss_ts, t_bs)), arr.ind = T)                  #
                                                                               #
    if(length(rm_ts) != 0){                                                    #
      poss_ts <- poss_ts[-rm_ts]                                               #
    }                                                                          #
                                                                               #
    t_bs <- c(t_bs, poss_ts)                                                   #
    rm(poss_ts, out_ts, neg_ts, rm_ts)                                         #
                                                                               #
    ts <- ts + 1                                                               #
  }                                                                            #
  rm(max_ts)                                                                   #
                                                                               #
  ts <- 1                                                                      #
  while(!file.exists(astroPaths[bs]) & ts < nB){                               #
    t_f <- t_bs[ts]                                                            #
    astroPaths[bs] <- astroPaths[t_f]                                          #
    indexPaths[bs] <- indexPaths[t_f]                                          #
    rdlsPaths[bs] <- rdlsPaths[t_f]                                            #
    ts <- ts + 1                                                               #
  }                                                                            #
                                                                               #
  if(!file.exists(astroPaths[bs]) & ts == nB){                                 #
    print(paste0("ERROR: No Astrometry file available for this field. Correct",#
                 " the astrometry.net input file (may require changes to this",#
                 "code) and try again. Returning NULL."))                      #
    return(NULL)                                                               #
  }                                                                            #
}                                                                              #
################################################################################

################################################################################
#-_----------------------------------_---_-----------------------------------_-#
#/ \--------------------------------/ \|/ \---------------------------------/ \#
#\_/--------------------------------\_/|\_/---------------------------------\_/#
#------------------------------------------------------------------------------#
################################################################################

################################################################################
################## Creation of Pixel Mask of Central Object ####################
################################################################################
print("### Creating pixel masks for Galaxy in merged offset maps ###")         #
gal_mask_paths <- array("", dim = nB, dimnames = list(band_str))               #
                                                                               #
# Arrays to hold necessary correction parameters in RA,DEC to X,Y conversion   #
cosd <- array(NA, dim = nB, dimnames = list(band_str))                         #
iflux_pix <- array(NA, dim = c(2, nB), dimnames = list(NULL, band_str))        #
iflux_coord <- array(NA, dim = c(2, nB), dimnames = list(NULL, band_str))      #
                                                                               #
# Array to hold coordinates to be operated with now and plotted later          #
pp2p <- array(NA, dim = c(2, 2, nB), dimnames = list(NULL, NULL, band_str))    #
                                                                               #
for(b in bands){                                                               #
  bs <- band_str[b]                                                            #
                                                                               #
  # Load iflux file to extract header info                                     #
  iflux_file <- Astro_Ipaths[bs]                                               #
                                                                               #
  if(iflux_file == ""){                                                        #
    next                                                                       #
  }                                                                            #
                                                                               #
  print(paste0("* Processing files in ", bs, " band *"))                       #
  iflux_fits <- readFITS(iflux_file)                                           #
  iflux_head <- iflux_fits$header                                              #
  rm(iflux_fits)                                                               #
                                                                               #
  # Get reference coords for iflux                                             #
  iflux_refx <- get_fits_header_num(iflux_head, REFX)                          #
  iflux_refy <- get_fits_header_num(iflux_head, REFY)                          #
  iflux_ref_coords <- c(get_fits_header_num(iflux_head, REFRA),                #
                        get_fits_header_num(iflux_head, REFDEC))               #
  binx <- get_fits_header_num(iflux_head, BINX)                                #
  biny <- get_fits_header_num(iflux_head, BINY)                                #
  pxscl <- get_fits_header_num(iflux_head, PIXSCALE)                           #
  ascl <- pxscl * binx                                                         #
  scl_x <- pxscl * binx / 3600                                                 #
  scl_y <- pxscl * biny / 3600                                                 #
  rm(iflux_refx, iflux_refy, pxscl, binx, biny, iflux_head)                    #
                                                                               #
  # Get Coordinates of object within iflux                                     #
  tgt_coords <- get_NED_coords(og_obj, 6)                                      #
  temp_a <- get_NED_median_diam(og_obj)                                        #
  tgt_a <- round(temp_a / (2 * ascl))                                          #
  b_over_a <- get_NED_largest_axis_ratio(og_obj)                               #
  tgt_pa <- (get_NED_newest_posang(og_obj) - 90) / 180 * pi                    #
                                                                               #
  # Get the band b astrometry file                                             #
  astro_file <- astroPaths[bs]                                                 #
                                                                               #
  # Getting Galaxy center pixel coordinates                                    #
  iflux_tgt_xy <- convert_wcs_to_pixel(array(tgt_coords, dim = c(1,2)),        #
                                       astro_file)                             #
  iflux_crpix <- convert_wcs_to_pixel(array(iflux_ref_coords, dim = c(1,2)),   #
                                      astro_file)                              #
                                                                               #
  # Save coordinates to holder arrays                                          #
  iflux_pix[,b] <- iflux_crpix                                                 #
  iflux_coord[,b] <- iflux_ref_coords                                          #
  cosd[b] <- cos(mean(c(tgt_coords[2], iflux_ref_coords[2]) * pi / 180))       #
  rm(iflux_crpix, iflux_ref_coords)                                            #
                                                                               #
  pp2p[1, 1, b] <- iflux_tgt_xy[1,1]                                           #
  pp2p[1, 2, b] <- iflux_tgt_xy[1,2]                                           #
                                                                               #
  t_split <- strsplit(astro_file, "/")                                         #
  basename <- strsplit(t_split[[1]][length(t_split[[1]])], "Astro")[[1]][1]    #
                                                                               #
  gal_mask_paths[b] <- paste0(gal_mask_folder, basename,"_Gal-Mask.fits")      #
  maskflag <- file.exists(gal_mask_paths[b])                                   #
                                                                               #
  if(maskflag){                                                                #
    print(paste0("* Gal mask for ", bs, " band already exists... *"))          #
  }else{                                                                       #
    mask <- array(NA, dim = dimsObj)                                           #
                                                                               #
    if(is.na(b_over_a) | is.na(tgt_pa)){                                       #
      # Loading iflux image and target info to python                          #
      load_cmd <- paste0("ifluxF = fits.open('", iflux_file, "')")             #
      py_run_string(load_cmd)                                                  #
      rm(load_cmd)                                                             #
                                                                               #
      py_run_string("img_dat = ifluxF[0].data")                                #
      py_run_string("ifluxF.close()")                                          #
      py_run_string("X = r.iflux_tgt_xy[0][0]")                                #
      py_run_string("Y = r.iflux_tgt_xy[0][1]")                                #
      py_run_string("SMA = r.tgt_a")                                           #
      rm(iflux_tgt_xy)                                                         #
                                                                               #
      try_a <- round(tgt_a * mult)                                             #
      py_run_string("t_sma = r.try_a")                                         #
                                                                               #
      # Python commands to get shape of target mask                            #
      get_cmd <- paste0("iso = sexD.get_closest_iso_for_target(img_dat, X, Y,",#
                        " t_sma, 0, 90)")                                      #
      py_run_string(get_cmd)                                                   #
      rm(get_cmd)                                                              #
                                                                               #
      try_eps <- py$iso$eps                                                    #
      try_pa <- py$iso$pa                                                      #
      try_a <- py$iso$sma                                                      #
      try_b <- try_a * sqrt(1 - try_eps)                                       #
                                                                               #
      pp2p[2, 1, b] <- py$iso$x0                                               #
      pp2p[2, 2, b] <- py$iso$y0                                               #
      py_run_string("del(img_dat, X, Y, SMA, t_sma, iso)")                     #
                                                                               #
      for(xa in 1:dimsObj[1]){                                                 #
        for(ya in 1:dimsObj[2]){                                               #
          rx <- ((xa - pp2p[2, 1, b]) * cos(try_pa) + (ya - pp2p[2, 2, b]) *   #
                   sin(try_pa))^2 / try_a^2                                    #
          ry <- ((xa - pp2p[2, 1, b]) * sin(try_pa) - (ya - pp2p[2, 2, b]) *   #
                   cos(try_pa))^2 / try_b^2                                    #
                                                                               #
          if(rx + ry <= 1){                                                    #
            mask[xa, ya] <- 1                                                  #
          }                                                                    #
        }                                                                      #
      }                                                                        #
      rm(try_pa, try_b, try_a)                                                 #
    }else{                                                                     #
      tgt_a <- round((temp_a + sqrt(temp_a)) / (2 * ascl))                     #
      tgt_b <- b_over_a * tgt_a                                                #
                                                                               #
      pp2p[2, 1, b] <- iflux_tgt_xy[1,1]                                       #
      pp2p[2, 2, b] <- iflux_tgt_xy[1,2]                                       #
                                                                               #
      for(xa in 1:dimsObj[1]){                                                 #
        for(ya in 1:dimsObj[2]){                                               #
          rx <- ((xa - pp2p[1, 1, b]) * cos(tgt_pa) + (ya - pp2p[1, 2, b]) *   #
                   sin(tgt_pa))^2 / tgt_a^2                                    #
          ry <- ((xa - pp2p[1, 1, b]) * sin(tgt_pa) - (ya - pp2p[1, 2, b]) *   #
                   cos(tgt_pa))^2 / tgt_b^2                                    #
                                                                               #
          if(rx + ry <= 1){                                                    #
            mask[xa, ya] <- 1                                                  #
          }                                                                    #
        }                                                                      #
      }                                                                        #
      rm(tgt_a, tgt_b)                                                         #
    }                                                                          #
                                                                               #
    writeFITSim(mask, gal_mask_paths[bs])                                      #
    rm(xa, ya, rx, ry, mask)                                                   #
  }                                                                            #
  rm(tgt_pa, b_over_a, temp_a)                                                 #
                                                                               #
  # Converting operation coordinates to plot scale coordinates                 #
  pp2p[,1, b] <- pp2p[,1, b] / dimsObj[1]                                      #
  pp2p[,2, b] <- pp2p[,2, b] / dimsObj[2]                                      #
}                                                                              #
rm(iflux_file, astro_file)                                                     #
################################################################################

################################################################################
#-_----------------------------------_---_-----------------------------------_-#
#/ \--------------------------------/ \|/ \---------------------------------/ \#
#\_/--------------------------------\_/|\_/---------------------------------\_/#
#------------------------------------------------------------------------------#
################################################################################

################################################################################
############################ MW POL FILE PROCESSOR #############################
################################################################################
print(paste0("### Processing individual offsets detected objects to estimate ",#
             "MW polarization within target field ###"))                       #
                                                                               #
# Lists of files and filepaths to MW polarization files                        #
pathList <- list.files(mwpol_folder, full.names = T)                           #
fileList <- list.files(mwpol_folder, full.names = F)                           #
                                                                               #
# Lists of MW stars polarization files in different bands                      #
starsList <- list.files(MWstars_folder, full.names = T)                        #
                                                                               #
# Intersecting former list with object of interest                             #
obj_ind <- grep(og_obj, starsList)                                             #
starsList <- starsList[obj_ind]                                                #
                                                                               #
# Intersect resulting list with bands                                          #
star_bands <- band_str                                                         #
no_files <- NULL                                                               #
for(b in 1:nB){                                                                #
  if(length(grep(band_str[b], starsList)) == 0){                               #
    no_files <- append(no_files, b)                                            #
  }                                                                            #
}                                                                              #
rm(b)                                                                          #
                                                                               #
# Filter those that don't                                                      #
star_bands <- star_bands[-no_files]                                            #
rm(no_files)                                                                   #
                                                                               #
nB_stars <- length(star_bands)                                                 #
                                                                               #
totalN_stars  <- 0                                                             #
                                                                               #
# Counting the maximum amount of stars within one given band                   #
for(b in 1:nB_stars){                                                          #
  bs <- star_bands[b]                                                          #
  print(paste0("* Counting amount of sources saved for ", bs, " band... *"))   #
                                                                               #
  # Find files for band b                                                      #
  b_ind <- grep(star_bands[b], starsList)                                      #
                                                                               #
  # Determine how many files those are                                         #
  nO <- length(b_ind)                                                          #
                                                                               #
  # Create list of only those files                                            #
  offsetList <- starsList[b_ind]                                               #
                                                                               #
  temp_total <- 0                                                              #
  # Iterate over those files to check which band has the most sources          #
  for(o in 1:nO){                                                              #
    t_stars <- read.table(offsetList[o], sep = ";", header = T)                #
                                                                               #
    temp_total <- temp_total + dim(t_stars)[1]                                 #
  }                                                                            #
                                                                               #
  if(temp_total > totalN_stars){                                               #
    totalN_stars <- temp_total                                                 #
  }                                                                            #
}                                                                              #
rm(temp_total, t_stars, o_ind, b_ind, nO, offsetList, b, o)                    #
                                                                               #
# Creating arrays to hold all the information in the files selected above      #
gaia_pars <- c("Gaia.DR3.ID", "RA (2016)", "DEC (2016)", "u_RA (2016)",        #
               "u_DEC (2016)")                                                 #
gPar <- length(gaia_pars)                                                      #
                                                                               #
pol_pars <- c("Q", "U", "u_Q", "u_U" , "P", "X", "u_P", "u_X")                 #
pPar <- length(pol_pars)                                                       #
                                                                               #
par_list <- c(gaia_pars, pol_pars)                                             #
nPar <- length(par_list)                                                       #
                                                                               #
pos_par <- c("x", "y", "offset")                                               #
nPos <- length(pos_par)                                                        #
                                                                               #
stars_pol <- array(NA, dim = c(totalN_stars, nB, nPar),                        #
                   dimnames = list(NULL, band_str, par_list))                  #
stars_xy <- array(NA, dim = c(totalN_stars, nB, nPos),                         #
                  dimnames = list(NULL, band_str, pos_par))                    #
                                                                               #
############ Merging sources pol info across all targets offsets ############# #
print("### Merging sources pol info across all targets offsets ###")           #
for(b in 1:nB){                                                                #
  # Find files for band b                                                      #
  bs <- band_str[b]                                                            #
  b_ind <- grep(bs, starsList)                                                 #
                                                                               #
  # Determine how many files those are                                         #
  nO <- length(b_ind)                                                          #
                                                                               #
  if(nO == 0){                                                                 #
    next                                                                       #
  }                                                                            #
                                                                               #
  # Create list of only those files                                            #
  offsetList <- starsList[b_ind]                                               #
                                                                               #
  numbs_O <- NULL                                                              #
                                                                               #
  for(o in 1:nO){                                                              #
    t_numb <- str_split(str_split(offsetList[o], "#")[[1]][2], "-")[[1]][1]    #
    numbs_O <- c(numbs_O, as.integer(t_numb))                                  #
  }                                                                            #
                                                                               #
  print(paste0("* Merging sources pol info for band ", bs, ". *"))             #
                                                                               #
  switch(b,                                                                    #
         ref_pixs <- ref_pixs_U, ref_pixs <- ref_pixs_B,                       #
         ref_pixs <- ref_pixs_V, ref_pixs <- ref_pixs_R,                       #
         ref_pixs <- ref_pixs_I)                                               #
                                                                               #
  last_tNstars <- 0                                                            #
  # Iterate over those files to load the information                           #
  for(o in 1:nO){                                                              #
    o_ind <- grep(paste0("#", numbs_O[o]), offsetList)                         #
                                                                               #
    t_stars <- read.table(offsetList[o_ind], sep = ";", header = T)            #
                                                                               #
    tNstars <- dim(t_stars)[1]                                                 #
                                                                               #
    dx <- ref_pixs[1, "x"] - ref_pixs[o, "x"]                                  #
    dy <- ref_pixs[1, "y"] - ref_pixs[o, "y"]                                  #
    dra_x <- (ref_pixs[1, "ra"] - ref_pixs[o, "ra"]) * cosd[b] / scl_x         #
    ddec_y <- (ref_pixs[1, "dec"] - ref_pixs[o, "dec"]) / scl_y                #
                                                                               #
    difx <- round(dx + dra_x)                                                  #
    dify <- round(dy + ddec_y)                                                 #
    rm(dx, dy, dra_x, ddec_y)                                                  #
                                                                               #
    for(s in 1:tNstars){                                                       #
      stars_pol[s + last_tNstars, bs, 'Q'] <- t_stars[s, 'Q']                  #
      stars_pol[s + last_tNstars, bs, 'U'] <- t_stars[s, 'U']                  #
      stars_pol[s + last_tNstars, bs, 'P'] <- t_stars[s, 'P']                  #
      stars_pol[s + last_tNstars, bs, 'X'] <- t_stars[s, 'X']                  #
      stars_pol[s + last_tNstars, bs, 'u_Q'] <- t_stars[s, 'uQ']               #
      stars_pol[s + last_tNstars, bs, 'u_U'] <- t_stars[s, 'uU']               #
      stars_pol[s + last_tNstars, bs, 'u_P'] <- t_stars[s, 'uP']               #
      stars_pol[s + last_tNstars, bs, 'u_X'] <- t_stars[s, 'uX']               #
                                                                               #
      # Adjust pixel coordinates for stars in offsets != 1                     #
      stars_xy[s + last_tNstars, bs, "x"] <- t_stars[s, 'x'] - difx            #
      stars_xy[s + last_tNstars, bs, "y"] <- t_stars[s, 'y'] - dify            #
      stars_xy[s + last_tNstars, bs, "offset"] <- o                            #
    }                                                                          #
                                                                               #
    last_tNstars <- last_tNstars + tNstars                                     #
  }                                                                            #
}                                                                              #
rm(t_stars, o_ind, b_ind, nO, offsetList, last_tNstars, tNstars, starsList,    #
   b, bs, o, s)                                                                #
                                                                               #
# Opening file with information regarding stars psf                            #
stars_apt_r <- read.table(starsApt_R_path, header = T, sep = ";", dec = ".")   #
rm(starsApt_R_path)                                                            #
                                                                               #
# Locking rows, where band information is stored                               #
# Result of apply should have dimnames equal to the entries of "band_str"      #
apt_R <- ceiling(2 * apply(stars_apt_r, 1, max, na.rm = T))                    #
apt_R[which(is.infinite(apt_R), arr.ind = T)] <- NA                            #
                                                                               #
sync_R <- ceiling(apt_R / 2)                                                   #
rm(stars_apt_r)                                                                #
                                                                               #
######## Rejecting sources that don't have proper motions in Gaia DR3 ######## #
print("### Rejecting sources that don't have proper motions in Gaia DR3 ###")  #
for(b in 1:nB){                                                                #
  bs <- band_str[b]                                                            #
  R <- apt_R[bs]                                                               #
                                                                               #
  # Skipping bands without info                                                #
  if(length(which(!is.na(stars_xy[,bs,]), arr.ind = T)) == 0){                 #
    next                                                                       #
  }                                                                            #
                                                                               #
  print(paste0("* Rejecting sources without proper motion information in Gaia",#
               " DR3 for band ", bs, ". *"))                                   #
                                                                               #
  srcs_xy <- stars_xy[, bs, 1:2]                                               #
  n_srcs <- dim(srcs_xy)[1]                                                    #
                                                                               #
  # Checking Gaia for Stars in the FOV                                         #
  dx <- 2 * max(abs(c(1862 - iflux_pix[1, bs], iflux_pix[1, bs] - 192)))       #
  dy <- 2 * max(abs(c(1980 - iflux_pix[2, bs], iflux_pix[2, bs] - 355)))       #
  dra <- abs(dx * scl_x / cosd[bs])                                            #
  ddec <- dy * scl_y                                                           #
  rm(dx, dy)                                                                   #
                                                                               #
  gaia_stars <- get_GAIA_stars_coords(iflux_coord[1, bs], iflux_coord[2, bs],  #
                                      dra, ddec)                               #
  rm(dra, ddec)                                                                #
                                                                               #
  Ngaia <- dim(gaia_stars)[1]                                                  #
                                                                               #
  # Gaia data is in J2016, adjust to J2000 just in case                        #
  gaia_J2000 <- convert_yr_to_J2000(gaia_stars[,2:10])                         #
                                                                               #
  # Get the band b astrometry file                                             #
  astro_file <- astroPaths[bs]                                                 #
                                                                               #
  # Converting sources X, Y to RA, DEC using astrometry based on iflux file    #
  srcs_rd <- convert_pixel_to_wcs(srcs_xy[,c("x","y")], astro_file)            #
                                                                               #
  # Converting Gaia RA, DEC to X, Y using astrometry based on iflux file       #
  gaia_xy <- convert_wcs_to_pixel(gaia_J2000[,c("ra","dec")], astro_file)      #
                                                                               #
  # Rejecting Gaia Stars outside of OBS field                                  #
  out_ind <- NULL                                                              #
                                                                               #
  for(ng in 1:Ngaia){                                                          #
    if(gaia_xy[ng, 1] < 192 || gaia_xy[ng, 1] > 1862 ||                        #
       gaia_xy[ng, 2] < 355 || gaia_xy[ng, 2] > 1980){                         #
      out_ind <- c(out_ind, ng)                                                #
    }                                                                          #
  }                                                                            #
                                                                               #
  if(length(out_ind) != 0){                                                    #
    gaia_xy <- gaia_xy[-out_ind,]                                              #
    gaia_stars <- gaia_stars[-out_ind,]                                        #
    gaia_J2000 <- gaia_J2000[-out_ind,]                                        #
    Ngaia <- dim(gaia_stars)[1]                                                #
  }                                                                            #
  rm(out_ind)                                                                  #
                                                                               #
  ind_to_rm <- NULL                                                            #
  tol_dist <- sync_R[bs]                                                       #
                                                                               #
  # Checking for matching sources between Gaia stars and extracted sources     #
  for(ns in 1:n_srcs){                                                         #
    sx <- srcs_xy[ns, "x"]                                                     #
    sy <- srcs_xy[ns, "y"]                                                     #
                                                                               #
    # No point in comparing coordinates that are not there                     #
    # (required due to the way "stars_xy" is constructed)                      #
    if(is.na(sx) || is.na(sy)){                                                #
      ind_to_rm <- c(ind_to_rm, ns)                                            #
      next                                                                     #
    }                                                                          #
                                                                               #
    # Checking gaia pixel coordinates list                                     #
    found_flag <- F                                                            #
    ind_match <- NULL                                                          #
    match_dist <- NULL                                                         #
                                                                               #
    for(ng in 1:Ngaia){                                                        #
      gx <- gaia_xy[ng, "x"]                                                   #
      gy <- gaia_xy[ng, "y"]                                                   #
                                                                               #
      sg_dist <- dist(rbind(c(sx, sy), c(gx, gy)))                             #
                                                                               #
      if(sg_dist <= tol_dist){                                                 #
        ind_match <- c(ind_match, ng)                                          #
        match_dist <- c(match_dist, sg_dist)                                   #
        found_flag <- T                                                        #
      }                                                                        #
    }                                                                          #
                                                                               #
    if(!found_flag){                                                           #
      ind_to_rm <- c(ind_to_rm, ns)                                            #
    }else{                                                                     #
      if(length(ind_match) == 1){                                              #
        stars_pol[ns, bs, "Gaia.DR3.ID"] <- gaia_stars[ind_match,"source_id"]  #
        stars_pol[ns, bs, "RA (2016)"] <-  gaia_stars[ind_match, "ra"]         #
        stars_pol[ns, bs, "DEC (2016)"] <-  gaia_stars[ind_match, "dec"]       #
        stars_pol[ns, bs, "u_RA (2016)"] <-  gaia_stars[ind_match, "u_ra"]     #
        stars_pol[ns, bs, "u_DEC (2016)"] <-  gaia_stars[ind_match, "u_dec"]   #
      }else{                                                                   #
        best_match <- which(match_dist == min(match_dist), arr.ind = T)        #
        best_ind <- ind_match[best_match]                                      #
                                                                               #
        stars_pol[ns, bs, "Gaia.DR3.ID"] <- gaia_stars[best_ind,"source_id"]   #
        stars_pol[ns, bs, "RA (2016)"] <-  gaia_stars[best_ind, "ra"]          #
        stars_pol[ns, bs, "DEC (2016)"] <-  gaia_stars[best_ind, "dec"]        #
        stars_pol[ns, bs, "u_RA (2016)"] <-  gaia_stars[best_ind, "u_ra"]      #
        stars_pol[ns, bs, "u_DEC (2016)"] <-  gaia_stars[best_ind, "u_dec"]    #
        rm(best_match, best_ind)                                               #
      }                                                                        #
    }                                                                          #
  }                                                                            #
  rm(sx, sy, gx, gy)                                                           #
                                                                               #
  if(length(ind_to_rm) != 0){                                                  #
    stars_xy[ind_to_rm, bs,] <- NA                                             #
    stars_pol[ind_to_rm, bs,] <- NA                                            #
  }                                                                            #
                                                                               #
  switch(b,                                                                    #
         gaia_xy_U <- gaia_xy, gaia_xy_B <- gaia_xy, gaia_xy_V <- gaia_xy,     #
         gaia_xy_R <- gaia_xy, gaia_xy_I <- gaia_xy)                           #
}                                                                              #
rm(srcs_xy, ind_to_rm, ind_match, match_dist, found_flag, gaia_stars, bs, b,   #
   gaia_xy, Ngaia, ns, ng, n_srcs)                                             #
                                                                               #
# Creating a Catalog of MW Stars with Gaia Astrometry and FORS2 Polarimetry ## #
print(paste0("### Creating Catalog of MW Stars with Gaia astrometry and FORS2",#
             " Polarimetry  ###"))                                             #
                                                                               #
all_par_list <- gaia_pars                                                      #
                                                                               #
for(b in band_str){                                                            #
  all_par_list <- c(all_par_list, paste0(pol_pars, " (", b,")"))               #
}                                                                              #
                                                                               #
aPar <- length(all_par_list)                                                   #
MW_star_catlg <- array(NA, dim = c(0, aPar),                                   #
                       dimnames = list(NULL, all_par_list))                    #
                                                                               #
gaia_inds <- 1:gPar                                                            #
pol_inds <- gPar + 1:pPar                                                      #
                                                                               #
for(b in bands){                                                               #
  bs <- band_str[b]                                                            #
                                                                               #
  catlg_pol_inds <- pol_inds + (b - 1) * pPar                                  #
                                                                               #
  tmp_stars <- stars_pol[,bs,]                                                 #
                                                                               #
  # Clear NA entries                                                           #
  NA_rows <- which(is.na(apply(tmp_stars, 1, sum)), arr.ind = T)               #
  tmp_stars <- tmp_stars[-NA_rows,]                                            #
  rm(NA_rows)                                                                  #
                                                                               #
  n_ts <- dim(tmp_stars)[1]                                                    #
  if(n_ts == 0){                                                               #
    next                                                                       #
  }                                                                            #
                                                                               #
  # Determining outliers (likely intrinsically polarized sources)              #
  medP <- median(tmp_stars[,"P"])                                              #
  uncP <- unc_median(tmp_stars[,"P"], tmp_stars[,"u_P"], 0)                    #
                                                                               #
  outP_inds <- which(tmp_stars[,"P"] > (medP + 3 * uncP), arr.ind = T)         #
                                                                               #
  # Estimating MW Stokes parameters                                            #
  mwQ <- median(tmp_stars[-outP_inds,"Q"])                                     #
  mwQ_unc <- unc_median(tmp_stars[-outP_inds,"Q"],                             #
                        tmp_stars[-outP_inds,"u_Q"], 0)                        #
  mwU <- median(tmp_stars[-outP_inds,"U"])                                     #
  mwU_unc <- unc_median(tmp_stars[-outP_inds,"U"],                             #
                        tmp_stars[-outP_inds,"u_U"], 0)                        #
                                                                               #
  # Subtracting MW Stokes estimate from all sources                            #
  tmp_stars[,"Q"] <- tmp_stars[,"Q"] - mwQ                                     #
  tmp_stars[,"U"] <- tmp_stars[,"U"] - mwU                                     #
  for(s in 1:n_ts){                                                            #
    tmp_stars[s,"u_Q"] <- unc_add(tmp_stars[s,"u_Q"], mwQ_unc)                 #
    tmp_stars[s,"u_U"] <- unc_add(tmp_stars[s,"u_U"], mwU_unc)                 #
  }                                                                            #
  rm(mwQ, mwU, mwQ_unc, mwU_unc)                                               #
                                                                               #
  # And recalculating P and X                                                  #
  tmp_Ps <- debiased_P_from_QU(tmp_stars[,"Q"], tmp_stars[,"U"],               #
                               tmp_stars[,"u_Q"], tmp_stars[,"u_U"])           #
  tmp_stars[,"P"] <- tmp_Ps$Data                                               #
  tmp_stars[,"u_P"] <- tmp_Ps$Unc                                              #
  rm(tmp_Ps)                                                                   #
                                                                               #
  tmp_Xs <- X_from_QU(tmp_stars[,"Q"], tmp_stars[,"U"],                        #
                      tmp_stars[,"u_Q"], tmp_stars[,"u_U"])                    #
  tmp_stars[,"X"] <- tmp_Xs$Data                                               #
  tmp_stars[,"u_X"] <- tmp_Xs$Unc                                              #
  rm(tmp_Xs)                                                                   #
                                                                               #
  for(s in 1:n_ts){                                                            #
    # Checking if there is an entry for this star                              #
    this_id <- tmp_stars[s, "Gaia.DR3.ID"]                                     #
    this_ind <- which(MW_star_catlg[,"Gaia.DR3.ID"] == this_id, arr.ind = T)   #
                                                                               #
    if(length(this_ind) == 0){                                                 #
      # Create a new entry                                                     #
      MW_star_catlg <- rbind(MW_star_catlg, rep(NA, aPar))                     #
                                                                               #
      # Index to that entry's row                                              #
      catlg_row <- dim(MW_star_catlg)[1]                                       #
                                                                               #
      # Fill new entry                                                         #
      MW_star_catlg[catlg_row, gaia_inds] <- tmp_stars[s, gaia_inds]           #
      MW_star_catlg[catlg_row, catlg_pol_inds] <- tmp_stars[s, pol_inds]       #
      rm(catlg_row)                                                            #
    }else{                                                                     #
      # Complete existing entry                                                #
      MW_star_catlg[this_ind, catlg_pol_inds] <- tmp_stars[s, pol_inds]        #
    }                                                                          #
  }                                                                            #
  rm(s, this_id, this_ind)                                                     #
}                                                                              #
rm(b, bs, tmp_stars, n_ts)                                                     #
                                                                               #
# Saving Catalog file                                                          #
write.table(MW_star_catlg, paste0(MWstars_folder, "MW_Stars-Field_", og_obj,   #
                                  ".csv"), sep = ";", dec = ".")               #
rm(MW_star_catlg)                                                              #
                                                                               #
############## Creating pixel coordinate mask template for stars ############# #
print(paste0("### Creating pixel coordinate mask template for MW stars  ###")) #
apt_R <- 1.5 * apt_R                                                           #
                                                                               #
star_mask_paths <- array("", dim = nB, dimnames = list(band_str))              #
                                                                               #
# Creating template for each band                                              #
for(b in 1:nB){                                                                #
  bs <- band_str[b]                                                            #
  R <- apt_R[bs]                                                               #
                                                                               #
  # Skipping bands without info                                                #
  if(is.na(R)){                                                                #
    next                                                                       #
  }                                                                            #
                                                                               #
  # Get the band b astrometry file                                             #
  astro_file <- astroPaths[bs]                                                 #
                                                                               #
  # Extract the common name                                                    #
  t_split <- strsplit(astro_file, "/")                                         #
  basename <- strsplit(t_split[[1]][length(t_split[[1]])], "Astro")[[1]][1]    #
                                                                               #
  star_mask_paths[b] <- paste0(star_mask_folder, basename,"_MW-Star-Mask.fits")#
  maskflag <- file.exists(star_mask_paths[b])                                  #
                                                                               #
  if(maskflag){                                                                #
    print(paste0("* Star mask for ", bs, " band already exists... *"))         #
  }else{                                                                       #
    # Defining template coordinates for pixels inside a circle                 #
    deltas <- -R:R                                                             #
    std_pts <- NULL                                                            #
    for(x in deltas){                                                          #
      for(y in deltas){                                                        #
        if(dist(rbind(c(0, 0), c(x, y))) <= R){                                #
          std_pts <- rbind(std_pts, c(x, y))                                   #
        }                                                                      #
      }                                                                        #
    }                                                                          #
                                                                               #
    mask_stars <- array(1, dim = dimsObj)                                      #
                                                                               #
    switch(b,                                                                  #
           gaia_xy <- gaia_xy_U, gaia_xy <- gaia_xy_B, gaia_xy <- gaia_xy_V,   #
           gaia_xy <- gaia_xy_R, gaia_xy <- gaia_xy_I)                         #
                                                                               #
    N_pts <- dim(std_pts)[1]                                                   #
    gaiaN <- dim(gaia_xy)[1]                                                   #
                                                                               #
    print(paste0("* Creating star mask for ", bs, " band... *"))               #
                                                                               #
    for(nas in 1:gaiaN){                                                       #
      xc <- gaia_xy[nas, 'x']                                                  #
      yc <- gaia_xy[nas, 'y']                                                  #
                                                                               #
      pts <- std_pts                                                           #
      pts[,1] <- pts[,1] + xc                                                  #
      pts[,2] <- pts[,2] + yc                                                  #
                                                                               #
      for(p in 1:N_pts){                                                       #
        x <- pts[p, 1]                                                         #
        y <- pts[p, 2]                                                         #
        mask_stars[x, y] <- NA                                                 #
        mask_stars[x, y] <- NA                                                 #
      }                                                                        #
    }                                                                          #
    # Saving MW star mask fits                                                 #
    writeFITSim(mask_stars, star_mask_paths[b])                                #
    rm(mask_stars, deltas, x, y, N_pts, gaiaN, xc, yc)                         #
  }                                                                            #
}                                                                              #
rm(R, b, bs)                                                                   #
                                                                               #
############# Crossing information with Sources from each Offset ############# #
print("### Merging MW sources info across all targets offsets ###")            #
                                                                               #
mw_labels <- c("raw", "instC")                                                 #
stokes_par <- c("Q","U")                                                       #
pol_par <- c("P", "X")                                                         #
mw_dim <- c(length(stokes_par), length(munc_par), nB, length(mw_labels))       #
mw_dname <- list(stokes_par, munc_par, band_str, mw_labels)                    #
mw_out_labels <- c(mw_labels, "allC")                                          #
mw_typesN <- length(mw_out_labels)                                             #
                                                                               #
mwQU <- array(0, dim = mw_dim, dimnames = mw_dname)                            #
                                                                               #
for(b in bands){                                                               #
  bs <- band_str[b]                                                            #
  b_ind <- grep(bs, fileList)                                                  #
                                                                               #
  if(length(b_ind) == 0){                                                      #
    next                                                                       #
  }                                                                            #
                                                                               #
  print(paste0("* Merging list of MW sources Stokes parameters in ", bs,       #
               " band *"))                                                     #
                                                                               #
  tpathList <- pathList[b_ind]                                                 #
  tfileList <- fileList[b_ind]                                                 #
                                                                               #
  # Needed to compare raw values between offsets                               #
  rawind <- grep("raw", tfileList)                                             #
  # Needed because chromatic correction of galaxy polarization must be         #
  # applied after all corrections other                                        #
  instind <- grep("instC", tfileList)                                          #
  # Needed to build the corrected values sample                                #
  allind <- grep("allC", tfileList)                                            #
                                                                               #
  deb <- T                                                                     #
                                                                               #
  for(type in 1:mw_typesN){                                                    #
    switch(type,                                                               #
           c(ttpathList <- tpathList[rawind], ttfileList <- tfileList[rawind], #
             basename <- strsplit(ttfileList, "_#")[[1]][1],                   #
             typename <- "raw"),                                               #
           c(ttpathList <- tpathList[instind], typename <- "instC",            #
             ttfileList <- tfileList[instind],                                 #
             basename <- strsplit(ttfileList, "_#")[[1]][1]),                  #
           c(ttpathList <- NA, ttfileList <- tfileList[allind],                #
             basename <- strsplit(ttfileList, "_#")[[1]][1],                   #
             typename <- "allC", deb <- T))                                    #
                                                                               #
    N <- length(ttfileList)                                                    #
                                                                               #
    if(type == mw_typesN){                                                     #
      mwpol <- array(NA, dim = c(dim(stars_pol)[1], 5),                        #
                     dimnames = list(NULL, c("Q.Data", "U.Data", "Q.Unc",      #
                                             "U.Unc", "offset")))              #
                                                                               #
      mwpol[,1:4] <- stars_pol[, bs, c("Q", "U", "u_Q", "u_U")]                #
      mwpol[,5] <- stars_xy[, bs, "offset"]                                    #
                                                                               #
      NA_rows <- which(is.na(apply(mwpol, 1, sum)), arr.ind = T)               #
                                                                               #
      mwpol <- mwpol[-NA_rows,]                                                #
      rm(NA_rows)                                                              #
    }else{                                                                     #
      offset <- NULL                                                           #
                                                                               #
      for(n in 1:N){                                                           #
        if(n == 1){                                                            #
          mwpol <- read.table(ttpathList[n], sep = ";", header = T)            #
          offset <- c(offset, rep(n, dim(mwpol)[1]))                           #
        }else{                                                                 #
          temp <- read.table(ttpathList[n], sep = ";", header = T)             #
          mwpol <- rbind(mwpol, temp)                                          #
          offset <- c(offset, rep(n, dim(temp)[1]))                            #
        }                                                                      #
      }                                                                        #
                                                                               #
      mwpol <- cbind(mwpol, offset)                                            #
      rm(n, offset)                                                            #
    }                                                                          #
                                                                               #
    oqu_dim <- c(length(stokes_par), length(munc_par), N)                      #
    oqu_dname <- list(stokes_par, munc_par, NULL)                              #
    opx_dim <- c(length(pol_par), length(munc_par), N)                         #
    opx_dname <- list(pol_par, munc_par, NULL)                                 #
    qu_dim <- c(length(stokes_par), length(munc_par))                          #
    qu_dname <- list(stokes_par, munc_par)                                     #
    px_dim <- c(length(pol_par), length(munc_par))                             #
    px_dname <- list(pol_par, munc_par)                                        #
                                                                               #
    mw_o_QU <- array(NA, dim = oqu_dim, dimnames = oqu_dname)                  #
    mw_o_PX <- array(NA, dim = opx_dim, dimnames = opx_dname)                  #
    mw_QU <- array(NA, dim = qu_dim, dimnames = qu_dname)                      #
    mw_QUf <- array(NA, dim = qu_dim, dimnames = qu_dname)                     #
    mw_PX <- array(NA, dim = px_dim, dimnames = px_dname)                      #
    mw_Pf <- array(NA, dim = c(1, px_dim[2]),                                  #
                   dimnames = list(c('P'), munc_par))                          #
                                                                               #
    for(n in 1:N){                                                             #
      ns <- which(mwpol[, "offset"] == n, arr.ind = T)                         #
                                                                               #
      mw_o_QU["Q", "Data", n] <- median(mwpol[ns, "Q.Data"])                   #
      mw_o_QU["U", "Data", n] <- median(mwpol[ns, "U.Data"])                   #
      mw_o_QU["Q", "Unc", n] <- unc_median(mwpol[ns, "Q.Data"],                #
                                           mwpol[ns, "Q.Unc"], 0)              #
      mw_o_QU["U", "Unc", n] <- unc_median(mwpol[ns, "U.Data"],                #
                                           mwpol[ns, "U.Unc"], 0)              #
    }                                                                          #
                                                                               #
    mw_QU["Q", "Data"] <- median(mwpol[,"Q.Data"])                             #
    mw_QU["U", "Data"] <- median(mwpol[,"U.Data"])                             #
    mw_QU["Q", "Unc"] <- unc_median(mwpol[,"Q.Data"], mwpol[,"Q.Unc"], 0)      #
    mw_QU["U", "Unc"] <- unc_median(mwpol[,"U.Data"], mwpol[,"U.Unc"], 0)      #
                                                                               #
    # These will be used for reduction of merged maps                          #
    if(type == 1 || type == 2){                                                #
      mwQU["Q", "Data", b, type] <- mw_QU["Q", "Data"]                         #
      mwQU["Q", "Unc", b, type] <- mw_QU["Q", "Unc"]                           #
      mwQU["U", "Data", b, type] <- mw_QU["U", "Data"]                         #
      mwQU["U", "Unc", b, type] <- mw_QU["U", "Unc"]                           #
    }                                                                          #
                                                                               #
    if(deb){                                                                   #
      temp1 <- debiased_P_from_QU(mw_QU['Q', "Data"], mw_QU['U', "Data"],      #
                                  mw_QU['Q', "Unc"], mw_QU['U', "Unc"])        #
      temp2 <- debiased_P_from_QU(mw_o_QU['Q', "Data",], mw_o_QU['U', "Data",],#
                                  mw_o_QU['Q', "Unc",], mw_o_QU['U', "Unc",])  #
      temp3 <- debiased_P_from_QU(mwpol[,"Q.Data"], mwpol[,"U.Data"],          #
                                  mwpol[,"Q.Unc"], mwpol[,"U.Unc"])            #
    }else{                                                                     #
      temp1 <- P_from_QU(mw_QU['Q', "Data"], mw_QU['U', "Data"],               #
                         mw_QU['Q', "Unc"], mw_QU['U', "Unc"])                 #
      temp2 <- P_from_QU(mw_o_QU['Q', "Data",], mw_o_QU['U', "Data",],         #
                         mw_o_QU['Q', "Unc",], mw_o_QU['U', "Unc",])           #
      temp3 <- P_from_QU(mwpol[,"Q.Data"], mwpol[,"U.Data"],                   #
                         mwpol[,"Q.Unc"], mwpol[,"U.Unc"])                     #
    }                                                                          #
    mw_PX['P', "Data"] <- temp1$Data                                           #
    mw_PX['P', "Unc"] <- temp1$Unc                                             #
    mw_o_PX['P', "Data",] <- temp2$Data                                        #
    mw_o_PX['P', "Unc",] <- temp2$Unc                                          #
    mwpolP <- temp3$Data                                                       #
    mwpolP_unc <- temp3$Unc                                                    #
                                                                               #
    # Determining outliers (likely intrinsically polarized sources)            #
    medP <- median(mwpolP)                                                     #
    uncP <- unc_median(mwpolP, mwpolP_unc, 0)                                  #
                                                                               #
    outP_inds <- which(mwpolP > (medP + 3 * uncP), arr.ind = T)                #
                                                                               #
    if(length(outP_inds) != 0){                                                #
      mwpolPf <- mwpolP[-outP_inds]                                            #
      mwpolf <- mwpol[-outP_inds,]                                             #
                                                                               #
      mw_QUf["Q", "Data"] <- median(mwpolf[,"Q.Data"])                         #
      mw_QUf["U", "Data"] <- median(mwpolf[,"U.Data"])                         #
      mw_QUf["Q", "Unc"] <- unc_median(mwpolf[,"Q.Data"], mwpolf[,"Q.Unc"], 0) #
      mw_QUf["U", "Unc"] <- unc_median(mwpolf[,"U.Data"], mwpolf[,"U.Unc"], 0) #
                                                                               #
      if(deb){                                                                 #
        temp1 <- debiased_P_from_QU(mw_QUf['Q', "Data"], mw_QUf['U', "Data"],  #
                                    mw_QUf['Q', "Unc"], mw_QUf['U', "Unc"])    #
      }else{                                                                   #
        temp1 <- P_from_QU(mw_QUf['Q', "Data"], mw_QUf['U', "Data"],           #
                           mw_QUf['Q', "Unc"], mw_QUf['U', "Unc"])             #
      }                                                                        #
                                                                               #
      mw_Pf['P', "Data"] <- temp1$Data                                         #
      mw_Pf['P', "Unc"] <- temp1$Unc                                           #
    }                                                                          #
                                                                               #
    temp1 <- X_from_QU(mw_QU['Q', "Data"], mw_QU['U', "Data"],                 #
                       mw_QU['Q', "Unc"], mw_QU['U', "Unc"])                   #
    temp2 <- X_from_QU(mw_o_QU['Q', "Data",], mw_o_QU['U', "Data",],           #
                       mw_o_QU['Q', "Unc",], mw_o_QU['U', "Unc",])             #
    mw_PX['X', "Data"] <- temp1$Data                                           #
    mw_PX['X', "Unc"] <- temp1$Unc                                             #
    mw_o_PX['X', "Data",] <- temp2$Data                                        #
    mw_o_PX['X', "Unc",] <- temp2$Unc                                          #
    rm(temp1, temp2, temp3)                                                    #
                                                                               #
    write.table(mwpol, paste0(mwPol_path, basename, "-Merged-Sources_",        #
                              typename, ".csv"), sep = ";", dec = ".")         #
    write.table(mw_QU, paste0(mwPol_path, basename, "-Merged-QU_", typename,   #
                              ".csv"), sep = ";", dec = ".")                   #
    write.table(mw_PX, paste0(mwPol_path, basename, "-Merged-PX_", typename,   #
                              ".csv"), sep = ";", dec = ".")                   #
                                                                               #
    pdf_mw_det <- paste0(mwPol_path, basename, "-Merged-QvsU_", typename,      #
                         "_detailed.pdf")                                      #
                                                                               #
    if(length(outP_inds) != 0){                                                #
      pdf(pdf_mw_det, width = 31, height = 15)                                 #
      par(mfrow = c(1, 2))                                                     #
    }else{                                                                     #
      pdf(pdf_mw_det, width = 15, height = 15)                                 #
    }                                                                          #
    rm(pdf_mw_det)                                                             #
                                                                               #
    med_cols <- c("cornflowerblue", "darkturquoise", "darkgreen", "darkcyan",  #
                  "limegreen")                                                 #
    pt_cols <- gray.colors(N, 0, 0.6, 1.5)                                     #
    os_pch <- 20                                                               #
    mo_pch <- 15                                                               #
    m_pch <- 14                                                                #
    m_col <- "red"                                                             #
                                                                               #
    print(paste0("** Plotting seperate offset unfiltered sources for type ",   #
                 typename, " **"))                                             #
    ## Plotting unfiltered sources                                             #
    limmin <- min(c(mwpol[,'Q.Data'] - mwpol[,'Q.Unc'],                        #
                    -(mwpol[,'Q.Data'] + mwpol[,'Q.Unc']),                     #
                    mwpol[,'U.Data'] - mwpol[,'U.Unc'],                        #
                    -(mwpol[,'U.Data'] + mwpol[,'U.Unc']), -max(mwpolP)))      #
    limmax <- max(c(mwpol[,'Q.Data'] + mwpol[,'Q.Unc'],                        #
                    -(mwpol[,'Q.Data'] - mwpol[,'Q.Unc']),                     #
                    mwpol[,'U.Data'] + mwpol[,'U.Unc'],                        #
                    -(mwpol[,'U.Data'] - mwpol[,'U.Unc']), max(mwpolP)))       #
                                                                               #
    limabs <- max(abs(c(limmin, limmax)))                                      #
    axlim <- c(-limabs, limabs) * 100                                          #
    rm(limmin, limmax, limabs)                                                 #
                                                                               #
    oi <- 0                                                                    #
    oi_ind <- NULL                                                             #
    while(length(oi_ind) == 0 & oi <= N){                                      #
      oi <- oi + 1                                                             #
      oi_ind <- which(mwpol[,"offset"] == oi, arr.ind = T)                     #
    }                                                                          #
                                                                               #
    if(oi <= N){                                                               #
      # Stars in offset 1                                                      #
      plot(mwpol[oi_ind, 'Q.Data'] * 100, mwpol[oi_ind, 'U.Data'] * 100,       #
           xlab = 'q (%)', ylab = 'u (%)', xlim = axlim, ylim = axlim,         #
           pch = os_pch, cex = 1, col = pt_cols[oi], cex.axis = 2,             #
           main = "q vs u of MW stars at all offsets")                         #
      arrows(x0 = mwpol[oi_ind, 'Q.Data'] * 100,                               #
             y0 = (mwpol[oi_ind, 'U.Data'] - mwpol[oi_ind, 'U.Unc']) * 100,    #
             y1 = (mwpol[oi_ind, 'U.Data'] + mwpol[oi_ind, 'U.Unc']) * 100,    #
             code = 3, col = pt_cols[oi], length = 0.02, angle = 90)           #
      arrows(x0 = (mwpol[oi_ind, 'Q.Data'] - mwpol[oi_ind, 'Q.Unc']) * 100,    #
             y0 = mwpol[oi_ind, 'U.Data'] * 100,                               #
             x1 = (mwpol[oi_ind, 'Q.Data'] + mwpol[oi_ind, 'Q.Unc']) * 100,    #
             code = 3, col = pt_cols[oi], length = 0.02, angle = 90)           #
                                                                               #
      # Plotting Pol Magnitude Circle                                          #
      for(ns in oi_ind){                                                       #
        curve(sqrt((mwpolP[ns] * 100)^2 - x^2),                                #
              from = -mwpolP[ns] * 100, to = mwpolP[ns] * 100,                 #
              add = TRUE, col = pt_cols[oi], lty = 3, lwd=1)                   #
        curve(-sqrt((mwpolP[ns] * 100)^2 - x^2),                               #
              from = -mwpolP[ns] * 100, to = mwpolP[ns] * 100,                 #
              add = TRUE, col = pt_cols[oi], lty = 3, lwd=1)                   #
      }                                                                        #
                                                                               #
      # Median of stars in offset 1                                            #
      points(mw_o_QU['Q', "Data", oi] * 100, mw_o_QU['U', "Data", oi] * 100,   #
             pch = mo_pch, cex = 1.5, col = med_cols[oi])                      #
      arrows(x0 = mw_o_QU['Q', "Data", oi] * 100,                              #
             y0 = (mw_o_QU['U', "Data", oi] - mw_o_QU['U', "Unc", oi]) * 100,  #
             y1 = (mw_o_QU['U', "Data", oi] + mw_o_QU['U', "Unc", oi]) * 100,  #
             code = 3, length = 0.03, angle = 90, col = med_cols[oi])          #
      arrows(x0 = (mw_o_QU['Q', "Data", oi] - mw_o_QU['Q', "Unc", oi]) * 100,  #
             y0 = mw_o_QU['U', "Data", oi] * 100,                              #
             x1 = (mw_o_QU['Q', "Data", oi] + mw_o_QU['Q', "Unc", oi]) * 100,  #
             code = 3, length = 0.03, angle = 90, col = med_cols[oi])          #
                                                                               #
      # Plotting Pol Magnitude Circle                                          #
      curve(sqrt((mw_o_PX['P', "Data", oi] * 100)^2 - x^2),                    #
            from = - mw_o_PX['P', "Data", oi] * 100,                           #
            to = mw_o_PX['P', "Data", oi] * 100, add = TRUE,                   #
            col = med_cols[oi], lty = 1, lwd = 1.5)                            #
      curve(-sqrt((mw_o_PX['P', "Data", oi] * 100)^2 - x^2),                   #
            from = - mw_o_PX['P', "Data", oi] * 100,                           #
            to = mw_o_PX['P', "Data", oi] * 100, add = TRUE,                   #
            col = med_cols[oi], lty = 1, lwd = 1.5)                            #
                                                                               #
      for(n in (oi + 1):N){                                                    #
        o_ind <- which(mwpol[,"offset"] == n, arr.ind = T)                     #
                                                                               #
        if(length(o_ind) != 0){                                                #
          # Stars in offset n                                                  #
          points(mwpol[o_ind, 'Q.Data'] * 100, mwpol[o_ind, 'U.Data'] * 100,   #
                 pch = os_pch, cex = 1, col = pt_cols[n])                      #
          arrows(x0 = mwpol[o_ind, 'Q.Data'] * 100,                            #
                 y0 = (mwpol[o_ind, 'U.Data'] - mwpol[o_ind, 'U.Unc']) * 100,  #
                 y1 = (mwpol[o_ind, 'U.Data'] + mwpol[o_ind, 'U.Unc']) * 100,  #
                 code = 3, col = pt_cols[n], length = 0.02, angle = 90)        #
          arrows(x0 = (mwpol[o_ind, 'Q.Data'] - mwpol[o_ind, 'Q.Unc']) * 100,  #
                 y0 = mwpol[o_ind, 'U.Data'] * 100,                            #
                 x1 = (mwpol[o_ind, 'Q.Data'] + mwpol[o_ind, 'Q.Unc']) * 100,  #
                 code = 3, col = pt_cols[n], length = 0.02, angle = 90)        #
                                                                               #
          # Plotting Pol Magnitude Circle                                      #
          for(ns in o_ind){                                                    #
            curve(sqrt((mwpolP[ns] * 100)^2 - x^2),                            #
                  from = -mwpolP[ns] * 100, to = mwpolP[ns] * 100,             #
                  add = TRUE, col = pt_cols[n], lty = 3, lwd = 1)              #
            curve(-sqrt((mwpolP[ns] * 100)^2 - x^2),                           #
                  from = -mwpolP[ns] * 100, to = mwpolP[ns] * 100,             #
                  add = TRUE, col = pt_cols[n], lty = 3, lwd = 1)              #
          }                                                                    #
                                                                               #
          # Median of stars in offset n                                        #
          points(mw_o_QU['Q', "Data", n] * 100, mw_o_QU['U', "Data", n] * 100, #
                 pch = mo_pch, cex = 1.5, col = med_cols[n])                   #
          arrows(x0 = mw_o_QU['Q', "Data", n] * 100,                           #
                 y0 = (mw_o_QU['U', "Data", n] - mw_o_QU['U', "Unc", n]) * 100,#
                 y1 = (mw_o_QU['U', "Data", n] + mw_o_QU['U', "Unc", n]) * 100,#
                 code = 3, length = 0.03, angle = 90, col = med_cols[n])       #
          arrows(x0 = (mw_o_QU['Q', "Data", n] - mw_o_QU['Q', "Unc", n]) * 100,#
                 y0 = mw_o_QU['U', "Data", n] * 100,                           #
                 x1 = (mw_o_QU['Q', "Data", n] + mw_o_QU['Q', "Unc", n]) * 100,#
                 code = 3, length = 0.03, angle = 90, col = med_cols[n])       #
                                                                               #
          # Plotting Pol Magnitude Circle                                      #
          curve(sqrt((mw_o_PX['P', "Data", n] * 100)^2 - x^2),                 #
                from = - mw_o_PX['P', "Data", n] * 100,                        #
                to = mw_o_PX['P', "Data", n] * 100, add = TRUE,                #
                col = med_cols[n], lty = 1, lwd = 1.5)                         #
          curve(-sqrt((mw_o_PX['P', "Data", n] * 100)^2 - x^2),                #
                from = - mw_o_PX['P', "Data", n] * 100,                        #
                to = mw_o_PX['P', "Data", n] * 100, add = TRUE,                #
                col = med_cols[n], lty = 1, lwd = 1.5)                         #
        }                                                                      #
      }                                                                        #
                                                                               #
      # Median of all stars                                                    #
      points(mw_QU['Q', "Data"] * 100, mw_QU['U', "Data"] * 100, pch = m_pch,  #
             cex = 2, col = "red")                                             #
      arrows(x0 = mw_QU['Q', "Data"] * 100,                                    #
             y0 = (mw_QU['U', "Data"] - mw_QU['U', "Unc"]) * 100,              #
             y1 = (mw_QU['U', "Data"] + mw_QU['U', "Unc"]) * 100, code = 3,    #
             length = 0.03, angle = 90, col = "red")                           #
      arrows(x0 = (mw_QU['Q', "Data"] - mw_QU['Q', "Unc"]) * 100,              #
             y0 = mw_QU['U', "Data"] * 100,                                    #
             x1 = (mw_QU['Q', "Data"] + mw_QU['Q', "Unc"]) * 100, code = 3,    #
             length = 0.03, angle = 90, col = "red")                           #
                                                                               #
      # Plotting Pol Magnitude Circle                                          #
      curve(sqrt((mw_PX['P', "Data"] * 100)^2 - x^2),                          #
            from = - mw_PX['P', "Data"] * 100, to = mw_PX['P', "Data"] * 100,  #
            add = TRUE, col = "red", lty = 1, lwd = 2)                         #
      curve(-sqrt((mw_PX['P', "Data"] * 100)^2 - x^2),                         #
            from = - mw_PX['P', "Data"] * 100, to = mw_PX['P', "Data"] * 100,  #
            add = TRUE, col = "red", lty = 1, lwd = 2)                         #
                                                                               #
      # Setting up legend                                                      #
      leg_labs <- NULL                                                         #
      leg_pchs <- NULL                                                         #
      leg_cols <- NULL                                                         #
                                                                               #
      for(n in 1:N){                                                           #
        leg_labs <- c(leg_labs, paste0("Stars in Offset #", n),                #
                      paste0("Median of Offset #", n))                         #
        leg_pchs <- c(leg_pchs, os_pch, mo_pch)                                #
        leg_cols <- c(leg_cols, pt_cols[n], med_cols[n])                       #
      }                                                                        #
                                                                               #
      leg_labs <- c(leg_labs, "Median of Field")                               #
      leg_pchs <- c(leg_pchs, m_pch)                                           #
      leg_cols <- c(leg_cols, m_col)                                           #
                                                                               #
      # Legend to points                                                       #
      legend(axlim[2] * 0.5, axlim[2], legend = leg_labs, pch = leg_pchs,      #
             col = leg_cols, cex = 2)                                          #
    }                                                                          #
                                                                               #
    print(paste0("** Plotting seperate offset filtered sources for type ",     #
                 typename, " **"))                                             #
    if(length(outP_inds) != 0){                                                #
      ## Plotting filtered sources                                             #
      limmin <- min(c(mwpolf[,'Q.Data'] - mwpolf[,'Q.Unc'],                    #
                      -(mwpolf[,'Q.Data'] + mwpolf[,'Q.Unc']),                 #
                      mwpolf[,'U.Data'] - mwpolf[,'U.Unc'],                    #
                      -(mwpolf[,'U.Data'] + mwpolf[,'U.Unc']), -max(mwpolPf))) #
      limmax <- max(c(mwpolf[,'Q.Data'] + mwpolf[,'Q.Unc'],                    #
                      -(mwpolf[,'Q.Data'] - mwpolf[,'Q.Unc']),                 #
                      mwpolf[,'U.Data'] + mwpolf[,'U.Unc'],                    #
                      -(mwpolf[,'U.Data'] - mwpolf[,'U.Unc']), max(mwpolPf)))  #
                                                                               #
      limabs <- max(abs(c(limmin, limmax)))                                    #
      axlim <- c(-limabs, limabs) * 100                                        #
      rm(limmin, limmax, limabs)                                               #
                                                                               #
      oi <- 0                                                                  #
      oi_ind <- NULL                                                           #
      while(length(oi_ind) == 0 & oi <= N){                                    #
        oi <- oi + 1                                                           #
        oi_ind <- which(mwpolf[,"offset"] == oi, arr.ind = T)                  #
      }                                                                        #
                                                                               #
      if(oi <= N){                                                             #
        # Stars in offset 1                                                    #
        plot(mwpolf[oi_ind, 'Q.Data'] * 100, mwpolf[oi_ind, 'U.Data'] * 100,   #
             xlab = 'q (%)', ylab = 'u (%)', xlim = axlim, ylim = axlim,       #
             pch = 20, cex = 1, col = pt_cols[oi], cex.axis = 2,               #
             main = "q vs u of MW stars at all offsets (non-Intrs-Pol)")       #
        arrows(x0 = mwpolf[oi_ind, 'Q.Data'] * 100,                            #
               y0 = (mwpolf[oi_ind, 'U.Data'] - mwpolf[oi_ind, 'U.Unc']) * 100,#
               y1 = (mwpolf[oi_ind, 'U.Data'] + mwpolf[oi_ind, 'U.Unc']) * 100,#
               code = 3, col = pt_cols[oi], length = 0.02, angle = 90)         #
        arrows(x0 = (mwpolf[oi_ind, 'Q.Data'] - mwpolf[oi_ind, 'Q.Unc']) * 100,#
               y0 = mwpolf[oi_ind, 'U.Data'] * 100,                            #
               x1 = (mwpolf[oi_ind, 'Q.Data'] + mwpolf[oi_ind, 'Q.Unc']) * 100,#
               code = 3, col = pt_cols[oi], length = 0.02, angle = 90)         #
                                                                               #
        # Plotting Pol Magnitude Circle                                        #
        for(ns in oi_ind){                                                     #
          curve(sqrt((mwpolPf[ns] * 100)^2 - x^2),                             #
                from = -mwpolPf[ns] * 100, to = mwpolPf[ns] * 100,             #
                add = TRUE, col = pt_cols[oi], lty = 3, lwd = 1)               #
          curve(-sqrt((mwpolPf[ns] * 100)^2 - x^2),                            #
                from = -mwpolPf[ns] * 100, to = mwpolPf[ns] * 100,             #
                add = TRUE, col = pt_cols[oi], lty = 3, lwd = 1)               #
        }                                                                      #
                                                                               #
        # Median of stars in offset 1                                          #
        points(mw_o_QU['Q', "Data", oi] * 100, mw_o_QU['U', "Data", oi] * 100, #
               pch = 15, cex = 1.5, col = med_cols[oi])                        #
        arrows(x0 = mw_o_QU['Q', "Data", oi] * 100,                            #
               y0 = (mw_o_QU['U', "Data", oi] - mw_o_QU['U', "Unc", oi]) * 100,#
               y1 = (mw_o_QU['U', "Data", oi] + mw_o_QU['U', "Unc", oi]) * 100,#
               code = 3, length = 0.03, angle = 90, col = med_cols[oi])        #
        arrows(x0 = (mw_o_QU['Q', "Data", oi] - mw_o_QU['Q', "Unc", oi]) * 100,#
               y0 = mw_o_QU['U', "Data", oi] * 100,                            #
               x1 = (mw_o_QU['Q', "Data", oi] + mw_o_QU['Q', "Unc", oi]) * 100,#
               code = 3, length = 0.03, angle = 90, col = med_cols[oi])        #
                                                                               #
        # Plotting Pol Magnitude Circle                                        #
        curve(sqrt((mw_o_PX['P', "Data", oi] * 100)^2 - x^2),                  #
              from = - mw_o_PX['P', "Data", oi] * 100,                         #
              to = mw_o_PX['P', "Data", oi] * 100, add = TRUE,                 #
              col = med_cols[oi], lty = 1, lwd = 1.5)                          #
        curve(-sqrt((mw_o_PX['P', "Data", oi] * 100)^2 - x^2),                 #
              from = - mw_o_PX['P', "Data", oi] * 100,                         #
              to = mw_o_PX['P', "Data", oi] * 100, add = TRUE,                 #
              col = med_cols[oi], lty = 1, lwd = 1.5)                          #
                                                                               #
        for(n in (oi + 1):N){                                                  #
          o_ind <- which(mwpolf[,"offset"] == n, arr.ind = T)                  #
                                                                               #
          if(length(o_ind) != 0){                                              #
            # Stars in offset n                                                #
            points(mwpolf[o_ind,'Q.Data'] * 100, mwpolf[o_ind,'U.Data'] * 100, #
                   pch = 20, cex = 1, col = pt_cols[n])                        #
            arrows(x0 = mwpolf[o_ind, 'Q.Data'] * 100,                         #
                   y0 = (mwpolf[o_ind,'U.Data'] - mwpolf[o_ind,'U.Unc']) * 100,#
                   y1 = (mwpolf[o_ind,'U.Data'] + mwpolf[o_ind,'U.Unc']) * 100,#
                   code = 3, col = pt_cols[n], length = 0.02, angle = 90)      #
            arrows(x0 = (mwpolf[o_ind,'Q.Data'] - mwpolf[o_ind,'Q.Unc']) * 100,#
                   y0 = mwpolf[o_ind, 'U.Data'] * 100,                         #
                   x1 = (mwpolf[o_ind,'Q.Data'] + mwpolf[o_ind,'Q.Unc']) * 100,#
                   code = 3, col = pt_cols[n], length = 0.02, angle = 90)      #
                                                                               #
            # Plotting Pol Magnitude Circle                                    #
            for(ns in o_ind){                                                  #
              curve(sqrt((mwpolPf[ns] * 100)^2 - x^2),                         #
                    from = -mwpolPf[ns] * 100, to = mwpolPf[ns] * 100,         #
                    add = TRUE, col = pt_cols[n], lty = 3, lwd = 1)            #
              curve(-sqrt((mwpolPf[ns] * 100)^2 - x^2),                        #
                    from = -mwpolPf[ns] * 100, to = mwpolPf[ns] * 100,         #
                    add = TRUE, col = pt_cols[n], lty = 3, lwd = 1)            #
            }                                                                  #
                                                                               #
            # Median of stars in offset n                                      #
            points(mw_o_QU['Q', "Data",n] * 100, mw_o_QU['U', "Data",n] * 100, #
                   pch = 15, cex = 1.5, col = med_cols[n])                     #
            arrows(x0 = mw_o_QU['Q', "Data", n] * 100,                         #
                   y0 = (mw_o_QU['U', "Data",n] - mw_o_QU['U', "Unc",n]) * 100,#
                   y1 = (mw_o_QU['U', "Data",n] + mw_o_QU['U', "Unc",n]) * 100,#
                   code = 3, length = 0.03, angle = 90, col = med_cols[n])     #
            arrows(x0 = (mw_o_QU['Q', "Data",n] - mw_o_QU['Q', "Unc",n]) * 100,#
                   y0 = mw_o_QU['U', "Data", n] * 100,                         #
                   x1 = (mw_o_QU['Q', "Data",n] + mw_o_QU['Q', "Unc",n]) * 100,#
                   code = 3, length = 0.03, angle = 90, col = med_cols[n])     #
                                                                               #
            # Plotting Pol Magnitude Circle                                    #
            curve(sqrt((mw_o_PX['P', "Data", n] * 100)^2 - x^2),               #
                  from = - mw_o_PX['P', "Data", n] * 100,                      #
                  to = mw_o_PX['P', "Data", n] * 100, add = TRUE,              #
                  col = med_cols[n], lty = 1, lwd = 1.5)                       #
            curve(-sqrt((mw_o_PX['P', "Data", n] * 100)^2 - x^2),              #
                  from = - mw_o_PX['P', "Data", n] * 100,                      #
                  to = mw_o_PX['P', "Data", n] * 100, add = TRUE,              #
                  col = med_cols[n], lty = 1, lwd = 1.5)                       #
          }                                                                    #
        }                                                                      #
                                                                               #
        # Median of all stars                                                  #
        points(mw_QUf['Q', "Data"] * 100, mw_QUf['U', "Data"] * 100, pch = 14, #
               cex = 2, col = "red")                                           #
        arrows(x0 = mw_QUf['Q', "Data"] * 100,                                 #
               y0 = (mw_QUf['U', "Data"] - mw_QUf['U', "Unc"]) * 100,          #
               y1 = (mw_QUf['U', "Data"] + mw_QUf['U', "Unc"]) * 100, code = 3,#
               length = 0.03, angle = 90, col = "red")                         #
        arrows(x0 = (mw_QUf['Q', "Data"] - mw_QUf['Q', "Unc"]) * 100,          #
               y0 = mw_QUf['U', "Data"] * 100,                                 #
               x1 = (mw_QUf['Q', "Data"] + mw_QUf['Q', "Unc"]) * 100, code = 3,#
               length = 0.03, angle = 90, col = "red")                         #
                                                                               #
        # Plotting Pol Magnitude Circle                                        #
        curve(sqrt((mw_Pf['P', "Data"] * 100)^2 - x^2),                        #
              from = - mw_Pf['P', "Data"] * 100, to = mw_Pf['P', "Data"] * 100,#
              add = TRUE, col = "red", lty = 1, lwd = 2)                       #
        curve(-sqrt((mw_Pf['P', "Data"] * 100)^2 - x^2),                       #
              from = - mw_Pf['P', "Data"] * 100, to = mw_Pf['P', "Data"] * 100,#
              add = TRUE, col = "red", lty = 1, lwd = 2)                       #
                                                                               #
        # Legend to points                                                     #
        legend(axlim[2] * 0.5, axlim[2], legend = leg_labs, pch = leg_pchs,    #
               col = leg_cols, cex = 2)                                        #
      }                                                                        #
    }                                                                          #
    dev.off()                                                                  #
                                                                               #
    pdf_mw_simp <- paste0(mwPol_path, basename, "-Merged-QvsU_", typename,     #
                          "_simple.pdf")                                       #
                                                                               #
    if(length(outP_inds) != 0){                                                #
      pdf(pdf_mw_simp, width = 31, height = 15)                                #
      par(mfrow = c(1, 2))                                                     #
    }else{                                                                     #
      pdf(pdf_mw_simp, width = 15, height = 15)                                #
    }                                                                          #
    rm(pdf_mw_simp)                                                            #
                                                                               #
    med_cols <- c("cornflowerblue", "darkturquoise", "darkgreen", "darkcyan",  #
                  "limegreen")                                                 #
    pt_cols <- "black"                                                         #
                                                                               #
    print(paste0("** Plotting all filtered sources with uncertainty bars for ",#
                 "type ", typename, " **"))                                    #
    if(length(outP_inds) != 0){                                                #
      ## Plotting filtered sources with uncertainty bars                       #
      limmin <- min(c(mwpolf[,'Q.Data'] - mwpolf[,'Q.Unc'],                    #
                      -(mwpolf[,'Q.Data'] + mwpolf[,'Q.Unc']),                 #
                      mwpolf[,'U.Data'] - mwpolf[,'U.Unc'],                    #
                      -(mwpolf[,'U.Data'] + mwpolf[,'U.Unc']), -max(mwpolPf))) #
      limmax <- max(c(mwpolf[,'Q.Data'] + mwpolf[,'Q.Unc'],                    #
                      -(mwpolf[,'Q.Data'] - mwpolf[,'Q.Unc']),                 #
                      mwpolf[,'U.Data'] + mwpolf[,'U.Unc'],                    #
                      -(mwpolf[,'U.Data'] - mwpolf[,'U.Unc']), max(mwpolPf)))  #
                                                                               #
      limabs <- max(abs(c(limmin, limmax)))                                    #
      axlim <- c(-limabs, limabs) * 100                                        #
      rm(limmin, limmax, limabs)                                               #
                                                                               #
      # Filtered stars                                                         #
      plot(mwpolf[,'Q.Data'] * 100, mwpolf[,'U.Data'] * 100,                   #
           xlab = 'q (%)', ylab = 'u (%)', xlim = axlim, ylim = axlim,         #
           pch = 20, cex = 1, col = pt_cols, cex.axis = 2,                     #
           main = "q vs u of MW stars (non-Intrs-Pol) with Unc bars")          #
      arrows(x0 = mwpolf[,'Q.Data'] * 100,                                     #
             y0 = (mwpolf[,'U.Data'] - mwpolf[,'U.Unc']) * 100,                #
             y1 = (mwpolf[,'U.Data'] + mwpolf[,'U.Unc']) * 100,                #
             code = 3, col = pt_cols, length = 0.02, angle = 90)               #
      arrows(x0 = (mwpolf[,'Q.Data'] - mwpolf[,'Q.Unc']) * 100,                #
             y0 = mwpolf[,'U.Data'] * 100,                                     #
             x1 = (mwpolf[,'Q.Data'] + mwpolf[,'Q.Unc']) * 100,                #
             code = 3, col = pt_cols, length = 0.02, angle = 90)               #
                                                                               #
      # Median of all stars                                                    #
      points(mw_QUf['Q', "Data"] * 100, mw_QUf['U', "Data"] * 100, pch = 14,   #
             cex = 2, col = "red")                                             #
      arrows(x0 = mw_QUf['Q', "Data"] * 100,                                   #
             y0 = (mw_QUf['U', "Data"] - mw_QUf['U', "Unc"]) * 100,            #
             y1 = (mw_QUf['U', "Data"] + mw_QUf['U', "Unc"]) * 100, code = 3,  #
             length = 0.03, angle = 90, col = "red")                           #
      arrows(x0 = (mw_QUf['Q', "Data"] - mw_QUf['Q', "Unc"]) * 100,            #
             y0 = mw_QUf['U', "Data"] * 100,                                   #
             x1 = (mw_QUf['Q', "Data"] + mw_QUf['Q', "Unc"]) * 100, code = 3,  #
             length = 0.03, angle = 90, col = "red")                           #
                                                                               #
      # Plotting Pol Magnitude Circle                                          #
      curve(sqrt((mw_Pf['P', "Data"] * 100)^2 - x^2),                          #
            from = - mw_Pf['P', "Data"] * 100, to = mw_Pf['P', "Data"] * 100,  #
            add = TRUE, col = "red", lty = 1, lwd = 2)                         #
      curve(-sqrt((mw_Pf['P', "Data"] * 100)^2 - x^2),                         #
            from = - mw_Pf['P', "Data"] * 100, to = mw_Pf['P', "Data"] * 100,  #
            add = TRUE, col = "red", lty = 1, lwd = 2)                         #
                                                                               #
      # Setting up legend                                                      #
      leg_labs <- c("Stars", "Median of Field")                                #
      leg_pchs <- c(os_pch, m_pch)                                             #
      leg_cols <- c(pt_cols, m_col)                                            #
                                                                               #
      # Legend to points                                                       #
      legend(axlim[2] * 0.5, axlim[2], legend = leg_labs, pch = leg_pchs,      #
             col = leg_cols, cex = 2)                                          #
                                                                               #
      print(paste0("** Plotting all filtered sources without uncertainty bars",#
                   " for type ", typename, " **"))                             #
      ## Plotting filtered sources without uncertainty bars                    #
      # Filtered stars                                                         #
      plot(mwpolf[,'Q.Data'] * 100, mwpolf[,'U.Data'] * 100,                   #
           xlab = 'q (%)', ylab = 'u (%)', xlim = axlim, ylim = axlim,         #
           pch = 20, cex = 1, col = pt_cols, cex.axis = 2,                     #
           main = "q vs u of MW stars (non-Intrs-Pol) without Unc bars")       #
                                                                               #
      # Median of all stars                                                    #
      points(mw_QUf['Q', "Data"] * 100, mw_QUf['U', "Data"] * 100, pch = 14,   #
             cex = 2, col = "red")                                             #
      arrows(x0 = mw_QUf['Q', "Data"] * 100,                                   #
             y0 = (mw_QUf['U', "Data"] - mw_QUf['U', "Unc"]) * 100,            #
             y1 = (mw_QUf['U', "Data"] + mw_QUf['U', "Unc"]) * 100, code = 3,  #
             length = 0.03, angle = 90, col = "red")                           #
      arrows(x0 = (mw_QUf['Q', "Data"] - mw_QUf['Q', "Unc"]) * 100,            #
             y0 = mw_QUf['U', "Data"] * 100,                                   #
             x1 = (mw_QUf['Q', "Data"] + mw_QUf['Q', "Unc"]) * 100, code = 3,  #
             length = 0.03, angle = 90, col = "red")                           #
      rm(mwpolf, mw_QUf, mw_Pf)                                                #
                                                                               #
      # Legend to points                                                       #
      legend(axlim[2] * 0.5, axlim[2], legend = leg_labs, pch = leg_pchs,      #
             col = leg_cols, cex = 2)                                          #
    }                                                                          #
    dev.off()                                                                  #
    rm(axlim, mwpol, mw_QU, mw_PX, leg_labs, leg_pchs, leg_cols)               #
  }                                                                            #
}                                                                              #
################################################################################

################################################################################
#-_----------------------------------_---_-----------------------------------_-#
#/ \--------------------------------/ \|/ \---------------------------------/ \#
#\_/--------------------------------\_/|\_/---------------------------------\_/#
#------------------------------------------------------------------------------#
################################################################################

################################################################################
############################ STOKES FILE PROCESSOR #############################
################################################################################
# From FORS2 User Manual Tab. 4.7, values in deg                               #
X_err <- array(NA, dim = nB, dimnames = list(band_str))                        #
X_err[1] <- -2.07 * pi / 180                                                   #
X_err[2] <- +1.54 * pi / 180                                                   #
X_err[3] <- +1.80 * pi / 180                                                   #
X_err[4] <- -1.19 * pi / 180                                                   #
X_err[5] <- -2.89 * pi / 180                                                   #
uX_err <- 0.005 * pi / 180                                                     #
                                                                               #
# Will be used within the next main loop                                       #
delt <- dimsObj[1] * .8 / dimsObj[2]                                           #
                                                                               #
print("### Merging all available offsets of targets Q, U, P and X maps ###")   #
                                                                               #
for(this in 1:length(Ipaths_labels)){                                          #
  ############################ Listing input files ########################### #
  ##############################################################################
  switch(this,                                                               # #
         this_folder <- obspol_folder, this_folder <- sky_instCpol_folder,   # #
         this_folder <- skypol_folder)                                       # #
                                                                             # #
  this_I <- Ipaths_labels[this]                                              # #
  tin_type <- file_type_in[this]                                             # #
  t_type <- file_type_out[this]                                              # #
  t_label <- file_label[this]                                                # #
  tout_fold <- folder_t[this]                                                # #
                                                                             # #
  print(paste0("## Processing files of type ", t_label, " ##"))              # #
                                                                             # #
  pathList <- list.files(this_folder, full.names = T)                        # #
  fileList <- list.files(this_folder, full.names = F)                        # #
                                                                             # #
  # Exclusion of non-fits files from the previous lists                      # #
  is_fits <- grep(".fits", pathList)                                         # #
                                                                             # #
  if(length(which(is_fits == FALSE, arr.ind = TRUE)) != 0){                  # #
    pathList <- pathList[is_fits]                                            # #
    fileList <- fileList[is_fits]                                            # #
  }                                                                          # #
  rm(is_fits)                                                                # #
                                                                             # #
  fitsCount <- length(fileList)                                              # #
                                                                             # #
  # Determining type of fits files present                                   # #
  # Possible types 1 of fits files must match 'file_types_1'                 # #
  # Possible types 2 of fits files must can be:                              # #
  # "STOKES Q", "STOKES U"                                                   # #
  types_arr <- get_fits_header_list_str(pathList, TYP1)                      # #
  stokes_arr <- get_fits_header_list_str(pathList, TYP2)                     # #
  auth_arr <- get_fits_header_list_str(pathList, AUTHOR)                     # #
                                                                             # #
  print("* Checking for unexpected type files *")                            # #
                                                                             # #
  for(f in 1:fitsCount){                                                     # #
                                                                             # #
    if((types_arr[f] != tin_type)){                                          # #
      types_arr[f] <- "NOT EXP"                                              # #
    }                                                                        # #
    if((stokes_arr[f] != 'STOKES Q' && stokes_arr[f] != 'STOKES U')){        # #
      stokes_arr[f] <- "NOT STOKES"                                          # #
    }                                                                        # #
                                                                             # #
    if(auth_arr[f] != "JRS2023"){                                            # #
      auth_arr[f] <- "NOT EXP"                                               # #
    }                                                                        # #
  }                                                                          # #
  rm(fitsCount, f)                                                           # #
                                                                             # #
  print("* Indexing stokes files *")                                         # #
                                                                             # #
  Q_ind <- as.vector(which(stokes_arr == "STOKES Q" & auth_arr == "JRS2023" &# #
                             types_arr == tin_type, arr.ind = TRUE))         # #
  qCount <- length(Q_ind)                                                    # #
                                                                             # #
  U_ind <- as.vector(which(stokes_arr == "STOKES U" & auth_arr == "JRS2023" &# #
                             types_arr == tin_type, arr.ind = TRUE))         # #
  uCount <- length(U_ind)                                                    # #
  rm(stokes_arr, auth_arr, types_arr)                                        # #
                                                                             # #
  print("* Checking if Stokes Q and U files are all paired *")               # #
                                                                             # #
  err_no_pair <- "At least one file is not paired."                          # #
  stopifnot(err_no_pair = qCount == uCount)                                  # #
  print("* Stokes Q and U are all paired *")                                 # #
                                                                             # #
  tCount <- qCount                                                           # #
  rm(qCount, uCount, err_no_pair)                                            # #
                                                                             # #
  print("* Tagging Stokes Q and U files by band filter *")                   # #
  tagDim <- c(2, tCount)                                                     # #
  tagDimNames <- list(c("Q", "U"), NULL)                                     # #
                                                                             # #
  bands_arr <- array(NA, dim = tagDim, dimnames = tagDimNames)               # #
  bands_arr["Q",] <- get_fits_header_list_str(pathList[Q_ind], FILTER)       # #
  bands_arr["U",] <- get_fits_header_list_str(pathList[U_ind], FILTER)       # #
                                                                             # #
  numbs_arr <- array(0, dim = tagDim, dimnames = tagDimNames)                # #
  rm(tagDim, tagDimNames)                                                    # #
                                                                             # #
  print("* Listing Stokes Q and U files *")                                  # #
  for(t in 1:tCount){                                                        # #
    name_q <- fileList[Q_ind[t]]                                             # #
    name_u <- fileList[U_ind[t]]                                             # #
                                                                             # #
    numbs_arr["Q", t] <- as.numeric(                                         # #
      str_extract(str_extract(name_q, "#[0-9]+"), "[0-9]+"))                 # #
    numbs_arr["U", t] <- as.numeric(                                         # #
      str_extract(str_extract(name_u, "#[0-9]+"), "[0-9]+"))                 # #
                                                                             # #
    print(paste0("** Stokes Q #",t, " is in ", bands_arr["Q",t], " band **"))# #
    print(paste0("** Stokes U #",t, " is in ", bands_arr["U",t], " band **"))# #
    rm(name_q, name_u)                                                       # #
  }                                                                          # #
  rm(t, tCount)                                                              # #
                                                                             # #
  print("* Creating Stokes Q and U files band filter indexes *")             # #
  max_N <- array(0, dim = nB, dimnames = list(band_str))                     # #
  for(b in bands){                                                           # #
    bs <- band_str[b]                                                        # #
                                                                             # #
    t_U_ind <- which(bands_arr["U",] == bs, arr.ind = TRUE)                  # #
    t_Q_ind <- which(bands_arr["Q",] == bs, arr.ind = TRUE)                  # #
                                                                             # #
    switch(b,                                                                # #
           c(u_ind_U <- t_U_ind, u_ind_Q <- t_Q_ind),                        # #
           c(b_ind_U <- t_U_ind, b_ind_Q <- t_Q_ind),                        # #
           c(v_ind_U <- t_U_ind, v_ind_Q <- t_Q_ind),                        # #
           c(r_ind_U <- t_U_ind, r_ind_Q <- t_Q_ind),                        # #
           c(i_ind_U <- t_U_ind, i_ind_Q <- t_Q_ind))                        # #
                                                                             # #
    print(paste0("** Checking how many measurement sets there are in ", bs,  # #
                 " band **"))                                                # #
    max_N[bs] <- length(numbs_arr["U", t_U_ind])                             # #
                                                                             # #
    if(is.infinite(max_N[bs]) || max_N[bs] == 0){                            # #
      max_N[bs] <- 0                                                         # #
      bands <- bands[-which(bands == b, arr.ind = T)]                        # #
                                                                             # #
      switch(b,                                                              # #
             rm(u_ind_U, u_ind_Q), rm(b_ind_U, b_ind_Q), rm(v_ind_U,v_ind_Q),# #
             rm(r_ind_U, r_ind_Q), rm(i_ind_U, i_ind_Q))                     # #
    }                                                                        # #
  }                                                                          # #
  rm(b, bands_arr)                                                           # #
  ##############################################################################
                                                                               #
  for(b in bands){                                                             #
    bs <- band_str[b]                                                          #
                                                                               #
    ############################# Loading and Merging ######################## #
    ############################################################################
    print(paste0("* Starting to process ", t_label, " files in ", bs, " ban",# #
                 "d *"))                                                     # #
                                                                             # #
    print("** Loading relevant band filter file indexes **")                 # #
    switch(b,                                                                # #
           c(band_ind_Q <- u_ind_Q, band_ind_U <- u_ind_U,                   # #
             rm(u_ind_Q, u_ind_U)),                                          # #
           c(band_ind_Q <- b_ind_Q, band_ind_U <- b_ind_U,                   # #
             rm(b_ind_Q, b_ind_U)),                                          # #
           c(band_ind_Q <- v_ind_Q, band_ind_U <- v_ind_U,                   # #
             rm(v_ind_Q, v_ind_U)),                                          # #
           c(band_ind_Q <- r_ind_Q, band_ind_U <- r_ind_U,                   # #
             rm(r_ind_Q, r_ind_U)),                                          # #
           c(band_ind_Q <- i_ind_Q, band_ind_U <- i_ind_U,                   # #
             rm(i_ind_Q, i_ind_U)))                                          # #
                                                                             # #
    c2X <- cos(2 * X_err[b])                                                 # #
    s2X <- sin(2 * X_err[b])                                                 # #
    u_c2X <- sqrt((2 * s2X * uX_err)^2)                                      # #
    u_s2X <- sqrt((2 * c2X * uX_err)^2)                                      # #
                                                                             # #
    N <- max_N[bs]                                                           # #
                                                                             # #
    temp_dim <- c(dimsObj, length(munc_par), N)                              # #
    temp_dname <- list(NULL, NULL, c("Data", "Unc"), NULL)                   # #
                                                                             # #
    temp_Q <- array(NA, dim = temp_dim, dimnames = temp_dname)               # #
    temp_U <- array(NA, dim = temp_dim, dimnames = temp_dname)               # #
    hold_crpix_N <- array(NA, dim = c(N, 3))                                 # #
                                                                             # #
    if(N != 1){                                                              # #
      # Iterating over the different offset runs at the same band filter     # #
      for(n in 1:N){                                                         # #
        no <- numbs_arr["U", band_ind_U[n]]                                  # #
                                                                             # #
        print(paste0("** Starting to process ", t_label, " file set #", n,   # #
                     "/", N, " in ", bs, " band **"))                        # #
                                                                             # #
        print("*** Loading relevant set number file indexes ***")            # #
        numb_ind_Q <- which(numbs_arr["Q",] == no, arr.ind = TRUE)           # #
        numb_ind_U <- which(numbs_arr["U",] == no, arr.ind = TRUE)           # #
                                                                             # #
        # Intersect indexes of band with set_numb                            # #
        print(paste0("*** Getting relevant indexes from intersection of set",# #
                     " number indexes with band indexes ***"))               # #
        temp_ind_Q <- intersect(band_ind_Q, numb_ind_Q)                      # #
        length_ind_Q <- length(temp_ind_Q)                                   # #
        temp_ind_U <- intersect(band_ind_U, numb_ind_U)                      # #
        length_ind_U <- length(temp_ind_U)                                   # #
                                                                             # #
        print("*** Checking that there is 1 Q and U files ***")              # #
        if(length_ind_Q != 1 || length_ind_U != 1){                          # #
          print(paste0("#SKIP! Expected 2 data files, 1 Q and 1 U, ",        # #
                       length_ind_Q, " were found for Q, and ", length_ind_U,# #
                       " were found for U instead for the filter & number s",# #
                       "ubset: band = ", bs, ", number = ", n))              # #
          next                                                               # #
        }                                                                    # #
        rm(length_ind_Q, length_ind_U)                                       # #
                                                                             # #
        ######################### Loading Relevant Files ################### # #
        ########################################################################
        if(n == N){                                                        # # #
          basename <- strsplit(fileList[Q_ind[temp_ind_Q]], "_#")[[1]][1]  # # #
        }                                                                  # # #
                                                                           # # #
        print(paste0("*** Loading Q from set #", n, "/", N, " in ", bs,    # # #
                     " band ***"))                                         # # #
        Q_fits <- readFITS(pathList[Q_ind[temp_ind_Q]])                    # # #
        Q_data <- Q_fits$imDat                                             # # #
        Q_head <- Q_fits$header[-c(8:10)]                                  # # #
        rm(temp_ind_Q, Q_fits)                                             # # #
                                                                           # # #
        print(paste0("*** Loading U from set #", n, "/", N, " in ", bs,    # # #
                     " band ***"))                                         # # #
        U_fits <- readFITS(pathList[U_ind[temp_ind_U]])                    # # #
        U_data <- U_fits$imDat                                             # # #
        U_head <- U_fits$header[-c(8:10)]                                  # # #
        rm(temp_ind_U, U_fits)                                             # # #
        ########################################################################
                                                                             # #
        #------------------------------------------------------------------# # #
                                                                             # #
        ################## Placing Offsets in holder frame ################# # #
        ########################################################################
        refx <- get_fits_header_num(Q_head, REFX)                          # # #
        refy <- get_fits_header_num(Q_head, REFY)                          # # #
        refz <- 1                                                          # # #
        valx <- get_fits_header_num(Q_head, VALX)                          # # #
        valy <- get_fits_header_num(Q_head, VALY)                          # # #
        valz <- NA                                                         # # #
        typex <- "  'PIXEL     '             / Coordinate system of x-axis"# # #
        typey <- "  'PIXEL     '             / Coordinate system of y-axis"# # #
        typez <- "  'DATA & UNC'             / Coordinate system of z-axis"# # #
        crpix <- c(refx, refy, refz)                                       # # #
        crval <- c(valx, valy, valz)                                       # # #
        ctype <- c(typex, typey, typez)                                    # # #
        binx <- get_fits_header_num(Q_head, BINX)                          # # #
        biny <- get_fits_header_num(Q_head, BINY)                          # # #
        pxscl <- get_fits_header_num(Q_head, PIXSCALE)                     # # #
        rm(refx, refy, refz, valx, valy, valz, typex, typey, typez)        # # #
                                                                           # # #
        Q_head <- delKwv('SRC_1', Q_head)                                  # # #
        Q_head <- delKwv('SRC_2', Q_head)                                  # # #
        Q_head <- delKwv('SRC_3', Q_head)                                  # # #
        Q_head <- delKwv('SRC_4', Q_head)                                  # # #
        Q_head <- delKwv('SRC_5', Q_head)                                  # # #
        Q_head <- delKwv('SRC_6', Q_head)                                  # # #
        Q_head <- delKwv('SRC_7', Q_head)                                  # # #
        Q_head <- delKwv('SRC_8', Q_head)                                  # # #
        Q_head <- modVal('EXTNAME', 'MERGED', "Extension name", Q_head)    # # #
                                                                           # # #
        U_head <- delKwv('SRC_1', U_head)                                  # # #
        U_head <- delKwv('SRC_2', U_head)                                  # # #
        U_head <- delKwv('SRC_3', U_head)                                  # # #
        U_head <- delKwv('SRC_4', U_head)                                  # # #
        U_head <- delKwv('SRC_5', U_head)                                  # # #
        U_head <- delKwv('SRC_6', U_head)                                  # # #
        U_head <- delKwv('SRC_7', U_head)                                  # # #
        U_head <- delKwv('SRC_8', U_head)                                  # # #
        U_head <- modVal('EXTNAME', 'MERGED', "Extension name", U_head)    # # #
                                                                           # # #
        if(n == 1){                                                        # # #
          ref_mrg <- array(NA, dim = c(7, N),                              # # #
                           dimnames = list(c("x", "y", "ra", "dec",        # # #
                                             "scl_x", "scl_y", "scl_flg"), # # #
                                           NULL))                          # # #
          hold_crval <- crval                                              # # #
          hold_ctype <- ctype                                              # # #
        }                                                                  # # #
                                                                           # # #
        hold_crpix_N[n, ] <- crpix                                         # # #
                                                                           # # #
        ref_mrg["x", n] <- crpix[1]                                        # # #
        ref_mrg["y", n] <- crpix[2]                                        # # #
        ref_mrg["ra", n] <- get_fits_header_num(Q_head, REFRA)             # # #
        ref_mrg["dec", n] <- get_fits_header_num(Q_head, REFDEC)           # # #
        ref_mrg["scl_x", n] <- pxscl * binx                                # # #
        ref_mrg["scl_y", n] <- pxscl * biny                                # # #
        ref_mrg["scl_flg", n] <- 1                                         # # #
                                                                           # # #
        if(n == 1){                                                        # # #
          scl_x <- ref_mrg["scl_x", 1] / 3600                              # # #
          scl_y <- ref_mrg["scl_y", 1] / 3600                              # # #
        }                                                                  # # #
                                                                           # # #
        if(n != 1 && (ref_mrg["scl_x", n] != ref_mrg["scl_x", 1] ||        # # #
                      ref_mrg["scl_y", n] != ref_mrg["scl_y", 1])){        # # #
          ref_mrg["scl_flg", n] <- 0                                       # # #
          print(paste0("WARNING: Offset #", n, " has a different scale t", # # #
                       "han that of Offset #1, and as such it will not me",# # #
                       "rged with the remaining offsets."))                # # #
        }                                                                  # # #
        rm(binx, biny, pxscl)                                              # # #
                                                                           # # #
        if(n == 1){                                                        # # #
          temp_Q[,, "Data", n] <- Q_data[,,1]                              # # #
          temp_Q[,, "Unc", n] <- Q_data[,,2]                               # # #
          temp_U[,, "Data", n] <- U_data[,,1]                              # # #
          temp_U[,, "Unc", n] <- U_data[,,2]                               # # #
        }else{                                                             # # #
          if(ref_mrg["scl_flg", n]){                                       # # #
            dx <- ref_mrg["x", 1] - ref_mrg["x", n]                        # # #
            dy <- ref_mrg["y", 1] - ref_mrg["y", n]                        # # #
                                                                           # # #
            dra_x <- (ref_mrg["ra",1] - ref_mrg["ra",n]) * cosd[b] / scl_x # # #
            ddec_y <- (ref_mrg["dec",1] - ref_mrg["dec",n]) / scl_y        # # #
            difx <- round(dx + dra_x)                                      # # #
            dify <- round(dy + ddec_y)                                     # # #
            rm(dx, dy, dra_x, ddec_y)                                      # # #
                                                                           # # #
            if(difx >= 0){                                                 # # #
              hold_crpix_N[n, 1] <- hold_crpix_N[n, 1] + difx              # # #
              hold_xs <- 1:(dimsObj[1] - difx)                             # # #
              off_xs <- (1 + difx):dimsObj[1]                              # # #
            }else{                                                         # # #
              hold_crpix_N[n, 1] <- hold_crpix_N[n, 1] - difx              # # #
              hold_xs <- (1 - difx):dimsObj[1]                             # # #
              off_xs <- 1:(dimsObj[1] + difx)                              # # #
            }                                                              # # #
            if(dify <= 0){                                                 # # #
              hold_crpix_N[n, 2] <- hold_crpix_N[n, 2] + dify              # # #
              hold_ys <- (1 - dify):dimsObj[2]                             # # #
              off_ys <- 1:(dimsObj[2] + dify)                              # # #
            }else{                                                         # # #
              hold_crpix_N[n, 2] <- hold_crpix_N[n, 2] - dify              # # #
              hold_ys <- 1:(dimsObj[2] - dify)                             # # #
              off_ys <- (1 + dify):dimsObj[2]                              # # #
            }                                                              # # #
            rm(difx, dify)                                                 # # #
          }                                                                # # #
                                                                           # # #
          temp_Q[hold_xs, hold_ys,, n] <- Q_data[off_xs, off_ys, 1:2]      # # #
                                                                           # # #
          print(paste0("*** Q ", t_label," from set #", n, "/", N, " in ", # # #
                       bs, " band has been loaded to holder frame ***"))   # # #
                                                                           # # #
          temp_U[hold_xs, hold_ys,, n] <- U_data[off_xs, off_ys, 1:2]      # # #
                                                                           # # #
          print(paste0("*** U ", t_label," from set #", n, "/", N, " in ", # # #
                       bs, " band has been loaded to holder frame ***"))   # # #
        }                                                                  # # #
      }                                                                    # # #
      ##########################################################################
    }else{                                                                   # #
      print(paste0("** Starting to process ", t_label, " files in ", bs,     # #
                   " band **"))                                              # #
                                                                             # #
      # Intersect indexes of band with set_numb                              # #
      print(paste0("*** Getting relevant band indexes ***"))                 # #
      temp_ind_Q <- band_ind_Q                                               # #
      length_ind_Q <- length(temp_ind_Q)                                     # #
      temp_ind_U <- band_ind_U                                               # #
      length_ind_U <- length(temp_ind_U)                                     # #
                                                                             # #
      print("*** Checking that there is 1 Q and U files ***")                # #
      if(length_ind_Q != 1 || length_ind_U != 1){                            # #
        print(paste0("#SKIP! Expected 2 data files, 1 Q and 1 U, ",          # #
                     length_ind_Q, " were found for Q, and ", length_ind_U,  # #
                     "were found for U instead for the filter & number su",  # #
                     "bset: band = ", bs, ", number = ", n))                 # #
        next                                                                 # #
      }                                                                      # #
      rm(length_ind_Q, length_ind_U)                                         # #
                                                                             # #
      ######################### Loading Relevant Files ##################### # #
      ##########################################################################
      basename <- strsplit(fileList[Q_ind[temp_ind_Q]], "_#")[[1]][1]      # # #
                                                                           # # #
      print(paste0("*** Loading Q in ", bs, " band ***"))                  # # #
      Q_fits <- readFITS(pathList[Q_ind[temp_ind_Q]])                      # # #
      temp_Q[,,"Data", 1] <- Q_fits$imDat[,,1]                             # # #
      temp_Q[,,"Unc", 1] <- Q_fits$imDat[,,2]                              # # #
      Q_head <- Q_fits$header[-c(8:10)]                                    # # #
      rm(temp_ind_Q, Q_fits)                                               # # #
                                                                           # # #
      print(paste0("*** Loading U in ", bs, " band ***"))                  # # #
      U_fits <- readFITS(pathList[U_ind[temp_ind_U]])                      # # #
      temp_U[,,"Data", 1] <- U_fits$imDat[,,1]                             # # #
      temp_U[,,"Unc", 1] <- U_fits$imDat[,,2]                              # # #
      U_head <- U_fits$header[-c(8:10)]                                    # # #
      rm(temp_ind_U, U_fits)                                               # # #
      ##########################################################################
                                                                             # #
      #--------------------------------------------------------------------# # #
                                                                             # #
      ####################### Setting up header info ####################### # #
      ##########################################################################
      refx <- get_fits_header_num(Q_head, REFX)                            # # #
      refy <- get_fits_header_num(Q_head, REFY)                            # # #
      refz <- 1                                                            # # #
      refra <- get_fits_header_num(Q_head, REFRA)                          # # #
      refdec <- get_fits_header_num(Q_head, REFDEC)                        # # #
      valx <- get_fits_header_num(Q_head, VALX)                            # # #
      valy <- get_fits_header_num(Q_head, VALY)                            # # #
      valz <- NA                                                           # # #
      typex <- "  'PIXEL     '              / Coordinate system of x-axis "# # #
      typey <- "  'PIXEL     '              / Coordinate system of y-axis "# # #
      typez <- "  'DATA & UNC'              / Coordinate system of z-axis "# # #
      hold_crpix_N[1,] <- c(refx, refy, refz)                              # # #
      hold_crval <- c(valx, valy, valz)                                    # # #
      hold_ctype <- c(typex, typey, typez)                                 # # #
                                                                           # # #
      scl_x <- pxscl * binx / 3600                                         # # #
      scl_y <- pxscl * biny / 3600                                         # # #
      rm(refx, refy, refz, valx, valy, valz, typex, typey, typez)          # # #
                                                                           # # #
      Q_head <- delKwv('SRC_1', Q_head)                                    # # #
      Q_head <- delKwv('SRC_2', Q_head)                                    # # #
      Q_head <- delKwv('SRC_3', Q_head)                                    # # #
      Q_head <- delKwv('SRC_4', Q_head)                                    # # #
      Q_head <- delKwv('SRC_5', Q_head)                                    # # #
      Q_head <- delKwv('SRC_6', Q_head)                                    # # #
      Q_head <- delKwv('SRC_7', Q_head)                                    # # #
      Q_head <- delKwv('SRC_8', Q_head)                                    # # #
      Q_head <- modVal('EXTNAME', 'MERGED', "Extension name", Q_head)      # # #
                                                                           # # #
      U_head <- delKwv('SRC_1', U_head)                                    # # #
      U_head <- delKwv('SRC_2', U_head)                                    # # #
      U_head <- delKwv('SRC_3', U_head)                                    # # #
      U_head <- delKwv('SRC_4', U_head)                                    # # #
      U_head <- delKwv('SRC_5', U_head)                                    # # #
      U_head <- delKwv('SRC_6', U_head)                                    # # #
      U_head <- delKwv('SRC_7', U_head)                                    # # #
      U_head <- delKwv('SRC_8', U_head)                                    # # #
      U_head <- modVal('EXTNAME', 'MERGED', "Extension name", U_head)      # # #
      ##########################################################################
    }                                                                        # #
    ############################################################################
                                                                               #
    #------------------------------------------------------------------------# #
                                                                               #
    ###################### Binning Stacked Q and U maps ###################### #
    ############################################################################
    print(paste0("** Binning stacked Q ", t_label, " map in ",bs," band **"))# #
                                                                             # #
    if(N != 1){                                                              # #
      binDimS <- array(NA, dim = c(N - 1, 2))                                # #
    }                                                                        # #
                                                                             # #
    for(n in 1:N){                                                           # #
      if(n == 1){                                                            # #
        Qbin_temp <- bin_shrink_map(temp_Q[,,,n], binS,                      # #
                                    hold_crpix_N[n, 1:2])                    # #
                                                                             # #
        binDim <- c(dim(Qbin_temp), N)                                       # #
        binNam <- list(NULL, NULL, c("Data", "Unc"), NULL)                   # #
                                                                             # #
        Qbin <- array(NA, dim = binDim, dimnames = binNam)                   # #
                                                                             # #
        Qbin[,,,n] <- Qbin_temp                                              # #
      }else{                                                                 # #
        Qbin_temp <- bin_shrink_map(temp_Q[,,,n], binS,                      # #
                                    hold_crpix_N[1, 1:2])                    # #
        dims_temp <- dim(Qbin_temp)                                          # #
        binDimS[n - 1, ] <- dims_temp[1:2]                                   # #
                                                                             # #
        Qbin[1:dims_temp[1], 1:dims_temp[2],, n] <- Qbin_temp                # #
      }                                                                      # #
    }                                                                        # #
    rm(Qbin_temp, Qbin_copy, dims_temp)                                      # #
                                                                             # #
    print(paste0("** Binning stacked U ", t_label, " map in ",bs," band **"))# #
                                                                             # #
    Ubin <- array(NA, dim = binDim, dimnames = binNam)                       # #
                                                                             # #
    for(n in 1:N){                                                           # #
      Ubin_temp <- bin_shrink_map(temp_U[,,,n], binS, hold_crpix_N[1, 1:2])  # #
                                                                             # #
      dims_temp <- dim(Ubin_temp)                                            # #
                                                                             # #
      Ubin[1:dims_temp[1], 1:dims_temp[2],, n] <- Ubin_temp                  # #
    }                                                                        # #
    rm(Ubin_temp, Qbin_copy, Ubin_copy, dims_temp)                           # #
                                                                             # #
    temp_I <- readFITS(Stacked_Ipaths[this_I, b])$imDat                      # #
    Ibin <- array(NA, dim = binDim, dimnames = binNam)                       # #
                                                                             # #
    print(paste0("** Binning stacked I ", t_label, " map in ",bs," band **"))# #
                                                                             # #
    for(n in 1:N){                                                           # #
        Ibin_temp <- bin_code_map(temp_I[,,,n], binS, hold_crpix_N[1, 1:2])  # #
        dims_temp <- dim(Ibin_temp)                                          # #
                                                                             # #
        Ibin[1:dims_temp[1], 1:dims_temp[2],, n] <- Ibin_temp                # #
    }                                                                        # #
    rm(temp_I, Ibin_temp, dims_temp)                                         # #
    ############################################################################
                                                                               #
    if(this != 1 && this != length(Ipaths_labels)){                            #
      t_out <- array(allCPol_paths, dim = dim(allCPol_paths),                  #
                     dimnames = dimnames(allCPol_paths))                       #
    }else{                                                                     #
      t_out <- tout_fold                                                       #
    }                                                                          #
                                                                               #
    print(paste0("** Calculating ", t_label," pol maps in ", bs, " band **"))  #
                                                                               #
    ################# Adapting Q and U to specific cases ##################### #
    ############################################################################
    if(t_label == "obs" || t_label == "skyC"){                               # #
      hld_Qbin <- Qbin                                                       # #
      hld_Ubin <- Ubin                                                       # #
    }                                                                        # #
                                                                             # #
    for(n in 1:N){                                                           # #
      if(t_label == "mwC"){                                                  # #
        Qbin[,,"Data", n] <- hld_Qbin[,,"Data", n] -                         # #
          mwQU["Q","Data", b, "raw"]                                         # #
        Qbin[,,"Unc", n] <- unc_add(hld_Qbin[,,"Unc", n],                    # #
                                    mwQU["Q", "Unc", b, "raw"])              # #
        Ubin[,,"Data", n] <- hld_Ubin[,,"Data", n] -                         # #
          mwQU["U","Data", b, "raw"]                                         # #
        Ubin[,,"Unc", n] <- unc_add(hld_Ubin[,,"Unc", n],                    # #
                                 mwQU["U", "Unc", b, "raw"])                 # #
      }                                                                      # #
      if(t_label == "instC"){                                                # #
        Qbin[,,"Data", n] <- hld_Qbin[,,"Data", n] - Qinst[,, "Data", b]     # #
        Qbin[,,"Unc", n] <- unc_add(hld_Qbin[,,"Unc", n],                    # #
                                    Qinst[,, "Unc", b])                      # #
        Ubin[,,"Data", n] <- hld_Ubin[,,"Data", n] - Uinst[,, "Data", b]     # #
        Ubin[,,"Unc", n] <- unc_add(hld_Ubin[,,"Unc", n],                    # #
                                    Uinst[,, "Unc", b])                      # #
      }                                                                      # #
      if(t_label == "chromC"){                                               # #
        Qc2X <- hld_Qbin[,,"Data", n] * c2X                                  # #
        Qs2X <- hld_Qbin[,,"Data", n] * s2X                                  # #
        Uc2X <- hld_Ubin[,,"Data", n] * c2X                                  # #
        Us2X <- hld_Ubin[,,"Data", n] * s2X                                  # #
        u_Qc2X <- unc_mult(hld_Qbin[,,"Data", n], c2X,                       # #
                           hld_Qbin[,,"Unc", n], u_c2X)                      # #
        u_Qs2X <- unc_mult(hld_Qbin[,,"Data", n], s2X,                       # #
                           hld_Qbin[,,"Unc", n], u_s2X)                      # #
        u_Uc2X <- unc_mult(hld_Ubin[,,"Data"], c2X,                          # #
                           hld_Ubin[,,"Unc", n], u_c2X)                      # #
        u_Us2X <- unc_mult(hld_Ubin[,,"Data", n], s2X,                       # #
                           hld_Ubin[,,"Unc", n], u_s2X)                      # #
                                                                             # #
        Qbin[,,"Data", n] <- Qc2X + Us2X                                     # #
        Ubin[,,"Data", n] <- Uc2X + Qs2X                                     # #
        Qbin[,,"Unc", n] <- unc_add(u_Qc2X, u_Us2X)                          # #
        Ubin[,,"Unc", n] <- unc_add(u_Uc2X, u_Qs2X)                          # #
        rm(Qc2X, Qs2X, Uc2X, Us2X, u_Qc2X, u_Qs2X, u_Uc2X, u_Us2X,           # #
           hld_Qbin, hld_Ubin)                                               # #
      }                                                                      # #
      if(t_label == "allC"){                                                 # #
        hld_Qbin <- Qbin                                                     # #
        hld_Ubin <- Ubin                                                     # #
                                                                             # #
        Qbin[,,"Data", n] <- hld_Qbin[,,"Data", n] -                         # #
          mwQU["Q", "Data", b, "instC"]                                      # #
        Ubin[,,"Data", n] <- hld_Ubin[,,"Data", n] -                         # #
          mwQU["U", "Data", b, "instC"]                                      # #
        Qbin[,,"Unc", n] <- unc_add(hld_Qbin[,,"Unc", n],                    # #
                                    mwQU["Q", "Unc", b, "instC"]^2)          # #
        Ubin[,,"Unc", n] <- unc_add(hld_Ubin[,,"Unc", n],                    # #
                                    mwQU["U", "Unc", b, "instC"]^2)          # #
                                                                             # #
        Qc2X <- Qbin[,,"Data", n] * c2X                                      # #
        Qs2X <- Qbin[,,"Data", n] * s2X                                      # #
        Uc2X <- Ubin[,,"Data", n] * c2X                                      # #
        Us2X <- Ubin[,,"Data", n] * s2X                                      # #
        u_Qc2X <- unc_mult(Qbin[,,"Data", n], c2X, Qbin[,,"Unc", n], u_c2X)  # #
        u_Qs2X <- unc_mult(Qbin[,,"Data", n], s2X, Qbin[,,"Unc", n], u_s2X)  # #
        u_Uc2X <- unc_mult(Ubin[,,"Data", n], c2X, Ubin[,,"Unc", n], u_c2X)  # #
        u_Us2X <- unc_mult(Ubin[,,"Data", n], s2X, Ubin[,,"Unc", n], u_s2X)  # #
                                                                             # #
        Qbin[,,"Data", n] <- Qc2X + Us2X                                     # #
        Ubin[,,"Data", n] <- Uc2X + Qs2X                                     # #
        Qbin[,,"Unc", n] <- unc_add(u_Qc2X, u_Us2X)                          # #
        Ubin[,,"Unc", n] <- unc_add(u_Uc2X, u_Qs2X)                          # #
      }                                                                      # #
      rm(Qc2X, Qs2X, Uc2X, Us2X, u_Qc2X, u_Qs2X, u_Uc2X, u_Us2X, hld_Qbin,   # #
         hld_Ubin)                                                           # #
                                                                             # #
      Q_head <- modVal('CORR-MW ','YES', "Was MW pol. subtracted?", Q_head)  # #
      U_head <- modVal('CORR-MW ','YES', "Was MW pol. subtracted?", U_head)  # #
      Q_head <- modVal('CORR-CHR','YES',"Was pol. angle corrected?",Q_head)  # #
      U_head <- modVal('CORR-CHR','YES',"Was pol. angle corrected?",U_head)  # #
      Q_head <- modVal('CORR-PBI','YES',"Was pol. bias corrected?", Q_head)  # #
      U_head <- modVal('CORR-PBI','YES',"Was pol. bias subtracted?",U_head)  # #
      Q_head <- modVal('SOURCE', 'merge_OFFSETS.R',                          # #
                       "Script used to generate this file ", Q_head)         # #
      U_head <- modVal('SOURCE', 'merge_OFFSETS.R',                          # #
                       "Script used to generate this file ", U_head)         # #
    }                                                                        # #
    ############################################################################
                                                                               #
    # Preping outputs                                                          #
    fits_out <- paste0("_", t_label, ".fits")                                  #
                                                                               #
    if(this != 1 && this != length(Ipaths_labels)){                            #
      temp_out <- array(str_remove(t_out[,1], "Dense/"), dim = dim(t_out)[1],  #
                        dimnames = list(dimnames(t_out)[[1]]))                 #
      stats_out_path <- array(paste0(temp_out, basename, "-Merged-Stats_",     #
                                     t_label, ".csv"),                         #
                              dim = dim(temp_out), dimnames=dimnames(temp_out))#
      pdf_fold <- temp_out                                                     #
      X_col_path <- array(paste0(temp_out, basename, "-Merged_Pol_Ang_X_",     #
                                 t_label, ".pdf"),                             #
                          dim = dim(temp_out), dimnames = dimnames(temp_out))  #
      P_col_path <- array(paste0(temp_out, basename, "-Merged_Pol_Deg_P_",     #
                                 t_label, ".pdf"),                             #
                          dim = dim(temp_out), dimnames = dimnames(temp_out))  #
      rm(temp_out)                                                             #
    }else{                                                                     #
      stats_out_path <- paste0(t_out, basename, "-Merged-Stats_", t_label,     #
                               ".csv")                                         #
      pdf_fold <- t_out                                                        #
      X_col_path <- paste0(t_out, basename, "-Merged_Pol_Ang_X_", t_label,     #
                           ".pdf")                                             #
      P_col_path <- paste0(t_out, basename, "-Merged_Pol_Deg_P_", t_label,     #
                           ".pdf")                                             #
    }                                                                          #
                                                                               #
    pol_stats <- array(NA, dim = c(4, 4, 3),                                   #
                       dimnames = list(c("Q", "U", "P", "X"),                  #
                                       c("min", "max", "median", "mad"),       #
                                       c("val", "x", "y")))                    #
                                                                               #
    Pbin <- array(NA, dim = dim(Qbin), dimnames = dimnames(Qbin))              #
    IPbin <- array(NA, dim = dim(Qbin), dimnames = dimnames(Qbin))             #
    Xbin <- array(NA, dim = dim(Qbin), dimnames = dimnames(Qbin))              #
                                                                               #
    ####################### Calculating P, IP and X ########################## #
    ############################################################################
    print(paste0("** Calculating ", t_label, " P and X in ", bs, " band ",   # #
                 "**"))                                                      # #
                                                                             # #
    for(n in 1:N){                                                           # #
      if(t_label == "allC"){                                                 # #
        P <- debiased_P_from_QU(Qbin[,,"Data", n], Ubin[,,"Data", n],        # #
                                Qbin[,,"Unc", n], Ubin[,,"Unc", n])          # #
      }else{                                                                 # #
        P <- P_from_QU(Qbin[,,"Data", n], Ubin[,,"Data", n],                 # #
                       Qbin[,,"Unc", n], Ubin[,,"Unc", n])                   # #
      }                                                                      # #
                                                                             # #
      X <- X_from_QU(Qbin[,,"Data", n], Ubin[,,"Data", n],                   # #
                     Qbin[,,"Unc", n], Ubin[,,"Unc", n])                     # #
                                                                             # #
      P_NA <- which(is.na(P$Data), arr.ind = T)                              # #
      P_flg <- length(P_NA) / 2 != length(P$Data)                            # #
                                                                             # #
      if(P_flg){                                                             # #
        P$Unc[P_NA] <- NA                                                    # #
        X$Data[P_NA] <- NA                                                   # #
        X$Unc[P_NA] <- NA                                                    # #
      }                                                                      # #
      rm(P_NA)                                                               # #
                                                                             # #
      Pbin[,,"Data", n] <- P$Data                                            # #
      Pbin[,,"Unc", n] <- P$Unc                                              # #
      rm(P)                                                                  # #
                                                                             # #
      Xbin[,,"Data", n] <- X$Data                                            # #
      Xbin[,,"Unc", n] <- X$Unc                                              # #
      rm(X)                                                                  # #
                                                                             # #
      IPbin[,,"Data", n] <- Pbin[,,"Data", n] * Ibin[,,"Data", n]            # #
      IPbin[,,"Unc", n] <- unc_mult(Pbin[,,"Data", n], Ibin[,,"Data", n],    # #
                                 Pbin[,,"Unc", n], Ibin[,,"Unc", n])         # #
    }                                                                        # #
    ############################################################################
                                                                               #
    #------------------------------------------------------------------------# #
                                                                               #
    #################### Expanding and Merging Offsets ####################### #
    ############################################################################
    Qmrg <- array(NA, dim = c(dimsObj, length(munc_par)),                    # #
                  dimnames = list(NULL, NULL, munc_par))                     # #
    mrgDim <- dim(Qmrg)                                                      # #
    mrgNam <- dimnames(Qmrg)                                                 # #
                                                                             # #
    Umrg <- array(NA, dim = mrgDim, dimnames = mrgNam)                       # #
    Pmrg <- array(NA, dim = mrgDim, dimnames = mrgNam)                       # #
    IPmrg <- array(NA, dim = mrgDim, dimnames = mrgNam)                      # #
    Xmrg <- array(NA, dim = mrgDim, dimnames = mrgNam)                       # #
                                                                             # #
    if(this != 1 && this != length(Ipaths_labels)){                          # #
      QmrgS <- array(NA, dim = mrgDim, dimnames = mrgNam)                    # #
      UmrgS <- array(NA, dim = mrgDim, dimnames = mrgNam)                    # #
    }                                                                        # #
                                                                             # #
    PmrgS <- array(NA, dim = mrgDim, dimnames = mrgNam)                      # #
    IPmrgS <- array(NA, dim = mrgDim, dimnames = mrgNam)                     # #
    XmrgS <- array(NA, dim = mrgDim, dimnames = mrgNam)                      # #
                                                                             # #
    if(N != 1){                                                              # #
      print(paste0("** Merging offsets of Stokes Q, U for ", t_label, " in ",# #
                   bs, " band **"))                                          # #
                                                                             # #
      Qmrg[,,"Data"] <- apply(temp_Q[,,"Data",], 1:2, median, na.rm = T)     # #
      Qmrg[,,"Unc"] <- unc_median(temp_Q[,,"Data",], temp_Q[,,"Unc",], 1:2)  # #
      Umrg[,,"Data"] <- apply(temp_U[,,"Data",], 1:2, median, na.rm = T)     # #
      Umrg[,,"Unc"] <- unc_median(temp_U[,,"Data",], temp_U[,,"Unc",], 1:2)  # #
      rm(temp_Q, temp_U)                                                     # #
                                                                             # #
      print(paste0("** Expanding and merging offsets of P, IP and X for ",   # #
                   t_label, " in ", bs, " band **"))                         # #
                                                                             # #
      temp_P <- array(NA, dim = c(mrgDim, N), dimnames = binNam)             # #
      temp_IP <- array(NA, dim = c(mrgDim, N), dimnames = binNam)            # #
      temp_X <- array(NA, dim = c(mrgDim, N), dimnames = binNam)             # #
                                                                             # #
      for(n in 1:N){                                                         # #
        if(n == 1){                                                          # #
          temp_P[,,,n] <- bin_expand_map(Pbin[,,,n], binS, mrgDim,           # #
                                         hold_crpix_N[1, 1:2])               # #
          temp_IP[,,,n] <- bin_expand_map(IPbin[,,,n], binS, mrgDim,         # #
                                          hold_crpix_N[1, 1:2])              # #
          temp_X[,,,n] <- bin_expand_map(Xbin[,,,n], binS, mrgDim,           # #
                                         hold_crpix_N[1, 1:2])               # #
        }else{                                                               # #
          rf <- binDimS[n - 1, 1]                                            # #
          cf <- binDimS[n - 1, 2]                                            # #
          temp_P[,,,n] <- bin_expand_map(Pbin[1:rf, 1:cf,,n], binS, mrgDim,  # #
                                         hold_crpix_N[1, 1:2])               # #
          temp_IP[,,,n] <- bin_expand_map(IPbin[1:rf, 1:cf,,n], binS, mrgDim,# #
                                          hold_crpix_N[1, 1:2])              # #
          temp_X[,,,n] <- bin_expand_map(Xbin[1:rf, 1:cf,,n], binS, mrgDim,  # #
                                         hold_crpix_N[1, 1:2])               # #
        }                                                                    # #
      }                                                                      # #
                                                                             # #
      Pmrg[,,"Data"] <- apply(temp_P[,,"Data",], 1:2, median, na.rm = T)     # #
      Pmrg[,,"Unc"] <- unc_median(temp_P[,,"Data",], temp_P[,,"Unc",], 1:2)  # #
      IPmrg[,,"Data"] <- apply(temp_IP[,,"Data",], 1:2, median, na.rm =T)    # #
      IPmrg[,,"Unc"] <- unc_median(temp_IP[,,"Data",], temp_IP[,,"Unc",],    # #
                                   1:2)                                      # #
      Xmrg[,,"Data"] <- apply(temp_X[,,"Data",], 1:2, median, na.rm = T)     # #
      Xmrg[,,"Unc"] <- unc_median(temp_X[,,"Data",], temp_X[,,"Unc",], 1:2)  # #
                                                                             # #
      temp_P <- array(NA, dim = c(mrgDim, N), dimnames = binNam)             # #
      temp_IP <- array(NA, dim = c(mrgDim, N), dimnames = binNam)            # #
      temp_X <- array(NA, dim = c(mrgDim, N), dimnames = binNam)             # #
                                                                             # #
      for(n in 1:N){                                                         # #
        if(n == 1){                                                          # #
          temp_P[,,,n] <- bin_expand_map(Pbin[,,,n], binS, mrgDim,           # #
                                         hold_crpix_N[1, 1:2], dense = F)    # #
          temp_IP[,,,n] <- bin_expand_map(IPbin[,,,n], binS, mrgDim,         # #
                                          hold_crpix_N[1, 1:2], dense = F)   # #
          temp_X[,,,n] <- bin_expand_map(Xbin[,,,n], binS, mrgDim,           # #
                                         hold_crpix_N[1, 1:2], dense = F)    # #
        }else{                                                               # #
          rf <- binDimS[n - 1, 1]                                            # #
          cf <- binDimS[n - 1, 2]                                            # #
          temp_P[,,,n] <- bin_expand_map(Pbin[1:rf, 1:cf,,n], binS, mrgDim,  # #
                                         hold_crpix_N[1, 1:2], dense = F)    # #
          temp_IP[,,,n] <- bin_expand_map(IPbin[1:rf, 1:cf,,n], binS, mrgDim,# #
                                          hold_crpix_N[1, 1:2], dense = F)   # #
          temp_X[,,,n] <- bin_expand_map(Xbin[1:rf, 1:cf,,n], binS, mrgDim,  # #
                                         hold_crpix_N[1, 1:2], dense = F)    # #
        }                                                                    # #
      }                                                                      # #
                                                                             # #
      PmrgS[,,"Data"] <- apply(temp_P[,,"Data",], 1:2, median, na.rm = T)    # #
      PmrgS[,,"Unc"] <- unc_median(temp_P[,,"Data",], temp_P[,,"Unc",], 1:2) # #
      IPmrgS[,,"Data"] <- apply(temp_IP[,,"Data",], 1:2, median,na.rm =T)    # #
      IPmrgS[,,"Unc"] <- unc_median(temp_IP[,,"Data",],temp_IP[,,"Unc",],1:2)# #
      XmrgS[,,"Data"] <- apply(temp_X[,,"Data",], 1:2, median, na.rm = T)    # #
      XmrgS[,,"Unc"] <- unc_median(temp_X[,,"Data",], temp_X[,,"Unc",], 1:2) # #
      rm(temp_P, temp_IP, temp_X)                                            # #
                                                                             # #
                                                                             # #
      if(this != 1 && this != length(Ipaths_labels)){                        # #
        temp_Q <- array(NA, dim = c(mrgDim, N), dimnames = binNam)           # #
        temp_U <- array(NA, dim = c(mrgDim, N), dimnames = binNam)           # #
                                                                             # #
        for(n in 1:N){                                                       # #
          if(n == 1){                                                        # #
            temp_Q[,,,n] <- bin_expand_map(Qbin[,,,n], binS, mrgDim,         # #
                                           hold_crpix_N[1, 1:2], dense = F)  # #
            temp_U[,,,n] <- bin_expand_map(Ubin[,,,n], binS, mrgDim,         # #
                                           hold_crpix_N[1, 1:2], dense = F)  # #
          }else{                                                             # #
            rf <- binDimS[n - 1, 1]                                          # #
            cf <- binDimS[n - 1, 2]                                          # #
            temp_Q[,,,n] <- bin_expand_map(Qbin[1:rf, 1:cf,,n], binS, mrgDim,# #
                                           hold_crpix_N[1, 1:2], dense = F)  # #
            temp_U[,,,n] <- bin_expand_map(Ubin[1:rf, 1:cf,,n], binS, mrgDim,# #
                                           hold_crpix_N[1, 1:2], dense = F)  # #
          }                                                                  # #
        }                                                                    # #
                                                                             # #
        QmrgS[,,"Data"] <- apply(temp_Q[,,"Data",], 1:2, median, na.rm = T)  # #
        QmrgS[,,"Unc"] <- unc_median(temp_Q[,,"Data",],                      # #
                                     temp_Q[,,"Unc",], 1:2)                  # #
        UmrgS[,,"Data"] <- apply(temp_U[,,"Data",], 1:2, median, na.rm = T)  # #
        UmrgS[,,"Unc"] <- unc_median(temp_U[,,"Data",],                      # #
                                     temp_U[,,"Unc",], 1:2)                  # #
      }                                                                      # #
    }else{                                                                   # #
      Qmrg[,,"Data"] <-temp_Q[,,"Data", 1]                                   # #
      Qmrg[,,"Unc"] <- temp_Q[,,"Unc", 1]                                    # #
      Umrg[,,"Data"] <- temp_U[,,"Data", 1]                                  # #
      Umrg[,,"Unc"] <- temp_U[,,"Unc", 1]                                    # #
      rm(temp_Q, temp_U)                                                     # #
                                                                             # #
      # Dense Maps (all pixels filled withing the pixel bin)                 # #
      Pmrg <- bin_expand_map(Pbin[,,,1], binS, mrgDim, hold_crpix_N[1,1:2])  # #
      IPmrg <- bin_expand_map(IPbin[,,,1],binS,mrgDim, hold_crpix_N[1,1:2])  # #
      Xmrg <- bin_expand_map(Xbin[,,,1], binS, mrgDim, hold_crpix_N[1,1:2])  # #
      PmrgS <- bin_expand_map(Pbin[,,,1], binS, mrgDim,                      # #
                              hold_crpix_N[1, 1:2], dense = F)               # #
      IPmrgS <- bin_expand_map(IPbin[,,,1], binS, mrgDim,                    # #
                               hold_crpix_N[1, 1:2], dense = F)              # #
      XmrgS <- bin_expand_map(Xbin[,,,1], binS, mrgDim,                      # #
                              hold_crpix_N[1, 1:2], dense = F)               # #
                                                                             # #
      if(this != 1 && this != length(Ipaths_labels)){                        # #
        # Sparse Maps (one pixel filled in the center of the pixel bin)      # #
        QmrgS <- bin_expand_map(Qbin[,,,1], binS, mrgDim,                    # #
                                hold_crpix_N[1, 1:2], dense = F)             # #
        UmrgS <- bin_expand_map(Ubin[,,,1], binS, mrgDim,                    # #
                                hold_crpix_N[1, 1:2], dense = F)             # #
      }                                                                      # #
    }                                                                        # #
    rm(IPbin)                                                                # #
    ############################################################################
                                                                               #
    #------------------------------------------------------------------------# #
                                                                               #
    ################### Clearing Pol Maps of I NA pixels ##################### #
    ############################################################################
    print(paste0("** Clearing maps of Iflux NA pixels for ", t_label, " in ",# #
                 bs, " band **"))                                            # #
                                                                             # #
    na_mask <- readFITS(na_mask_paths[bs])$imDat                             # #
                                                                             # #
    Qmrg[,,"Data"] <- Qmrg[,,"Data"] * na_mask                               # #
    Qmrg[,,"Unc"] <- Qmrg[,,"Unc"] * na_mask                                 # #
    Umrg[,,"Data"] <- Umrg[,,"Data"] * na_mask                               # #
    Umrg[,,"Unc"] <- Umrg[,,"Unc"] * na_mask                                 # #
    Pmrg[,,"Data"] <- Pmrg[,,"Data"] * na_mask                               # #
    Pmrg[,,"Unc"] <- Pmrg[,,"Unc"] * na_mask                                 # #
    IPmrg[,,"Data"] <- IPmrg[,,"Data"] * na_mask                             # #
    IPmrg[,,"Unc"] <- IPmrg[,,"Unc"] * na_mask                               # #
    Xmrg[,,"Data"] <- Xmrg[,,"Data"] * na_mask                               # #
    Xmrg[,,"Unc"] <- Xmrg[,,"Unc"] * na_mask                                 # #
    rm(na_mask)                                                              # #
    ############################################################################
                                                                               #
    #------------------------------------------------------------------------# #
                                                                               #
    ############ Masking Gaia stars within merged P, IP and X maps ########### #
    ############################################################################
    if(this_I == label_tag){                                                 # #
      print(paste0("** Copying merged Q, U, P, IP and X ", t_label, " maps ",# #
                   "for ", bs, " band and masking MW stars within them **")) # #
                                                                             # #
      Qmrg_skyCstarRM <- Qmrg                                                # #
      Umrg_skyCstarRM <- Umrg                                                # #
      Pmrg_skyCstarRM <- Pmrg                                                # #
      IPmrg_skyCstarRM <- IPmrg                                              # #
      Xmrg_skyCstarRM <- Xmrg                                                # #
                                                                             # #
      QmrgS_skyCstarRM <- QmrgS                                              # #
      UmrgS_skyCstarRM <- UmrgS                                              # #
      PmrgS_skyCstarRM <- PmrgS                                              # #
      IPmrgS_skyCstarRM <- IPmrgS                                            # #
      XmrgS_skyCstarRM <- XmrgS                                              # #
                                                                             # #
      # Loading MW stars mask                                                # #
      star_mask <- readFITS(star_mask_paths[bs])$imDat                       # #
                                                                             # #
      Qmrg_skyCstarRM[,,"Data"] <- Qmrg_skyCstarRM[,,"Data"] * star_mask     # #
      Qmrg_skyCstarRM[,,"Unc"] <- Qmrg_skyCstarRM[,,"Unc"] * star_mask       # #
      Umrg_skyCstarRM[,,"Data"] <- Umrg_skyCstarRM[,,"Data"] * star_mask     # #
      Umrg_skyCstarRM[,,"Unc"] <- Umrg_skyCstarRM[,,"Unc"] * star_mask       # #
      Pmrg_skyCstarRM[,,"Data"] <- Pmrg_skyCstarRM[,,"Data"] * star_mask     # #
      Pmrg_skyCstarRM[,,"Unc"] <- Pmrg_skyCstarRM[,,"Unc"] * star_mask       # #
      IPmrg_skyCstarRM[,,"Data"] <- IPmrg_skyCstarRM[,,"Data"] * star_mask   # #
      IPmrg_skyCstarRM[,,"Unc"] <- IPmrg_skyCstarRM[,,"Unc"] * star_mask     # #
      Xmrg_skyCstarRM[,,"Data"] <- Xmrg_skyCstarRM[,,"Data"] * star_mask     # #
      Xmrg_skyCstarRM[,,"Unc"] <- Xmrg_skyCstarRM[,,"Unc"] * star_mask       # #
                                                                             # #
      QmrgS_skyCstarRM[,,"Data"] <- QmrgS_skyCstarRM[,,"Data"] * star_mask   # #
      QmrgS_skyCstarRM[,,"Unc"] <- QmrgS_skyCstarRM[,,"Unc"] * star_mask     # #
      UmrgS_skyCstarRM[,,"Data"] <- UmrgS_skyCstarRM[,,"Data"] * star_mask   # #
      UmrgS_skyCstarRM[,,"Unc"] <- UmrgS_skyCstarRM[,,"Unc"] * star_mask     # #
      PmrgS_skyCstarRM[,,"Data"] <- PmrgS_skyCstarRM[,,"Data"] * star_mask   # #
      PmrgS_skyCstarRM[,,"Unc"] <- PmrgS_skyCstarRM[,,"Unc"] * star_mask     # #
      IPmrgS_skyCstarRM[,,"Data"] <- IPmrgS_skyCstarRM[,,"Data"] * star_mask # #
      IPmrgS_skyCstarRM[,,"Unc"] <- IPmrgS_skyCstarRM[,,"Unc"] * star_mask   # #
      XmrgS_skyCstarRM[,,"Data"] <- XmrgS_skyCstarRM[,,"Data"] * star_mask   # #
      XmrgS_skyCstarRM[,,"Unc"] <- XmrgS_skyCstarRM[,,"Unc"] * star_mask     # #
                                                                             # #
      rm(star_mask)                                                          # #
    }                                                                        # #
    ############################################################################
                                                                               #
    #------------------------------------------------------------------------# #
                                                                               #
    ################### Calculating Pol Maps Statistics ###################### #
    ############################################################################
    print(paste0("** Calculating statistics of Q, U, P, IP and X for ",      # #
                 t_label, " in ", bs, " band **"))                           # #
                                                                             # #
    # Stats for sky Cut case                                                 # #
    pol_stats["Q", c("min", "max", "median", "mad"), "val"] <-               # #
      c(min(Qmrg[,,"Data"], na.rm = T), max(Qmrg[,,"Data"], na.rm = T),      # #
        median(Qbin[,,"Data",], na.rm = T), mad(Qbin[,,"Data",], na.rm= T))  # #
    pol_stats["Q", "min", c("x", "y")] <-                                    # #
      which(Qmrg[,,"Data"] == pol_stats["Q", "min", "val"], arr.ind= T)[1,]  # #
    pol_stats["Q", "max", c("x", "y")] <-                                    # #
      which(Qmrg[,,"Data"] == pol_stats["Q", "max", "val"], arr.ind= T)[1,]  # #
                                                                             # #
    pol_stats["U", c("min", "max", "median", "mad"), "val"] <-               # #
      c(min(Umrg[,,"Data"], na.rm = T), max(Umrg[,,"Data"], na.rm = T),      # #
        median(Ubin[,,"Data",], na.rm = T), mad(Ubin[,,"Data",], na.rm= T))  # #
    pol_stats["U", "min", c("x", "y")] <-                                    # #
      which(Umrg[,,"Data"] == pol_stats["U", "min", "val"], arr.ind= T)[1,]  # #
    pol_stats["U", "max", c("x", "y")] <-                                    # #
      which(Umrg[,,"Data"] == pol_stats["U", "max", "val"], arr.ind= T)[1,]  # #
                                                                             # #
    pol_stats["P", c("min", "max", "median", "mad"), "val"] <-               # #
      c(min(Pmrg[,,"Data"], na.rm = T), max(Pmrg[,,"Data"], na.rm = T),      # #
        median(Pbin[,,"Data",], na.rm = T), mad(Pbin[,,"Data",], na.rm= T))  # #
    pol_stats["P", "min", c("x", "y")] <-                                    # #
      which(Pmrg[,,"Data"] == pol_stats["P", "min", "val"], arr.ind= T)[1,]  # #
    pol_stats["P", "max", c("x", "y")] <-                                    # #
      which(Pmrg[,,"Data"] == pol_stats["P", "max", "val"], arr.ind= T)[1,]  # #
                                                                             # #
    pol_stats["X", c("min", "max", "median", "mad"), "val"] <-               # #
      c(min(Xmrg[,,"Data"], na.rm = T), max(Xmrg[,,"Data"], na.rm = T),      # #
        median(Xbin[,,"Data",], na.rm = T), mad(Xbin[,,"Data",], na.rm= T))  # #
    pol_stats["X", "min", c("x", "y")] <-                                    # #
      which(Xmrg[,,"Data"] == pol_stats["X", "min", "val"], arr.ind= T)[1,]  # #
    pol_stats["X", "max", c("x", "y")] <-                                    # #
      which(Xmrg[,,"Data"] == pol_stats["X", "max", "val"], arr.ind= T)[1,]  # #
                                                                             # #
    switch(this,                                                             # #
           temp_out <- stats_out_path,                                       # #
           temp_out <- stats_out_path["skyCut"],                             # #
           temp_out <- stats_out_path)                                       # #
                                                                             # #
    write.table(pol_stats, temp_out, sep = ";", dec = ".")                   # #
    rm(Qbin, Ubin, Pbin, Xbin)                                               # #
                                                                             # #
    if(this_I == label_tag){                                                 # #
      pol_stats_skyCstarRM <- pol_stats                                      # #
                                                                             # #
      # Stats for sky and MW stars filtered case                             # #
      pol_stats_skyCstarRM["Q", c("min", "max", "median", "mad"), "val"] <-  # #
        c(min(QmrgS_skyCstarRM[,,"Data"], na.rm = T),                        # #
          max(QmrgS_skyCstarRM[,,"Data"], na.rm = T),                        # #
          median(QmrgS_skyCstarRM[,,"Data"], na.rm = T),                     # #
          mad(QmrgS_skyCstarRM[,,"Data"], na.rm= T))                         # #
      pol_stats_skyCstarRM["Q", "min", c("x", "y")] <-                       # #
        which(QmrgS_skyCstarRM[,,"Data"] ==                                  # #
                pol_stats_skyCstarRM["Q", "min", "val"], arr.ind= T)[1,]     # #
      pol_stats_skyCstarRM["Q", "max", c("x", "y")] <-                       # #
        which(QmrgS_skyCstarRM[,,"Data"] ==                                  # #
                pol_stats_skyCstarRM["Q", "max", "val"], arr.ind= T)[1,]     # #
                                                                             # #
      pol_stats_skyCstarRM["U", c("min", "max", "median", "mad"), "val"] <-  # #
        c(min(UmrgS_skyCstarRM[,,"Data"], na.rm = T),                        # #
          max(UmrgS_skyCstarRM[,,"Data"], na.rm = T),                        # #
          median(UmrgS_skyCstarRM[,,"Data"], na.rm = T),                     # #
          mad(UmrgS_skyCstarRM[,,"Data"], na.rm= T))                         # #
      pol_stats_skyCstarRM["U", "min", c("x", "y")] <-                       # #
        which(UmrgS_skyCstarRM[,,"Data"] ==                                  # #
                pol_stats_skyCstarRM["U", "min", "val"], arr.ind= T)[1,]     # #
      pol_stats_skyCstarRM["U", "max", c("x", "y")] <-                       # #
        which(UmrgS_skyCstarRM[,,"Data"] ==                                  # #
                pol_stats_skyCstarRM["U", "max", "val"], arr.ind= T)[1,]     # #
                                                                             # #
      pol_stats_skyCstarRM["P", c("min", "max", "median", "mad"), "val"] <-  # #
        c(min(PmrgS_skyCstarRM[,,"Data"], na.rm = T),                        # #
          max(PmrgS_skyCstarRM[,,"Data"], na.rm = T),                        # #
          median(PmrgS_skyCstarRM[,,"Data"], na.rm = T),                     # #
          mad(PmrgS_skyCstarRM[,,"Data"], na.rm= T))                         # #
      pol_stats_skyCstarRM["P", "min", c("x", "y")] <-                       # #
        which(PmrgS_skyCstarRM[,,"Data"] ==                                  # #
                pol_stats_skyCstarRM["P", "min", "val"], arr.ind= T)[1,]     # #
      pol_stats_skyCstarRM["P", "max", c("x", "y")] <-                       # #
        which(PmrgS_skyCstarRM[,,"Data"] ==                                  # #
                pol_stats_skyCstarRM["P", "max", "val"], arr.ind= T)[1,]     # #
                                                                             # #
      pol_stats_skyCstarRM["X", c("min", "max", "median", "mad"), "val"] <-  # #
        c(min(XmrgS_skyCstarRM[,,"Data"], na.rm = T),                        # #
          max(XmrgS_skyCstarRM[,,"Data"], na.rm = T),                        # #
          median(XmrgS_skyCstarRM[,,"Data"], na.rm = T),                     # #
          mad(XmrgS_skyCstarRM[,,"Data"], na.rm= T))                         # #
      pol_stats_skyCstarRM["X", "min", c("x", "y")] <-                       # #
        which(XmrgS_skyCstarRM[,,"Data"] ==                                  # #
                pol_stats_skyCstarRM["X", "min", "val"], arr.ind= T)[1,]     # #
      pol_stats_skyCstarRM["X", "max", c("x", "y")] <-                       # #
        which(XmrgS_skyCstarRM[,,"Data"] ==                                  # #
                pol_stats_skyCstarRM["X", "max", "val"], arr.ind= T)[1,]     # #
                                                                             # #
      temp_out <- paste0(stats_out_path["filtH&stars"], ".csv")              # #
      write.table(pol_stats_skyCstarRM, temp_out, sep = ";", dec = ".")      # #
      rm(temp_out)                                                           # #
    }                                                                        # #
    ############################################################################
                                                                               #
    #------------------------------------------------------------------------# #
                                                                               #
    ############# Creating Pol measurement vs uncertainty plots ############## #
    ############################################################################
    if(this != 1 && this != length(Ipaths_labels)){                          # #
      print(paste0("** Creating pdf files with plots of polarization measur",# #
                   "ements vs their uncertainty for ", t_label, " in ", bs,  # #
                   " band **"))                                              # #
                                                                             # #
      if(this_I == label_tag){                                               # #
        MxUnc_out <- paste0(t_out["skyCut", "MxUnc"], basename, "-MxUnc-")   # #
      }                                                                      # #
                                                                             # #
      for(pol_param in 1:5){                                                 # #
        switch(pol_param,                                                    # #
               c(tempPar <- QmrgS * 100, par_lab <- "Q", symb <- "(%)"),     # #
               c(tempPar <- UmrgS * 100, par_lab <- "U", symb <- "(%)"),     # #
               c(tempPar <- PmrgS * 100, par_lab <- "P", symb <- "(%)"),     # #
               c(tempPar <- IPmrgS, par_lab <- "IP", symb <- "(RU)"),        # #
               c(tempPar <- XmrgS, par_lab <- "X", symb <- "(º)")            # #
               )                                                             # #
                                                                             # #
        temp_out <- paste0(MxUnc_out, par_lab, "_", t_label, ".pdf")         # #
                                                                             # #
        temp_par <- as.vector(tempPar[,,1])                                  # #
        temp_par[which(is.infinite(temp_par), arr.ind = T)] <- NA            # #
        temp_NAs <- which(is.na(temp_par), arr.ind = T)                      # #
        temp_par <- temp_par[-temp_NAs]                                      # #
        temp_unc <- as.vector(tempPar[,,2])[-temp_NAs]                       # #
                                                                             # #
        pdf(temp_out, width = 15, height = 15)                               # #
        par(cex.main = 3)                                                    # #
        par(cex.axis = 2)                                                    # #
        par(cex.lab = 2)                                                     # #
                                                                             # #
        if(pol_param == 1 || pol_param == 2){                                # #
          xlabel <- paste0("abs(", par_lab, ") ", symb)                      # #
          mtitl <- paste0("abs(", par_lab, ") ", symb, " vs unc_", par_lab,  # #
                          " ", symb)                                         # #
                                                                             # #
          temp_par <- abs(temp_par)                                          # #
        }else{                                                               # #
          xlabel <- paste0(par_lab, " ", symb)                               # #
          mtitl <- paste0(par_lab, " ", symb, " vs unc_", par_lab, " ", symb)# #
        }                                                                    # #
        ylabel <- paste0("unc_", par_lab, " ", symb)                         # #
                                                                             # #
        y_max <- max(temp_unc)                                               # #
        y_min <- min(temp_unc)                                               # #
        x_max <- max(temp_par)                                               # #
        x_min <- min(temp_par)                                               # #
                                                                             # #
        matplot(temp_par, temp_unc, xlab = xlabel, ylab = ylabel,            # #
                xlim = c(x_min, x_max), ylim = c(y_min, y_max), cex = 2,     # #
                main = mtitl, pch = '+')                                     # #
        dev.off()                                                            # #
        rm(temp_par, temp_unc, temp_NAs)                                     # #
      }                                                                      # #
                                                                             # #
      if(this_I == label_tag){                                               # #
        MxUnc_out <- paste0(t_out["skyCut&starsRM", "MxUnc"], basename,      # #
                            "-MxUnc-")                                       # #
                                                                             # #
        for(pol_param in 1:5){                                               # #
          switch(pol_param,                                                  # #
                 c(tempPar <- QmrgS * 100, par_lab <- "Q", symb <- "(%)"),   # #
                 c(tempPar <- UmrgS * 100, par_lab <- "U", symb <- "(%)"),   # #
                 c(tempPar <- PmrgS * 100, par_lab <- "P", symb <- "(%)"),   # #
                 c(tempPar <- IPmrgS, par_lab <- "IP", symb <- "(RU)"),      # #
                 c(tempPar <- XmrgS, par_lab <- "X", symb <- "(º)")          # #
          )                                                                  # #
                                                                             # #
          temp_out <- paste0(MxUnc_out, par_lab, "_", t_label, ".pdf")       # #
                                                                             # #
          temp_par <- as.vector(tempPar[,,1])                                # #
          temp_par[which(is.infinite(temp_par), arr.ind = T)] <- NA          # #
          temp_NAs <- which(is.na(temp_par), arr.ind = T)                    # #
          temp_par <- temp_par[-temp_NAs]                                    # #
          temp_unc <- as.vector(tempPar[,,2])[-temp_NAs]                     # #
                                                                             # #
          pdf(temp_out, width = 15, height = 15)                             # #
          par(cex.main = 3)                                                  # #
          par(cex.axis = 2)                                                  # #
          par(cex.lab = 2)                                                   # #
                                                                             # #
          if(pol_param == 1 || pol_param == 2){                              # #
            xlabel <- paste0("abs(", par_lab, ") ", symb)                    # #
            mtitl <- paste0("abs(", par_lab, ") ", symb, " vs unc_", par_lab,# #
                            " ", symb)                                       # #
                                                                             # #
            temp_par <- abs(temp_par)                                        # #
          }else{                                                             # #
            xlabel <- paste0(par_lab, " ", symb)                             # #
            mtitl <- paste0(par_lab, " ", symb, " vs unc_", par_lab, " ",    # #
                            symb)                                            # #
          }                                                                  # #
          ylabel <- paste0("unc_", par_lab, " ", symb)                       # #
                                                                             # #
          y_max <- max(temp_unc)                                             # #
          y_min <- min(temp_unc)                                             # #
          x_max <- max(temp_par)                                             # #
          x_min <- min(temp_par)                                             # #
                                                                             # #
          matplot(temp_par, temp_unc, xlab = paste0(par_lab, " ", symb),     # #
                  ylab = paste0("unc_", par_lab, " ", symb),                 # #
                  xlim = c(x_min, x_max), ylim = c(y_min, y_max), cex = 2,   # #
                  main =paste0(par_lab," ",symb," vs unc_",par_lab," ",symb),# #
                  pch = '+')                                                 # #
          dev.off()                                                          # #
          rm(temp_par, temp_unc, temp_NAs)                                   # #
        }                                                                    # #
      }                                                                      # #
    }                                                                        # #
    ############################################################################
                                                                               #
    #------------------------------------------------------------------------# #
                                                                               #
    ###################### Creating Pol Maps FITS files ###################### #
    ############################################################################
    print(paste0("** Creating fits files for maps of Q, U, P, IP and X for ",# #
                 t_label, " in ", bs, " band **"))                           # #
                                                                             # #
    if(N != 1){                                                              # #
      hold_crval <- hold_crval[1:3]                                          # #
      hold_ctype <- hold_ctype[1:3]                                          # #
    }                                                                        # #
    hold_crpix <- c(iflux_pix[,bs], 1)                                       # #
                                                                             # #
    switch(this,                                                             # #
           temp_out <- paste0(t_out, basename),                              # #
           temp_out <- paste0(t_out["skyCut","Dense"], basename),            # #
           temp_out <- paste0(t_out, basename))                              # #
                                                                             # #
    P_head <- Q_head                                                         # #
    IP_head <- Q_head                                                        # #
    X_head <- Q_head                                                         # #
                                                                             # #
    # Modifying headers for fits file creation                               # #
    Q_head <- modVal('FILETYP1', t_type, "File type", Q_head)                # #
    U_head <- modVal('FILETYP1', t_type, "File type", U_head)                # #
    P_head <- modVal('FILETYP1', t_type, "File type", P_head)                # #
    IP_head <- modVal('FILETYP1', t_type, "File type", IP_head)              # #
    X_head <- modVal('FILETYP1', t_type, "File type", X_head)                # #
                                                                             # #
    P_head <- modVal('FILETYP2', 'POL DEG P', "File sub-type", P_head)       # #
    IP_head <- modVal('FILETYP2', 'POL INT IP', "File sub-type", IP_head)    # #
    X_head <- modVal('FILETYP2', 'POL ANG X', "File sub-type", X_head)       # #
                                                                             # #
    # Creating fits files for regular case (dense)                           # #
    writeFITSim(Qmrg, file = paste0(temp_out, "-Merged-Q", fits_out),        # #
                crvaln = hold_crval, crpixn = hold_crpix,                    # #
                ctypen = hold_ctype, header = Q_head)                        # #
    writeFITSim(Umrg, file = paste0(temp_out, "-Merged-U", fits_out),        # #
                crvaln = hold_crval, crpixn = hold_crpix,                    # #
                ctypen = hold_ctype, header = U_head)                        # #
    writeFITSim(Pmrg, file = paste0(temp_out, "-Merged-P", fits_out),        # #
                crvaln = hold_crval, crpixn = hold_crpix,                    # #
                ctypen = hold_ctype, header = P_head)                        # #
    writeFITSim(IPmrg, file=paste0(temp_out, "-Merged-IP", fits_out),        # #
                crvaln = hold_crval, crpixn = hold_crpix,                    # #
                ctypen = hold_ctype, header = IP_head)                       # #
    writeFITSim(Xmrg, file = paste0(temp_out, "-Merged-X", fits_out),        # #
                crvaln = hold_crval, crpixn = hold_crpix,                    # #
                ctypen = hold_ctype, header = X_head)                        # #
    rm(Qmrg, Umrg, IPmrg)                                                    # #
                                                                             # #
                                                                             # #
    if(this_I == label_tag){                                                 # #
      # Creating fits files for case: sky cut, sparse                        # #
      temp_out <- paste(t_out["skyCut","Sparse"], basename)                  # #
                                                                             # #
      writeFITSim(QmrgS, file = paste0(temp_out, "-Merged-Q", fits_out),     # #
                  crvaln = hold_crval, crpixn = hold_crpix,                  # #
                  ctypen = hold_ctype, header = Q_head)                      # #
      writeFITSim(UmrgS, file = paste0(temp_out, "-Merged-U", fits_out),     # #
                  crvaln = hold_crval, crpixn = hold_crpix,                  # #
                  ctypen = hold_ctype, header = U_head)                      # #
      writeFITSim(PmrgS, file = paste0(temp_out, "-Merged-P", fits_out),     # #
                  crvaln = hold_crval, crpixn = hold_crpix,                  # #
                  ctypen = hold_ctype, header = P_head)                      # #
      writeFITSim(IPmrgS, file = paste0(temp_out, "-Merged-IP", fits_out),   # #
                  crvaln = hold_crval, crpixn = hold_crpix,                  # #
                  ctypen = hold_ctype, header = IP_head)                     # #
      writeFITSim(XmrgS, file = paste0(temp_out, "-Merged-X", fits_out),     # #
                  crvaln = hold_crval, crpixn = hold_crpix,                  # #
                  ctypen = hold_ctype, header = X_head)                      # #
                                                                             # #
      # Creating fits files for case (sky cut & MW stars removed, dense)     # #
      temp_out <- paste0(t_out["skyCut&starsRM","Dense"], basename)          # #
                                                                             # #
      writeFITSim(Qmrg_skyCstarRM,                                           # #
                  file = paste0(temp_out, "-Merged-Q", fits_out),            # #
                  crvaln = hold_crval, crpixn = hold_crpix,                  # #
                  ctypen = hold_ctype, header = Q_head)                      # #
      writeFITSim(Umrg_skyCstarRM,                                           # #
                  file = paste0(temp_out, "-Merged-U", fits_out),            # #
                  crvaln = hold_crval, crpixn = hold_crpix,                  # #
                  ctypen = hold_ctype, header = U_head)                      # #
      writeFITSim(Pmrg_skyCstarRM,                                           # #
                  file = paste0(temp_out, "-Merged-P", fits_out),            # #
                  crvaln = hold_crval, crpixn = hold_crpix,                  # #
                  ctypen = hold_ctype, header = P_head)                      # #
      writeFITSim(IPmrg_skyCstarRM,                                          # #
                  file = paste0(temp_out, "-Merged-IP", fits_out),           # #
                  crvaln = hold_crval, crpixn = hold_crpix,                  # #
                  ctypen = hold_ctype, header = IP_head)                     # #
      writeFITSim(Xmrg_skyCstarRM,                                           # #
                  file = paste0(temp_out, "-Merged-X", fits_out),            # #
                  crvaln = hold_crval, crpixn = hold_crpix,                  # #
                  ctypen = hold_ctype, header = X_head)                      # #
      rm(Qmrg_skyCstarRM, Umrg_skyCstarRM, IPmrg_skyCstarRM)                 # #
                                                                             # #
      # Creating fits files for case (sky cut & MW stars removed, sparse)    # #
      temp_out <- paste0(t_out["skyCut&starsRM","Sparse"], basename)         # #
                                                                             # #
      writeFITSim(QmrgS_skyCstarRM,                                          # #
                  file = paste0(temp_out, "-Merged-Q", fits_out),            # #
                  crvaln = hold_crval, crpixn = hold_crpix,                  # #
                  ctypen = hold_ctype, header = Q_head)                      # #
      writeFITSim(UmrgS_skyCstarRM,                                          # #
                  file = paste0(temp_out, "-Merged-U", fits_out),            # #
                  crvaln = hold_crval, crpixn = hold_crpix,                  # #
                  ctypen = hold_ctype, header = U_head)                      # #
      writeFITSim(PmrgS_skyCstarRM,                                          # #
                  file = paste0(temp_out, "-Merged-P", fits_out),            # #
                  crvaln = hold_crval, crpixn = hold_crpix,                  # #
                  ctypen = hold_ctype, header = P_head)                      # #
      writeFITSim(IPmrgS_skyCstarRM,                                         # #
                  file = paste0(temp_out, "-Merged-IP", fits_out),           # #
                  crvaln = hold_crval, crpixn = hold_crpix,                  # #
                  ctypen = hold_ctype, header = IP_head)                     # #
      writeFITSim(XmrgS_skyCstarRM,                                          # #
                  file = paste0(temp_out, "-Merged-X", fits_out),            # #
                  crvaln = hold_crval, crpixn = hold_crpix,                  # #
                  ctypen = hold_ctype, header = X_head)                      # #
    }                                                                        # #
    rm(Q_head, U_head, P_head, IP_head, X_head, hold_crval, hold_ctype,      # #
       temp_out)                                                             # #
                                                                             # #
    print(paste0("** The FITS files for ", t_label," Q, U, P, IP and X in ", # #
                 bs, " have been created **"))                               # #
    ############################################################################
                                                                               #
    #------------------------------------------------------------------------# #
                                                                               #
    ##################### Creating Pol Maps PDF files ######################## #
    ############################################################################
    switch(this,                                                             # #
           c(temp_out_X <- X_col_path, temp_out_P <- P_col_path),            # #
           c(temp_out_X <- X_col_path["skyCut"],                             # #
             temp_out_P <- P_col_path["skyCut"]),                            # #
           c(temp_out_X <- X_col_path, temp_out_P <- P_col_path))            # #
                                                                             # #
    # Adjusting coordinates in array to plotting coordinates                 # #
    dif_x <- hold_crpix[1] - 0:dimsObj[1]                                    # #
    dif_y <- 0:dimsObj[2] - hold_crpix[2]                                    # #
    raseq <- hold_coords[1] + dif_x * scl_x / cosd[b]                        # #
    decseq <- hold_coords[2] + dif_y * scl_y                                 # #
    ra_bot <- round(raseq[round(0.1 * dimsObj[1])], 2)                       # #
    ra_top <- round(raseq[round(0.9 * dimsObj[1])], 2)                       # #
    dec_bot <- round(decseq[round(0.235 * dimsObj[2])], 2)                   # #
    dec_top <- round(decseq[round(0.935 * dimsObj[2])], 2)                   # #
                                                                             # #
    # Setting up axes tick marks and labels                                  # #
    xticks <- seq(0.1, .9, length.out = 5)                                   # #
    yticks <- seq(0.235, .935, length.out = 5)                               # #
    leg_X_ticks <- seq(-90, 90, 30)                                          # #
    leg_P_ticks <- seq(0, max(gP_crp[,,1], na.rm = T), length.out = 7) * 100 # #
    xlabs <- round(seq(ra_bot, ra_top, length.out = 5), 2)                   # #
    ylabs <- round(seq(dec_bot, dec_top, length.out = 5), 2)                 # #
    leg_X_labs <- paste0(leg_X_ticks, "°")                                   # #
    leg_P_labs <- paste0(leg_P_ticks, "%")                                   # #
                                                                             # #
    # Generate a circular color scale with 10 colors                         # #
    num_colors <- 181                                                        # #
    hue_range <- c(0, 360)                                                   # #
    hues <- seq(hue_range[1], hue_range[2], length.out = num_colors)         # #
    colors <- hcl(h = hues, c = 180, l = 65, alpha = 1)                      # #
    rm(num_colors, hue_range, hues)                                          # #
                                                                             # #
    # Polarization angle X PDF for regular case                              # #
    pdf(temp_out_X, width = 30, height = 30)                                 # #
    par(cex.main = 6)                                                        # #
    par(cex.axis = 4)                                                        # #
    par(cex.lab = 5)                                                         # #
                                                                             # #
    # Main Plot                                                              # #
    image.plot(Xmrg[,,1], xlab = ' ', ylab = ' ', xaxt = "n", yaxt = "n",    # #
               col = colors, bigplot = c(.09, .06 + delt, .09, .95),         # #
               smallplot = c(.07 + delt, .08 + delt, .09, .95),              # #
               legend.lab = "Pol. Ang. (deg)", legend.cex = 6,               # #
               legend.line = 15, axis.args = list(at = leg_X_ticks,          # #
                                                  labels = leg_X_labs))      # #
    # Main title                                                             # #
    mtext(paste0("Polarization Angle Map, ", t_label), cex = 4, line = 2)    # #
    # X-axis                                                                 # #
    axis(1, at = xticks, labels = paste0(xlabs, "°"), line = 3, lwd = 0)     # #
    axis(1, at = xticks, labels = rep("", length(xticks)), lwd = 0,          # #
         lwd.ticks = 1)                                                      # #
    mtext("Ra (deg)", side = 1, cex = 6, line = 11)                          # #
    # Y-axis                                                                 # #
    axis(2, at = yticks, labels = paste0(ylabs, "°"))                        # #
    mtext("Dec (deg)", side = 2, cex = 6, line = 8)                          # #
    # Target Center                                                          # #
    points(pp2p[1, 1, b], pp2p[1, 2, b], pch = 16, col = "black", cex = 3)   # #
    # N/E referencial                                                        # #
    arrows(x0 = delt + .1, y0 = 0.05, x1 = delt + .1, y1 = 0.1,              # #
           length = 0.25, col = "red", lwd = 5)                              # #
    arrows(x0 = delt + .1, y0 = 0.05, x1 = delt + .1 - 0.05, y1 = 0.05,      # #
           length = 0.25, col = "red", lwd = 5)                              # #
    text(x = delt + .1, y = 0.105, labels = "N", pos = 3, col = "red",       # #
         cex = 3)                                                            # #
    text(x = delt + .1 - 0.053, y = 0.05, labels = "E", pos = 2,             # #
         col = "red", cex = 3)                                               # #
    dev.off()                                                                # #
    rm(Xmrg)                                                                 # #
                                                                             # #
    # Generate a circular color scale with 10 colors                         # #
    num_colors <- 181                                                        # #
    hue_range <- c(0, 270)                                                   # #
    hues <- seq(hue_range[1], hue_range[2], length.out = num_colors)         # #
    colors <- hcl(h = hues, c = 180, l = 65, alpha = 1)                      # #
    rm(num_colors, hue_range, hues)                                          # #
                                                                             # #
    # Polarization angle P PDF for regular case                              # #
    pdf(temp_out_P, width = 30, height = 30)                                 # #
    par(cex.main = 6)                                                        # #
    par(cex.axis = 4)                                                        # #
    par(cex.lab = 5)                                                         # #
                                                                             # #
    # Main Plot                                                              # #
    image.plot(Pmrg[,,1] * 100, xlab = ' ', ylab = ' ', xaxt = "n", yaxt="n",# #
               col = colors, bigplot = c(.09, .06 + delt, .09, .95),         # #
               smallplot = c(.07 + delt, .08 + delt, .09, .95),              # #
               legend.lab = "Pol. Degree (%)", legend.cex = 6,               # #
               legend.line = 15, axis.args = list(at = leg_P_ticks,          # #
                                                  labels = leg_P_labs))      # #
    # Main title                                                             # #
    mtext(paste0("Polarization Degree Map, ", t_label), cex = 4, line = 2)   # #
    # X-axis                                                                 # #
    axis(1, at = xticks, labels = paste0(xlabs, "°"), line = 3, lwd = 0)     # #
    axis(1, at = xticks, labels = rep("", length(xticks)), lwd = 0,          # #
         lwd.ticks = 1)                                                      # #
    mtext("Ra (deg)", side = 1, cex = 6, line = 11)                          # #
    # Y-axis                                                                 # #
    axis(2, at = yticks, labels = paste0(ylabs, "°"))                        # #
    mtext("Dec (deg)", side = 2, cex = 6, line = 8)                          # #
    # Target Center                                                          # #
    points(pp2p[1, 1, b], pp2p[1, 2, b], pch = 16, col = "black", cex = 3)   # #
    # N/E referencial                                                        # #
    arrows(x0 = delt + .1, y0 = 0.05, x1 = delt + .1, y1 = 0.1,              # #
           length = 0.25, col = "red", lwd = 5)                              # #
    arrows(x0 = delt + .1, y0 = 0.05, x1 = delt + .1 - 0.05, y1 = 0.05,      # #
           length = 0.25, col = "red", lwd = 5)                              # #
    text(x = delt + .1, y = 0.105, labels = "N", pos = 3, col = "red",       # #
         cex = 3)                                                            # #
    text(x = delt + .1 - 0.053, y = 0.05, labels = "E", pos = 2,             # #
         col = "red", cex = 3)                                               # #
    dev.off()                                                                # #
    rm(Pmrg)                                                                 # #
                                                                             # #
    if(this_I == label_tag){                                                 # #
      temp_out_X <- X_col_path["skyCut&starsRM"]                             # #
      temp_out_P <- P_col_path["skyCut&starsRM"]                             # #
                                                                             # #
      # Generate a circular color scale with 10 colors                       # #
      num_colors <- 181                                                      # #
      hue_range <- c(0, 360)                                                 # #
      hues <- seq(hue_range[1], hue_range[2], length.out = num_colors)       # #
      colors <- hcl(h = hues, c = 180, l = 65, alpha = 1)                    # #
      rm(num_colors, hue_range, hues)                                        # #
                                                                             # #
      # Polarization angle X PDF for sky cut and MW stars removed case       # #
      pdf(temp_out_X, width = 30, height = 30)                               # #
      par(cex.main = 6)                                                      # #
      par(cex.axis = 4)                                                      # #
      par(cex.lab = 5)                                                       # #
                                                                             # #
      # Main Plot                                                            # #
      image.plot(Xmrg_skyCstarRM[,,1], xlab = ' ', ylab = ' ', xaxt = "n",   # #
                 yaxt = "n", legend.lab = "Pol. Ang. (deg)", legend.cex = 6, # #
                 bigplot = c(.09, .06 + delt, .09, .95), legend.line = 15,   # #
                 smallplot = c(.07 + delt, .08 + delt, .09, .95),            # #
                 col = colors, axis.args = list(at = leg_X_ticks,            # #
                                                labels = leg_X_labs))        # #
      # Main title                                                           # #
      mtext(paste0("Polarization Angle Map, ", t_label), cex = 4, line = 2)  # #
      # X-axis                                                               # #
      axis(1, at = xticks, labels = paste0(xlabs, "°"), line = 3, lwd = 0)   # #
      axis(1, at = xticks, labels = rep("", length(xticks)), lwd = 0,        # #
           lwd.ticks = 1)                                                    # #
      mtext("Ra (deg)", side = 1, cex = 6, line = 11)                        # #
      # Y-axis                                                               # #
      axis(2, at = yticks, labels = paste0(ylabs, "°"))                      # #
      mtext("Dec (deg)", side = 2, cex = 6, line = 8)                        # #
      # Target Center                                                        # #
      points(pp2p[1, 1, b], pp2p[1, 2, b], pch = 16, col = "black", cex = 3) # #
      # N/E referencial                                                      # #
      arrows(x0 = delt + .1, y0 = 0.05, x1 = delt + .1, y1 = 0.1,            # #
             length = 0.25, col = "red", lwd = 5)                            # #
      arrows(x0 = delt + .1, y0 = 0.05, x1 = delt + .1 - 0.05, y1 = 0.05,    # #
             length = 0.25, col = "red", lwd = 5)                            # #
      text(x = delt + .1, y = 0.105, labels = "N", pos = 3, col = "red",     # #
           cex = 3)                                                          # #
      text(x = delt + .1 - 0.053, y = 0.05, labels = "E", pos = 2,           # #
           col = "red", cex = 3)                                             # #
      dev.off()                                                              # #
      rm(Xmrg_skyCstarRM, temp_out_X)                                        # #
                                                                             # #
      # Generate a circular color scale with 10 colors                       # #
      num_colors <- 181                                                      # #
      hue_range <- c(0, 270)                                                 # #
      hues <- seq(hue_range[1], hue_range[2], length.out = num_colors)       # #
      colors <- hcl(h = hues, c = 180, l = 65, alpha = 1)                    # #
      rm(num_colors, hue_range, hues)                                        # #
                                                                             # #
      # Polarization angle P PDF for regular case                            # #
      pdf(temp_out_P, width = 30, height = 30)                               # #
      par(cex.main = 6)                                                      # #
      par(cex.axis = 4)                                                      # #
      par(cex.lab = 5)                                                       # #
                                                                             # #
      # Main Plot                                                            # #
      image.plot(Pmrg_skyCstarRM[,,1] * 100, xlab =' ', ylab =' ', xaxt ="n",# #
                 yaxt="n", col =colors, bigplot = c(.09, .06 + delt,.09,.95),# #
                 smallplot = c(.07 + delt, .08 + delt, .09, .95),            # #
                 legend.lab = "Pol. Degree (%)", legend.cex = 6,             # #
                 legend.line = 15, axis.args = list(at = leg_P_ticks,        # #
                                                    labels = leg_P_labs))    # #
      # Main title                                                           # #
      mtext(paste0("Polarization Degree Map, ", t_label), cex = 4, line = 2) # #
      # X-axis                                                               # #
      axis(1, at = xticks, labels = paste0(xlabs, "°"), line = 3, lwd = 0)   # #
      axis(1, at = xticks, labels = rep("", length(xticks)), lwd = 0,        # #
           lwd.ticks = 1)                                                    # #
      mtext("Ra (deg)", side = 1, cex = 6, line = 11)                        # #
      # Y-axis                                                               # #
      axis(2, at = yticks, labels = paste0(ylabs, "°"))                      # #
      mtext("Dec (deg)", side = 2, cex = 6, line = 8)                        # #
      # Target Center                                                        # #
      points(pp2p[1, 1, b], pp2p[1, 2, b], pch = 16, col = "black", cex = 3) # #
      # N/E referencial                                                      # #
      arrows(x0 = delt + .1, y0 = 0.05, x1 = delt + .1, y1 = 0.1,            # #
             length = 0.25, col = "red", lwd = 5)                            # #
      arrows(x0 = delt + .1, y0 = 0.05, x1 = delt + .1 - 0.05, y1 = 0.05,    # #
             length = 0.25, col = "red", lwd = 5)                            # #
      text(x = delt + .1, y = 0.105, labels = "N", pos = 3, col = "red",     # #
           cex = 3)                                                          # #
      text(x = delt + .1 - 0.053, y = 0.05, labels = "E", pos = 2,           # #
           col = "red", cex = 3)                                             # #
      dev.off()                                                              # #
      rm(Pmrg_skyCstarRM, temp_out_P)                                        # #
    }                                                                        # #
                                                                             # #
    if(N != 1){                                                              # #
      Ibkg <- readFITS(Merged_Ipaths[this_I, b])$imDat[,,1]                  # #
    }else{                                                                   # #
      Ibkg <- readFITS(Stacked_Ipaths[this_I, b])$imDat[,,1]                 # #
    }                                                                        # #
                                                                             # #
    gal_mask <- readFITS(gal_mask_paths[bs])$imDat                           # #
                                                                             # #
    switch(this,                                                             # #
           temp_out <- pdf_fold,                                             # #
           temp_out <- pdf_fold["skyCut"],                                   # #
           temp_out <- pdf_fold)                                             # #
                                                                             # #
    # Polarization Arrow Plots for sky Cut case                              # #
    create_arrow_plot_pdf(mag = PmrgS[,,"Data"], u_mag = PmrgS[,,"Unc"],     # #
                          ang = XmrgS[,,"Data"], u_ang = XmrgS[,,"Unc"],     # #
                          bkg = Ibkg, outpath = temp_out, binsize = binS,    # #
                          outname = paste0(basename, "-Merged_Arrow_PX_",    # #
                                           t_label),                         # #
                          mtitle = paste0("Polarization Map, ", t_label),    # #
                          points = pp2p[,,b], mask = gal_mask,               # #
                          x_ticks = xticks, y_ticks = yticks,                # #
                          x_labs = xlabs, y_labs = ylabs,                    # #
                          ref_pix = hold_crpix)                              # #
    rm(PmrgS)                                                                # #
                                                                             # #
    create_arrow_plot_pdf(mag = IPmrgS[,,"Data"],u_mag = IPmrgS[,,"Unc"],    # #
                          ang = XmrgS[,,"Data"], u_ang = XmrgS[,,"Unc"],     # #
                          bkg = Ibkg, outpath = temp_out, binsize = binS,    # #
                          outname = paste0(basename, "-Merged_Arrow_IPX_",   # #
                                           t_label),                         # #
                          mtitle = paste0("Polarization Flux Map, ",         # #
                                          t_label),                          # #
                          points = pp2p[,,b], mask = gal_mask,               # #
                          x_ticks = xticks, y_ticks = yticks,                # #
                          x_labs = xlabs, y_labs = ylabs, stat_flag = FALSE, # #
                          ref_pix = hold_crpix)                              # #
    rm(IPmrgS, XmrgS)                                                        # #
                                                                             # #
    if(this_I == label_tag){                                                 # #
      temp_out <- pdf_fold["skyCut&starsRM"]                                 # #
                                                                             # #
      # Polarization Arrow Plots for sky cut and MW stars removed case       # #
      create_arrow_plot_pdf(mag = PmrgS_skyCstarRM[,,"Data"],                # #
                            u_mag = PmrgS_skyCstarRM[,,"Unc"],               # #
                            ang = XmrgS_skyCstarRM[,,"Data"],                # #
                            u_ang = XmrgS_skyCstarRM[,,"Unc"],               # #
                            bkg = Ibkg, outpath = temp_out, binsize = binS,  # #
                            outname = paste0(basename, "-Merged_Arrow_PX_",  # #
                                             t_label),                       # #
                            mtitle = paste0("Polarization Map, ", t_label),  # #
                            points = pp2p[,,b], mask = gal_mask,             # #
                            x_ticks = xticks, y_ticks = yticks,              # #
                            x_labs = xlabs, y_labs = ylabs,                  # #
                            ref_pix = hold_crpix)                            # #
      rm(PmrgS_skyCstarRM)                                                   # #
                                                                             # #
      create_arrow_plot_pdf(mag = IPmrgS_skyCstarRM[,,"Data"],               # #
                            u_mag = IPmrgS_skyCstarRM[,,"Unc"],              # #
                            ang = XmrgS_skyCstarRM[,,"Data"],                # #
                            u_ang = XmrgS_skyCstarRM[,,"Unc"],               # #
                            bkg = Ibkg, outpath = temp_out, binsize = binS,  # #
                            outname = paste0(basename, "-Merged_Arrow_IPX_", # #
                                             t_label),                       # #
                            mtitle = paste0("Polarization Flux Map, ",       # #
                                            t_label),                        # #
                            points = pp2p[,,b], mask = gal_mask,             # #
                            x_ticks = xticks, y_ticks = yticks,              # #
                            x_labs = xlabs, y_labs = ylabs, stat_flag =FALSE,# #
                            ref_pix = hold_crpix)                            # #
      rm(IPmrgS_skyCstarRM, XmrgS_skyCstarRM)                                # #
    }                                                                        # #
    rm(Ibkg, gal_mask, temp_out)                                             # #
                                                                             # #
    print(paste0("** The pdfs for the ", t_label, " polarization maps in ",  # #
                 bs, " have been created **"))                               # #
    ############################################################################
                                                                               #
    print(paste0("** Finished processing ", t_label, " polarization maps in ", #
                 bs, " band **"))                                              #
  }                                                                            #
  print(paste0("## Finished processing all ", t_label, " polarization maps.")) #
}                                                                              #
rm(Merged_Ipaths, Stacked_Ipaths, t_label, band_str, b, b_ind, bs, fileList)   #
################################################################################

print(paste0("Finished processing OFFSETS with Binning of ", binS, "x", binS,
             " pix²."))