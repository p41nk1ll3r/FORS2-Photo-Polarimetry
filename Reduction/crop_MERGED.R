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
library(EBImage)

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

print(paste0("Cropping of target ", og_obj, " with bin size ", binS,
             " has started."))
  
home_folder <- "/media/joaomfras/GalaxyPol/Pol-Gal/"
libs_folder <- paste0(home_folder, "Polarimetric-Imaging-Reduction-Scripts-(FO",
                      "RS2)/Commit/")
lib_head_path <- paste0(libs_folder, "fetch_fits_header_info.R")
lib_pro_path <- paste0(libs_folder, "process_lib.R")
source(lib_head_path)
source(lib_pro_path)
rm(lib_head_path, lib_pro_path, libs_folder)

REFX <- "CRPIX1"
REFY <- "CRPIX2"
VALX <- "CRVAL1"
VALY <- "CRVAL2"
TYPEX <- "CTYPE1"
TYPEY <- "CTYPE2"
REFRA <- "RA     "
REFDEC <- "DEC    "
TARGET <- "OBJECT  ="
PIXSCALE <- "HIERARCH ESO INS PIXSCALE"
EXPT <- "EXPTIME"
BINX <- "HIERARCH ESO DET WIN1 BINX"
BINY <- "HIERARCH ESO DET WIN1 BINY"
BINS <- "BINNING ="


valz <- 1
valz2 <- 1
typex <- "  'PIXEL     '              / Coordinate system of x-axis "
typey <- "  'PIXEL     '              / Coordinate system of y-axis "
typez <- "  'DATA & UNC'              / Coordinate system of z-axis "
typez2 <- "  'I, P, IP & X'            / Coordinate system of z2-axis "


################################################################################
################################# User Prompts #################################
################################################################################
astro_cfg_path <- "/etc/astrometry.cfg"                                        #
bands <- c("u_HIGH", "b_HIGH", "v_HIGH", "R_SPECIAL", "I_BESS")                #
blen <- length(bands)                                                          #
verse <- paste0("bin", binS, "/")                                              #
main_folder <- paste0(home_folder, og_obj, "/Processed_OBJ/")                  #
in_folder <- paste0(main_folder, "Merged_Offsets_SMAx1.2/")                    #
                                                                               #
acq_folder <- paste0(home_folder, og_obj, "/Processed_ACQ/")                   #
acq_mrg_folder <- paste0(acq_folder, "Merged_ACQ/")                            #
acq_astro_folder <- paste0(acq_folder, "Astrometry_ACQ/")                      #
if(!dir.exists(acq_astro_folder)){                                             #
  mkdir(acq_astro_folder)                                                      #
}                                                                              #
                                                                               #
g_folds <- paste0(in_folder, "Stokes/", verse, "All_C/Sky_Cut-&-Stars_RM/",    #
                  "Dense/")                                                    #
g_folds_uncut <- paste0(in_folder, "Stokes/", verse, "All_C/Sky_Cut/Dense/")   #
s_folds <- paste0(in_folder, "Stokes/", verse, "Sky/")                         #
io_folds <- paste0(in_folder, "I_flux/Obs/")                                   #
is_folds <- paste0(in_folder, "I_flux/Sky/")                                   #
                                                                               #
# List and select acquisition file                                             #
acqList <- list.files(acq_mrg_folder, full.names = T)                          #
acqNames <- list.files(acq_mrg_folder, full.names = F)                         #
                                                                               #
is_acq <- grep("Acq-no-CRs.fits", acqList)                                     #
acqLen <- length(is_acq)                                                       #
                                                                               #
if(acqLen != 0){                                                               #
  acqList <- acqList[is_acq]                                                   #
  acqNames <- acqNames[is_acq]                                                 #
}else{                                                                         #
  print("ERROR: No Acquisition files found! Please run 'process_ACQ.R'.")      #
  print("Returning NULL.")                                                     #
  return(NULL)                                                                 #
}                                                                              #
rm(is_acq)                                                                     #
                                                                               #
if(acqLen != 1){                                                               #
  # Check exposure time of each file                                           #
  exptimes <- get_fits_header_list_num(acqList, EXPT)                          #
  best_acq <- which(exptimes == max(exptimes), arr.ind = T)                    #
  acq_file <- acqNames[best_acq]                                               #
  acq_path <- acqList[best_acq]                                                #
  rm(best_acq, exptimes)                                                       #
}else{                                                                         #
  acq_file <- acqNames[1]                                                      #
  acq_path <- acqList[1]                                                       #
}                                                                              #
rm(acqNames, acqList)                                                          #
                                                                               #
acq_expt <- get_fits_header_num(readFITS(acq_path)$header, EXPT)               #
acq_name <- paste0(strsplit(acq_file, "_#")[[1]][1], "_Astro_acq")             #
temp_acq_path <- paste0(acq_folder, acq_name, ".fits")                         #
astro_acq_path <- paste0(acq_astro_folder, acq_name, ".new")                   #
                                                                               #
astrometry_folder <- paste0(home_folder, og_obj, "/Astrometry/Results/")       #
astro_files <- array("", dim = blen, dimnames = list(bands))                   #
                                                                               #
mult <- 1.2                                                                    #####
                                                                               #
main_out <- paste0(main_folder, "Cropped/")                                    #
out_folder <- paste0(main_out, verse)                                          #
if(!dir.exists(out_folder)){                                                   #
  mkdir(out_folder)                                                            #
}                                                                              #
out_filt_folder <- paste0(out_folder, "Filtered/")                             #
out_unfilt_folder <- paste0(out_folder, "Unfiltered/")                         #
if(!dir.exists(out_filt_folder)){                                              #
  mkdir(out_filt_folder)                                                       #
}                                                                              #
if(!dir.exists(out_unfilt_folder)){                                            #
  mkdir(out_unfilt_folder)                                                     #
}                                                                              #
                                                                               #
#Wavelength grid                                                               #
wl_arr <- array(0, dim = c(blen, 2), dimnames = list(bands,                    #
                                                     c("filt_c", "filt_fwhm")))#
wl_arr["u_HIGH",] <- c(361, 50.5)                                              #
wl_arr["b_HIGH",] <- c(437, 102/2)                                             #
wl_arr["v_HIGH",] <- c(555, 123.2/2)                                           #
wl_arr["R_SPECIAL",] <- c(655, 165/2)                                          #
wl_arr["I_BESS",] <- c(768, 138/2)                                             #
wl_seq_nm <- seq(wl_arr[1,1] - 100, wl_arr[blen,1] + 100, 10)                  #
wl_seq <- wl_seq_nm / 1000                                                     #
                                                                               #
#Lists of Files                                                                #
galPs <- array(NA, dim = blen, dimnames = list(bands))                         #
galIPs <-array(NA, dim = blen, dimnames = list(bands))                         #
galXs <- array(NA, dim = blen, dimnames = list(bands))                         #
galQs <- array(NA, dim = blen, dimnames = list(bands))                         #
galUs <- array(NA, dim = blen, dimnames = list(bands))                         #
galPs_uncut <- array(NA, dim = blen, dimnames = list(bands))                   #
galIPs_uncut <-array(NA, dim = blen, dimnames = list(bands))                   #
galXs_uncut <- array(NA, dim = blen, dimnames = list(bands))                   #
galQs_uncut <- array(NA, dim = blen, dimnames = list(bands))                   #
galUs_uncut <- array(NA, dim = blen, dimnames = list(bands))                   #
obsIs <- array(NA, dim = blen, dimnames = list(bands))                         #
skyPs <- array(NA, dim = blen, dimnames = list(bands))                         #
skyXs <- array(NA, dim = blen, dimnames = list(bands))                         #
skyIPs <- array(NA, dim = blen, dimnames = list(bands))                        #
skyIs <- array(NA, dim = blen, dimnames = list(bands))                         #
                                                                               #
Gstats <- array(NA, dim = c(8, blen),                                          #
                dimnames = list(c("<P>(%)", "unc_<P>(%)", "<P/unc_P>", "<IP>", #
                                  "unc_<IP>", "<X>(deg)", "unc_<X>(deg)",      #
                                  "CF"), bands))                               #
G_uncut_stats <- array(NA, dim = c(8, blen),                                   #
                       dimnames = list(c("<P>(%)", "unc_<P>(%)", "<P/unc_P>",  #
                                         "<IP>", "unc_<IP>", "<X>(deg)",       #
                                         "unc_<X>(deg)", "CF"), bands))        #
Sstats <- array(NA, dim = c(6, blen),                                          #
                dimnames = list(c("<P>(%)", "unc_<P>(%)", "<IP>", "unc_<IP>",  #
                                  "<X>(deg)", "unc_<X>(deg)"), bands))         #
                                                                               #
print("Listing files paths...")                                                #
io_files <- list.files(io_folds, full.names = T)                               #
is_files <- list.files(is_folds, full.names = T)                               #
s_files <- list.files(s_folds, full.names = T)                                 #
g_files <- list.files(g_folds, full.names = T)                                 #
g_files_uncut <- list.files(g_folds_uncut, full.names = T)                     #
                                                                               #
astroList <- list.files(astrometry_folder, full.names = T)                     #
astroList <- astroList[grep(".new", astroList)]                                #
                                                                               #
b_files_exist <- array(FALSE, dim = blen, dimname = list(bands))               #
                                                                               #
# Flagging which bands have files available                                    #
# Only checks for one type of file because it assumes all files result from    #
# same procedures                                                              #
for(b in bands){                                                               #
  test_list <- grep(b, io_files, value = T )                                   #
                                                                               #
  if(length(test_list) != 0){                                                  #
    b_files_exist[b] <- TRUE                                                   #
  }                                                                            #
}                                                                              #
rm(test_list)                                                                  #
                                                                               #
# Listing available band files                                                 #
for(b in bands){                                                               #
                                                                               #
  if(!b_files_exist[b]){                                                       #
    next                                                                       #
  }                                                                            #
                                                                               #
  obsIs[b] <- grep("Merged_I-Flux_obs", grep(b, io_files, value =T), value = T)#
  skyIs[b] <- grep("Merged_I-Flux_sky", grep(b, is_files, value =T), value = T)#
  skyPs[b] <- grep("Merged-P", grep(b, s_files, value = T), value = T)         #
  skyXs[b] <- grep("Merged-X", grep(b, s_files, value = T), value = T)         #
  skyIPs[b] <- grep("Merged-IP", grep(b, s_files, value = T), value = T)       #
  galPs[b] <- grep("Merged-P", grep(b, g_files, value = T), value = T)         #
  galPs_uncut[b] <- grep("Merged-P", grep(b, g_files_uncut, value=T), value =T)#
  galIPs[b] <- grep("Merged-IP", grep(b, g_files, value = T), value = T)       #
  galIPs_uncut[b] <- grep("Merged-IP", grep(b, g_files_uncut, value=T),value=T)#
  galXs[b] <- grep("Merged-X", grep(b, g_files, value = T), value = T)         #
  galXs_uncut[b] <- grep("Merged-X", grep(b, g_files_uncut, value =T), value=T)#
  galQs[b] <- grep("Merged-Q", grep(b, g_files, value = T), value = T)         #
  galQs_uncut[b] <- grep("Merged-Q", grep(b, g_files_uncut, value=T), value =T)#
  galUs[b] <- grep("Merged-U", grep(b, g_files, value = T), value = T)         #
  galUs_uncut[b] <- grep("Merged-U", grep(b, g_files_uncut, value=T), value =T)#
                                                                               #
  is_b <- grep(b, astroList)                                                   #
                                                                               #
  if(length(is_b) != 0){                                                       #
    astro_files[b] <- astroList[is_b]                                          #
  }                                                                            #
}                                                                              #
                                                                               #
# Checking if all needed bands have an astrometry file associated to them      #
# If not, associate the astrometry file of the closest band                    #
wcs_adapt_flag <- rep(F, blen)                                                 #
bands_adapt <- rep(NA, blen)                                                   #
for(b in 1:blen){                                                              #
  bs <- bands[b]                                                               #
                                                                               #
  if(!b_files_exist[b]){                                                       #
    next                                                                       #
  }                                                                            #
                                                                               #
  max_ts <- blen - 1                                                           #
  t_bs <- NULL                                                                 #
  ts <- 1                                                                      #
                                                                               #
  while(length(t_bs) != max_ts){                                               #
    if(b == 1){                                                                #
      poss_ts <- b + ts                                                        #
    }                                                                          #
    if(b == blen){                                                             #
      poss_ts <- b - ts                                                        #
    }                                                                          #
    if(b > 1 & b < blen){                                                      #
      poss_ts <- c(b + ts, b - ts)                                             #
    }                                                                          #
                                                                               #
    neg_ts <- which(poss_ts <= 0, arr.ind = T)                                 #
                                                                               #
    if(length(neg_ts) != 0){                                                   #
      poss_ts <- poss_ts[-neg_ts]                                              #
    }                                                                          #
                                                                               #
    out_ts <- which(poss_ts > blen, arr.ind = T)                               #
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
  while(astro_files[b] == "" & ts < blen){                                     #
    wcs_adapt_flag[b] <- T                                                     #
    t_f <- t_bs[ts]                                                            #
    astro_files[b] <- astro_files[t_f]                                         #
    bands_adapt[b] <- bands[t_f]                                               #
    ts <- ts + 1                                                               #
  }                                                                            #
                                                                               #
  if(astro_files[b] == "" & ts == blen){                                       #
    print(paste0("ERROR: No Astrometry file available for this field. Correct",#
                 " the astrometry.net input file and try again. Returning NUL",#
                 "L."))                                                        #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
}                                                                              #
                                                                               #
temp_fits <- readFITS(acq_path)                                                #
dimsObj <- dim(temp_fits$imDat)[1:2]                                           #
                                                                               #
mask <- array(0, dim = dimsObj)                                                #
acq_mask <- array(0, dim = dimsObj)                                            #
                                                                               #
acq_head <- temp_fits$header                                                   #
rm(temp_fits)                                                                  #
                                                                               #
acq_name <- get_fits_header_str(acq_head, TARGET)                              #
basename <- paste0("FORS2_", acq_name, "_")                                    #
                                                                               #
# Get reference coords for acq                                                 #
binx <- get_fits_header_num(acq_head, BINX)                                    #
biny <- get_fits_header_num(acq_head, BINY)                                    #
pxscl <- get_fits_header_num(acq_head, PIXSCALE)                               #
ascl <- pxscl * binx                                                           #
scl_x <- pxscl * binx / 3600                                                   #
scl_y <- pxscl * biny / 3600                                                   #
acq_val <- c(get_fits_header_num(acq_head, VALX),                              #
             get_fits_header_num(acq_head, VALY))                              #
                                                                               #
# Get Coordinates object within acquisition                                    #
tgt_coords <- get_NED_coords(acq_name, 6)                                      #
temp_a <- get_NED_median_diam(acq_name)                                        #
tgt_a <- round((temp_a + sqrt(temp_a)) / (2 * ascl))                           #
                                                                               #
# Performing Astrometry on Acq to align it in PDF files                        #
if(!file.exists(astro_acq_path)){                                              #
  acqIso <- readFITS(acq_path)$imDat[,,1]                                      #
                                                                               #
  writeFITSim(acqIso, temp_acq_path)                                           #
                                                                               #
  astrometry_cmd <- paste0("solve-field --config ", astro_cfg_path, " --ra ",  #
                           tgt_coords[1], " --dec ", tgt_coords[2], " --radiu",#
                           "s .25 --dir ", acq_astro_folder," ", temp_acq_path)#
  system(astrometry_cmd)                                                       #
}                                                                              #
                                                                               #
acq_xy <- convert_wcs_to_pixel(array(tgt_coords, dim = c(1,2)), astro_acq_path)#
                                                                               #
b_over_a <- get_NED_largest_axis_ratio(acq_name)                               #
tgt_pa <- (get_NED_newest_posang(acq_name) - 90) / 180 * pi                    #
tgt_b <- b_over_a * tgt_a                                                      #
                                                                               #
if(is.na(b_over_a) | is.na(tgt_pa)){                                           #
  # Loading Acq image and target info to python                                #
  load_cmd <- paste0("acqF = fits.open('", acq_path, "')")                     #
  py_run_string(load_cmd)                                                      #
  rm(load_cmd)                                                                 #
                                                                               #
  py_run_string("img_dat = acqF[0].data[0]")                                   #
  py_run_string("acqF.close()")                                                #
  py_run_string("X = r.acq_xy[0][0]")                                          #
  py_run_string("Y = r.acq_xy[0][1]")                                          #
  py_run_string("SMA = r.tgt_a")                                               #
                                                                               #
  # Python commands to get shape of target mask                                #
  get_cmd <- paste0("iso = sexD.get_closest_iso_for_target(img_dat, X, Y, SMA",#
                    ", 0, 90)")                                                #
  py_run_string(get_cmd)                                                       #
  rm(get_cmd)                                                                  #
}                                                                              #
                                                                               #
if(is.na(b_over_a)){                                                           #
  tgt_eps <- py$iso$eps                                                        #
  tgt_b <- tgt_a * sqrt(1 - tgt_eps)                                           #
  rm(tgt_eps)                                                                  #
}                                                                              #
if(is.na(tgt_pa)){                                                             #
  tgt_pa <- py$iso$pa                                                          #
}                                                                              #
                                                                               #
tgt_a <- tgt_a * mult                                                          #
tgt_b <- tgt_b * mult                                                          #
                                                                               #
########### Creating cropped acquisition file, will be useful later ########## #
################################################################################
acq_data <- readFITS(acq_path)$imDat[,,1]                                    # #
                                                                             # #
for(xa in 1:dimsObj[1]){                                                     # #
  for(ya in 1:dimsObj[2]){                                                   # #
    rx <- ((xa - acq_xy[1]) * cos(tgt_pa) + (ya - acq_xy[2]) *               # #
             sin(tgt_pa))^2 / tgt_a^2                                        # #
    ry <- ((xa - acq_xy[1]) * sin(tgt_pa) - (ya - acq_xy[2]) *               # #
             cos(tgt_pa))^2 / tgt_b^2                                        # #
                                                                             # #
    if(rx + ry <= 1){                                                        # #
      acq_mask[xa, ya] <- 1                                                  # #
    }                                                                        # #
  }                                                                          # #
}                                                                            # #
rm(xa, ya, rx, ry)                                                           # #
                                                                             # #
round_acq_xy <- round(c(acq_xy[1], acq_xy[2]))                               # #
                                                                             # #
max_x <- max(which(acq_mask == 1, arr.ind = T)[,1])                          # #
min_x <- min(which(acq_mask == 1, arr.ind = T)[,1])                          # #
max_y <- max(which(acq_mask == 1, arr.ind = T)[,2])                          # #
min_y <- min(which(acq_mask == 1, arr.ind = T)[,2])                          # #
                                                                             # #
delX <- max(abs(c(max_x, min_x) - round_acq_xy[1]))                          # #
delY <- max(abs(c(max_y, min_y) - round_acq_xy[2]))                          # #
round_del_xy <- round(c(delX, delY))                                         # #
rm(delX, delY)                                                               # #
                                                                             # #
bkgI_crp_h <- crop_array(acq_data * acq_mask, round_acq_xy, round_del_xy)    # #
dims_bkg <- dim(bkgI_crp_h)                                                  # #
                                                                             # #
# This data will be used as background for arrow plots                       # #
# Check whats maximum flux of actual target and then mask other sources      # #
delta <- 50                                                                  # #
while(delta > dims_bkg[1] | delta > dims_bkg[2]){                            # #
  delta <- delta - 10                                                        # #
}                                                                            # #
delta <- delta / 2                                                           # #
                                                                             # #
center_xs <- -delta:delta + round_del_xy[1]                                  # #
center_ys <- -delta:delta + round_del_xy[2]                                  # #
max_center <- max(bkgI_crp_h[center_xs, center_ys])                          # #
                                                                             # #
bkgI_crp_h[which(bkgI_crp_h > max_center, arr.ind = T)] <- NA                # #
rm(delta, center_xs, center_ys, max_center)                                  # #
                                                                             # #
acq_crp_out_fits <- paste0(main_out, basename, "-Merged-ACQ.fits")           # #
if(!file.exists(acq_crp_out_fits)){                                          # #
  acq_mask[which(acq_mask == 1, arr.ind = T)] <- NA                          # #
  acq_mask[which(acq_mask == 0, arr.ind = T)] <- 1                           # #
                                                                             # #
  acq_masked_data <- acq_data * acq_mask                                     # #
  acq_data <- (acq_data - median(acq_masked_data, na.rm = T)) / acq_expt     # #
  rm(acq_masked_data)                                                        # #
                                                                             # #
  acq_mask[which(acq_mask == 1, arr.ind = T)] <- 0                           # #
  acq_mask[which(is.na(acq_mask), arr.ind = T)] <- 1                         # #
                                                                             # #
  acq_data <- acq_data * acq_mask                                            # #
  acq_crp <- crop_array(acq_data, round_acq_xy, round_del_xy)                # #
                                                                             # #
  acq_crp[which(acq_crp < 0, arr.ind = T)] <- 0                              # #
                                                                             # #
  cpix <- round_del_xy                                                       # #
  cval <- round(acq_val / dimsObj * dims_bkg[1:2])                           # #
  ctyp <- c(typex, typey)                                                    # #
                                                                             # #
  writeFITSim(acq_crp, acq_crp_out_fits, crpixn = cpix, crvaln = cval,       # #
              ctypen = ctyp, header = acq_head)                              # #
  rm(acq_crp)                                                                # #
}                                                                            # #
rm(acq_crp_out_fits, acq_expt, acq_val, acq_mask, acq_data)                  # #
################################################################################

############ Iterating over bands because there is vertical shift ############ #
################################################################################
for(b in bands){                                                             # #
  if(!b_files_exist[b]){                                                     # #
    next                                                                     # #
  }                                                                          # #
                                                                             # #
  galP_path <- galPs[b]                                                      # #
  galX_path <- galXs[b]                                                      # #
  galIP_path <- galIPs[b]                                                    # #
  galQ_path <- galQs[b]                                                      # #
  galU_path <- galUs[b]                                                      # #
  galP_uncut_path <- galPs_uncut[b]                                          # #
  galX_uncut_path <- galXs_uncut[b]                                          # #
  galIP_uncut_path <- galIPs_uncut[b]                                        # #
  galQ_uncut_path <- galQs_uncut[b]                                          # #
  galU_uncut_path <- galUs_uncut[b]                                          # #
  obsI_path <- obsIs[b]                                                      # #
  skyP_path <- skyPs[b]                                                      # #
  skyX_path <- skyXs[b]                                                      # #
  skyIP_path <- skyIPs[b]                                                    # #
  skyI_path <- skyIs[b]                                                      # #
                                                                             # #
  astro_path <- astro_files[b]                                               # #
                                                                             # #
  ############################ CREATING DATA MASK ########################## # #
  ##############################################################################
  print(paste0("Creating aperture mask for ", b," band files..."))         # # #
                                                                           # # #
  temp_fits <- readFITS(galP_path)                                         # # #
  temp_head <- temp_fits$header[-c(8:10)]                                  # # #
                                                                           # # #
  temp_head <- modVal('EXTNAME', 'CROPED', "Extension name", temp_head)    # # #
  temp_head <- modVal('FILETYP1', 'CROP BUNDLE', "File type", temp_head)   # # #
  temp_head <- modVal('SOURCE', 'crop_MERGED.fits',                        # # #
                      "Script used to generate this file ", temp_head)     # # #
  temp_head <- delKwv('FILETYP2', temp_head)                               # # #
  binSt <- as.numeric(get_fits_header_str(temp_head, "BINNING"))           # # #
                                                                           # # #
  #Calculate scale relationship between pixel unit and fov angular unit    # # #
  binx <- get_fits_header_num(temp_head, BINX)                             # # #
  biny <- get_fits_header_num(temp_head, BINY)                             # # #
  pxscl <- get_fits_header_num(temp_head, PIXSCALE)                        # # #
  scl_x <- pxscl * binx / 3600                                             # # #
  scl_y <- pxscl * biny / 3600                                             # # #
                                                                           # # #
  if(binSt != binS){                                                       # # #
    print("ERROR: Binning step of files does not match folder name.")      # # #
    print(paste0("Folder name: ", st_folder))                              # # #
    print(paste0("Binning step of files: ", binSt))                        # # #
    print("The script will be stopped.")                                   # # #
    stop()                                                                 # # #
  }                                                                        # # #
                                                                           # # #
  # Load Coordinates of Ref Pixel                                          # # #
  ref_xy <- c(get_fits_header_num(temp_head, REFX),                        # # #
              get_fits_header_num(temp_head, REFY))                        # # #
                                                                           # # #
  nb <- which(bands == b, arr.ind = T)                                     # # #
  if(wcs_adapt_flag[nb]){                                                  # # #
    ref_coords <- c(get_fits_header_num(temp_head, REFRA),                 # # #
                    get_fits_header_num(temp_head, REFDEC))                # # #
    t_b <- which(bands == bands_adapt[nb], arr.ind = T)                    # # #
                                                                           # # #
    t_fits_head <- readFITS(galPs[t_b])$header[-c(8:10)]                   # # #
                                                                           # # #
    t_ref_xy <- c(get_fits_header_num(t_fits_head, REFX),                  # # #
                  get_fits_header_num(t_fits_head, REFY))                  # # #
    t_ref_coords <- c(get_fits_header_num(t_fits_head, REFRA),             # # #
                      get_fits_header_num(t_fits_head, REFDEC))            # # #
                                                                           # # #
    dx <- t_ref_xy[1] - ref_xy[1]                                          # # #
    dy <- t_ref_xy[2] - ref_xy[2]                                          # # #
                                                                           # # #
    cosd <- cos(mean(c(t_ref_coords[2], ref_coords[2])) * pi / 180)        # # #
                                                                           # # #
    dra_x <- (t_ref_coords[1] - ref_coords[1]) * cosd / scl_x              # # #
    ddec_y <- (t_ref_coords[2] - ref_coords[2]) / scl_y                    # # #
    difx <- dx + dra_x                                                     # # #
    dify <- dy + ddec_y                                                    # # #
    if(og_obj == "NGC-3351"){                                              # # #
      difx <- difx + 17                                                    # # #
      dify <- dify - 45                                                    # # #
    }                                                                      # # #
                                                                           # # #
    if(og_obj == "NGC-3244"){                                              # # #
      difx <- difx - 4                                                     # # #
    }                                                                      # # #
    rm(dx, dy, dra_x, ddec_y)                                              # # #
                                                                           # # #
    ref_dif <- c(difx, dify)                                               # # #
  }else{                                                                   # # #
    ref_dif <- c(0, 0)                                                     # # #
  }                                                                        # # #
  ref_coords <- convert_pixel_to_wcs(array(ref_xy + ref_dif, dim = c(1,2)),# # #
                                     astro_path)                           # # #
  ref_val <- c(get_fits_header_num(temp_head, VALX),                       # # #
               get_fits_header_num(temp_head, VALY))                       # # #
                                                                           # # #
  tgt_xy <- convert_wcs_to_pixel(array(tgt_coords, dim = c(1,2)),          # # #
                                 astro_path) - ref_dif                     # # #
  cosd <- cos(mean(c(ref_coords[2], tgt_coords[2])) * pi / 180)            # # #
                                                                           # # #
  ctyp <- c(typex, typey, typez, typez2)                                   # # #
                                                                           # # #
  for(xa in 1:dimsObj[1]){                                                 # # #
    for(ya in 1:dimsObj[2]){                                               # # #
      rx <- ((xa - tgt_xy[1]) * cos(tgt_pa) + (ya - tgt_xy[2]) *           # # #
               sin(tgt_pa))^2 / tgt_a^2                                    # # #
      ry <- ((xa - tgt_xy[1]) * sin(tgt_pa) - (ya - tgt_xy[2]) *           # # #
               cos(tgt_pa))^2 / tgt_b^2                                    # # #
                                                                           # # #
      if(rx + ry <= 1){                                                    # # #
        mask[xa, ya] <- 1                                                  # # #
      }                                                                    # # #
    }                                                                      # # #
  }                                                                        # # #
  rm(xa, ya, rx, ry)                                                       # # #
                                                                           # # #
  mask[which(mask == 0, arr.ind = T)] <- NA                                # # #
                                                                           # # #
  round_tgt_xy <- round(c(tgt_xy[1], tgt_xy[2]))                           # # #
                                                                           # # #
  max_x <- max(which(mask == 1, arr.ind = T)[,1])                          # # #
  min_x <- min(which(mask == 1, arr.ind = T)[,1])                          # # #
  max_y <- max(which(mask == 1, arr.ind = T)[,2])                          # # #
  min_y <- min(which(mask == 1, arr.ind = T)[,2])                          # # #
                                                                           # # #
  delX <- max(abs(c(max_x, min_x) - round_tgt_xy[1]))                      # # #
  delY <- max(abs(c(max_y, min_y) - round_tgt_xy[2]))                      # # #
                                                                           # # #
  round_del_xy <- round(c(delX, delY))                                     # # #
  rm(delX, delY)                                                           # # #
                                                                           # # #
  ##############################################################################
                                                                             # #
  g_name <- paste0(basename, b)                                              # #
  s_name <- paste0(basename, b, "_Sky")                                      # #
                                                                             # #
  print(paste0("Loading ", b, " band files..."))                             # #
  gP <- readFITS(galP_path)                                                  # #
  gX <- readFITS(galX_path)                                                  # #
  gIP <- readFITS(galIP_path)                                                # #
  gQ <- readFITS(galQ_path)                                                  # #
  gU <- readFITS(galU_path)                                                  # #
  gP_uncut <- readFITS(galP_uncut_path)                                      # #
  gX_uncut <- readFITS(galX_uncut_path)                                      # #
  gIP_uncut <- readFITS(galIP_uncut_path)                                    # #
  gQ_uncut <- readFITS(galQ_uncut_path)                                      # #
  gU_uncut <- readFITS(galU_uncut_path)                                      # #
  oI <- readFITS(obsI_path)                                                  # #
  sP <- readFITS(skyP_path)                                                  # #
  sX <- readFITS(skyX_path)                                                  # #
  sIP <- readFITS(skyIP_path)                                                # #
  sI <- readFITS(skyI_path)                                                  # #
                                                                             # #
  gP_img <- gP$imDat                                                         # #
  gX_img <- gX$imDat                                                         # #
  gIP_img <- gIP$imDat                                                       # #
  gQ_img <- gQ$imDat                                                         # #
  gU_img <- gU$imDat                                                         # #
  gP_uncut_img <- gP_uncut$imDat                                             # #
  gX_uncut_img <- gX_uncut$imDat                                             # #
  gIP_uncut_img <- gIP_uncut$imDat                                           # #
  gQ_uncut_img <- gQ_uncut$imDat                                             # #
  gU_uncut_img <- gU_uncut$imDat                                             # #
  oI_img <- oI$imDat                                                         # #
  sP_img <- sP$imDat                                                         # #
  sX_img <- sX$imDat                                                         # #
  sIP_img <- sIP$imDat                                                       # #
  sI_img <- sI$imDat                                                         # #
                                                                             # #
  print(paste0("Applying aperture mask to ", b, " band files..."))           # #
  gP_img[,,1] <- gP_img[,,1] * mask                                          # #
  gX_img[,,1] <- gX_img[,,1] * mask                                          # #
  gIP_img[,,1] < gIP_img[,,1] * mask                                         # #
  gQ_img[,,1] <- gQ_img[,,1] * mask                                          # #
  gU_img[,,1] <- gU_img[,,1] * mask                                          # #
  gP_uncut_img[,,1] <- gP_uncut_img[,,1] * mask                              # #
  gX_uncut_img[,,1] <- gX_uncut_img[,,1] * mask                              # #
  gIP_uncut_img[,,1] < gIP_uncut_img[,,1] * mask                             # #
  gQ_uncut_img[,,1] <- gQ_uncut_img[,,1] * mask                              # #
  gU_uncut_img[,,1] <- gU_uncut_img[,,1] * mask                              # #
  oI_img[,,1] <- oI_img[,,1] * mask                                          # #
  gP_img[,,2] <- gP_img[,,2] * mask                                          # #
  gX_img[,,2] <- gX_img[,,2] * mask                                          # #
  gIP_img[,,2] < gIP_img[,,2] * mask                                         # #
  gQ_img[,,2] <- gQ_img[,,2] * mask                                          # #
  gU_img[,,2] <- gU_img[,,2] * mask                                          # #
  gP_uncut_img[,,2] <- gP_uncut_img[,,2] * mask                              # #
  gX_uncut_img[,,2] <- gX_uncut_img[,,2] * mask                              # #
  gIP_uncut_img[,,2] < gIP_uncut_img[,,2] * mask                             # #
  gQ_uncut_img[,,2] <- gQ_uncut_img[,,2] * mask                              # #
  gU_uncut_img[,,2] <- gU_uncut_img[,,2] * mask                              # #
  oI_img[,,2] <- oI_img[,,2] * mask                                          # #
  sP_img[,,1] <- sP_img[,,1] * mask                                          # #
  sX_img[,,1] <- sX_img[,,1] * mask                                          # #
  sIP_img[,,1] < sIP_img[,,1] * mask                                         # #
  sIM_img <- sI_img[,,1]                                                     # #
  sI_img[,,1] <- sI_img[,,1] * mask                                          # #
  sP_img[,,2] <- sP_img[,,2] * mask                                          # #
  sX_img[,,2] <- sX_img[,,2] * mask                                          # #
  sIP_img[,,2] < sIP_img[,,2] * mask                                         # #
  sI_img[,,2] <- sI_img[,,2] * mask                                          # #
                                                                             # #
  ############################### CROPING DATA ############################# # #
  ##############################################################################
  print(paste0("Cropping masked ", b, " band files..."))                   # # #
                                                                           # # #
  gP_crp <- crop_array(gP_img, round_tgt_xy, round_del_xy)                 # # #
  gX_crp <- crop_array(gX_img, round_tgt_xy, round_del_xy)                 # # #
  gIP_crp <- crop_array(gIP_img, round_tgt_xy, round_del_xy)               # # #
  gQ_crp <- crop_array(gQ_img, round_tgt_xy, round_del_xy)                 # # #
  gU_crp <- crop_array(gU_img, round_tgt_xy, round_del_xy)                 # # #
  gP_uncut_crp <- crop_array(gP_uncut_img, round_tgt_xy, round_del_xy)     # # #
  gX_uncut_crp <- crop_array(gX_uncut_img, round_tgt_xy, round_del_xy)     # # #
  gIP_uncut_crp <- crop_array(gIP_uncut_img, round_tgt_xy, round_del_xy)   # # #
  gQ_uncut_crp <- crop_array(gQ_uncut_img, round_tgt_xy, round_del_xy)     # # #
  gU_uncut_crp <- crop_array(gU_uncut_img, round_tgt_xy, round_del_xy)     # # #
  oI_crp <- crop_array(oI_img, round_tgt_xy, round_del_xy)                 # # #
  sP_crp <- crop_array(sP_img, round_tgt_xy, round_del_xy)                 # # #
  sX_crp <- crop_array(sX_img, round_tgt_xy, round_del_xy)                 # # #
  sIP_crp <- crop_array(sIP_img, round_tgt_xy, round_del_xy)               # # #
  sI_crp <- crop_array(sI_img, round_tgt_xy, round_del_xy)                 # # #
  msk_crp <- crop_array(mask, round_tgt_xy, round_del_xy)                  # # #
  rm(gP_img, gX_img, gIP_img, gP_uncut_img, gX_uncut_img, gIP_uncut_img,   # # #
     oI_img, sP_img, sX_img, sI_img)                                       # # #
                                                                           # # #
  astro_out <- paste0(main_out, g_name, "_WCS.fits")                       # # #
                                                                           # # #
  if(!file.exists(astro_out)){                                             # # #
    astro <- readFITS(astro_path)                                          # # #
    astro_img <- astro$imDat                                               # # #
    astro_head <- astro$header[-(10:30)]                                   # # #
                                                                           # # #
    astro_ref <- c(get_fits_header_num(astro_head, 'CRPIX1'),              # # #
                   get_fits_header_num(astro_head, 'CRPIX2'))              # # #
                                                                           # # #
    wcs_tgt_xy <- round(c(tgt_xy[1] + ref_dif[1], tgt_xy[2] + ref_dif[2])) # # #
    yi <- dimsObj[1] - wcs_tgt_xy[1] - 1 - round_del_xy[1]                 # # #
    yf <- dimsObj[1] - wcs_tgt_xy[1] - 1 + round_del_xy[1]                 # # #
    xi <- wcs_tgt_xy[2] - 1 - round_del_xy[2]                              # # #
    xf <- wcs_tgt_xy[2] - 1 + round_del_xy[2]                              # # #
                                                                           # # #
    dif_ref <- astro_ref - wcs_tgt_xy                                      # # #
    new_refx <- round_del_xy[1] + 1 + dif_ref[1]                           # # #
    new_refy <- round_del_xy[2] + 1 + dif_ref[2]                           # # #
                                                                           # # #
    py_run_string("output_path = r.astro_out")                             # # #
    py_run_string("hdu = fits.open(r.astro_path)")                         # # #
                                                                           # # #
    # Crop data                                                            # # #
    crop_cmd <- paste0("hdu[0].data = hdu[0].data[int(r.xi):int(r.xf),",   # # #
                          " int(r.yi):int(r.yf)]")                         # # #
    py_run_string(crop_cmd)                                                # # #
                                                                           # # #
    # Corrigir pixel de referência no header                               # # #
    py_run_string("hdu[0].header['CRPIX1'] = float(r.new_refx)")           # # #
    py_run_string("hdu[0].header['CRPIX2'] = float(r.new_refy)")           # # #
                                                                           # # #
    # Save                                                                 # # #
    save_cmd <- paste0("hdu.writeto(output_path, overwrite=True, ",        # # #
                       "output_verify='ignore')")                          # # #
    py_run_string(save_cmd)                                                # # #
    rm(save_cmd, crop_cmd, astro, astro_img, astro_head, astro_ref,        # # #
       xi, xf, yi, yf, dif_ref, new_refx, new_refy)                        # # #
  }                                                                        # # #
                                                                           # # #
  dims_crp <- dim(gP_crp[,,1:2])                                           # # #
  crp_tgt_xy <- round_del_xy + 1                                           # # #
  crp_tgt_val <- round(ref_val / dimsObj * dims_crp[1:2])                  # # #
  rm(round_del_xy, round_tgt_xy)                                           # # #
                                                                           # # #
  ghold <- array(NA, dim = c(dims_crp, 6),                                 # # #
                 dimnames = list(NULL, NULL, c("Data", "Unc"),             # # #
                                 c("I", "P", "<IP>", "X", "Q", "U")))      # # #
  ghold_uncut <- array(NA, dim = c(dims_crp, 6),                           # # #
                       dimnames = list(NULL, NULL, c("Data", "Unc"),       # # #
                                       c("I", "P", "<IP>", "X", "Q", "U")))# # #
  ghold[,,1,1] <- oI_crp[,,1] - sI_crp[,,1]                                # # #
  ghold[,,2,1] <- unc_add(oI_crp[,,2], sI_crp[,,2])                        # # #
  ghold_uncut[,,1,1] <- ghold[,,1,1]                                       # # #
  ghold_uncut[,,2,1] <- ghold[,,2,1]                                       # # #
  ##############################################################################
                                                                             # #
  #------------------------------------------------------------------------# # #
                                                                             # #
  ################## CHECKING FOR/CREATING CONTOUR ARRAY ################### # #
  ##############################################################################
  iso_com <- "contours.fits"                                                   #
  iso_paths <- list.files(main_out, full.names = T)                            #
  iso_file <- iso_paths[grep(iso_com, iso_paths)]                              #
  iso_flag <- length(iso_file) == 1                                            #
                                                                               #
  if(!iso_flag){                                                               #
    ctr_com <- "contours.ctr"                                                  #
    ctr_file <- iso_paths[grep(ctr_com, iso_paths)]                            #
    ctr_flag <- length(ctr_file) == 1                                          #
    iso_file <- paste0(main_out, basename, iso_com)                            #
                                                                               #
    if(ctr_flag){                                                              #
      py_run_string("max_x = r.dims_bkg[0]")                                   #
      py_run_string("max_y = r.dims_bkg[1]")                                   #
      py_run_string("ctr_path = r.ctr_file")                                   #
      py_run_string("iso_fits_path = r.iso_file")                              #
                                                                               #
      py_run_string("ctr_data = sexD.load_ctr_file(ctr_path)")                 #
      my_arr_cmd <- paste0("my_arr = sexD.create_2d_array_from_ctr(ctr_data,", #
                           "max_x,max_y)")                                     #
      py_run_string(my_arr_cmd)                                                #
      py_run_string("sexD.save_arr_as_fits(my_arr, iso_fits_path)")            #
      py_run_string("del(max_x, max_y, ctr_path, iso_fits_path, ctr_data)")    #
      py_run_string("del(my_arr)")                                             #
                                                                               #
      ctr_arr_h <- readFITS(iso_file)$imDat                                    #
      iso_flag <- TRUE                                                         #
    }else{                                                                     #
      print(paste0("There are no '_contours.fits' or '_contours.ctr' files fo",#
                   "r this object. Use the cropped acq file, outputted by thi",#
                   "s code to create '_contours.ctr' and '_contours.lev' file",#
                   "s in SAO DS9. Save them in the same folder you found the ",#
                   "cropped acq file and re-run this code."))                  #
    }                                                                          #
  }else{                                                                       #
    ctr_arr_h <- readFITS(iso_file)$imDat                                      #
  }                                                                            #
  ##############################################################################
                                                                             # #
  #------------------------------------------------------------------------# # #
                                                                             # #
  ######################### COMPRESSING BINNED DATA ######################## # #
  ##############################################################################
  print(paste0("Binning Gal ", b, " band data..."))                        # # #
  gP_bin <- bin_shrink_map(gP_crp, binS, crp_tgt_xy)                       # # #
  gX_bin <- bin_shrink_map(gX_crp, binS, crp_tgt_xy)                       # # #
  gIP_bin <- bin_shrink_map(gIP_crp, binS, crp_tgt_xy)                     # # #
  gQ_bin <- bin_shrink_map(gQ_crp, binS, crp_tgt_xy)                       # # #
  gU_bin <- bin_shrink_map(gU_crp, binS, crp_tgt_xy)                       # # #
  gP_uncut_bin <- bin_shrink_map(gP_uncut_crp, binS, crp_tgt_xy)           # # #
  gX_uncut_bin <- bin_shrink_map(gX_uncut_crp, binS, crp_tgt_xy)           # # #
  gIP_uncut_bin <- bin_shrink_map(gIP_uncut_crp, binS, crp_tgt_xy)         # # #
  gQ_uncut_bin <- bin_shrink_map(gP_uncut_crp, binS, crp_tgt_xy)           # # #
  gU_uncut_bin <- bin_shrink_map(gP_uncut_crp, binS, crp_tgt_xy)           # # #
  gI_bin <- bin_shrink_map(ghold[,,,1], binS, crp_tgt_xy)                  # # #
                                                                           # # #
  print(paste0("Binning Sky ", b, " band data..."))                        # # #
  sP_bin <- bin_shrink_map(sP_crp, binS, crp_tgt_xy)                       # # #
  sX_bin <- bin_shrink_map(sX_crp, binS, crp_tgt_xy)                       # # #
  sIP_bin <- bin_shrink_map(sIP_crp, binS, crp_tgt_xy)                     # # #
  sI_bin <- bin_shrink_map(sI_crp, binS, crp_tgt_xy)                       # # #
  sP_crpS <- bin_expand_map(sP_bin, binS, dims_crp, crp_tgt_xy, dense = F) # # #
  sX_crpS <- bin_expand_map(sX_bin, binS, dims_crp, crp_tgt_xy, dense = F) # # #
  sIP_crpS <- bin_expand_map(sIP_bin, binS, dims_crp, crp_tgt_xy, dense =F)# # #
                                                                           # # #
  if(dims_bkg[1] != dims_crp[1] | dims_bkg[2] != dims_crp[2]){             # # #
    bkgI_crp <- array(0, dim = dims_crp[1:2])                              # # #
    ctr_arr <- array(0, dim = dims_crp[1:2])                               # # #
                                                                           # # #
    delta_dim <- dims_bkg - dims_crp[1:2]                                  # # #
                                                                           # # #
    del_x <- floor(delta_dim[1] / 2)                                       # # #
    del_y <- floor(delta_dim[2] / 2)                                       # # #
                                                                           # # #
    # Cutting or padding depending on delta_dim sign                       # # #
    if(del_x > 0){                                                         # # #
      xs_h <- (1 + del_x):(dims_bkg[1] - del_x)                            # # #
                                                                           # # #
      if(delta_dim[1] %% 2 == 1){                                          # # #
        xs_h <- xs_h[-dims_crp[1]]                                         # # #
      }                                                                    # # #
                                                                           # # #
      xs_n <- 1:dims_crp[1]                                                # # #
    }else{                                                                 # # #
      xs_h <- 1:dims_bkg[1]                                                # # #
      xs_n <- (1 - del_x):(dims_crp[1] + del_x)                            # # #
                                                                           # # #
      if(delta_dim[1] %% 2 == 1){                                          # # #
        xs_n <- c(xs_n, xs_n[dims_crop[1]])                                # # #
      }                                                                    # # #
    }                                                                      # # #
    if(del_y > 0){                                                         # # #
      ys_h <- (1 + del_y):(dims_bkg[2] - del_y)                            # # #
                                                                           # # #
      if(delta_dim[2] %% 2 == 1){                                          # # #
        ys_h <- ys_h[-dims_crp[2]]                                         # # #
      }                                                                    # # #
                                                                           # # #
      ys_n <- 1:dims_crp[2]                                                # # #
    }else{                                                                 # # #
      ys_h <- 1:dims_bkg[2]                                                # # #
      ys_n <- (1 - del_y):(dims_crp[2] + del_y)                            # # #
                                                                           # # #
      if(delta_dim[2] %% 2 == 1){                                          # # #
        ys_n <- c(ys_n, ys_n[dims_crp[2]] + 1)                             # # #
      }                                                                    # # #
    }                                                                      # # #
                                                                           # # #
    bkgI_crp[xs_n, ys_n] <- bkgI_crp_h[xs_h, ys_h]                         # # #
    ctr_arr[xs_n, ys_n] <- ctr_arr_h[xs_h, ys_h]                           # # #
  }else{                                                                   # # #
    bkgI_crp <- bkgI_crp_h                                                 # # #
    ctr_arr <- ctr_arr_h                                                   # # #
  }                                                                        # # #
  ##############################################################################
                                                                             # #
  #------------------------------------------------------------------------# # #
                                                                             # #
  ################### CREATING HISTOGRAMS OF SNR AND NSR ################### # #
  ##############################################################################
  print(paste0("Creating SNR and NSR histograms for Gal ", b," band..."))  # # #
                                                                           # # #
  pdf_name <- paste0(out_folder, g_name, "_SNR_NSR_gal_hist.pdf")          # # #
  pdf(pdf_name, width = 31, height = 31)                                   # # #
  par(mfrow = c(2, 2))                                                     # # #
  par(cex.main = 3)                                                        # # #
  par(cex.axis = 2)                                                        # # #
  par(cex.lab = 2)                                                         # # #
                                                                           # # #
  # dI / I                                                                 # # #
  h_temp <- abs(gI_bin[,,2] / gI_bin[,,1])                                 # # #
  h_temp[which(is.infinite(h_temp), arr.ind = T)] <- NA                    # # #
  h_dat <- as.vector(h_temp[which(!is.na(h_temp), arr.ind = T)])           # # #
  rm(h_temp)                                                               # # #
                                                                           # # #
  last_brk <- ceiling(max(h_dat, na.rm = T)) + 0.5                         # # #
  hist_brks <- c(0, seq(0.05, 1.05, 0.1), last_brk)                        # # #
                                                                           # # #
  hist_dat <- hist(h_dat, breaks = hist_brks, plot = FALSE, right = FALSE) # # #
                                                                           # # #
  h_freq <- hist_dat$counts / sum(hist_dat$counts)                         # # #
  h_freq[length(h_freq) - 1] <-                                            # # #
    h_freq[length(h_freq) - 1] + h_freq[length(h_freq)]                    # # #
  h_freq <- h_freq[-length(h_freq)]                                        # # #
                                                                           # # #
  h_brks <- hist_dat$breaks[-length(hist_brks)]                            # # #
                                                                           # # #
  h_med <- median(h_dat, na.rm = T)                                        # # #
                                                                           # # #
  diffs <- abs(h_brks - h_med)                                             # # #
  diffs_min <- min(diffs)                                                  # # #
                                                                           # # #
  brk_i <- which(diffs == diffs_min, arr.ind = T)                          # # #
                                                                           # # #
  if(h_brks[brk_i] <= h_med){                                              # # #
    brk_w <- h_brks[brk_i + 1] - h_brks[brk_i]                             # # #
    med_pos <- (h_med - h_brks[brk_i]) / brk_w + brk_i - 1                 # # #
  }else{                                                                   # # #
    brk_w <- h_brks[brk_i] - h_brks[brk_i - 1]                             # # #
    med_pos <- brk_i - (h_brks[brk_i] - h_med) / brk_w - 1                 # # #
  }                                                                        # # #
                                                                           # # #
  barplot(h_freq, space = 0, xlab = "sI/I", ylab = "Frequency",            # # #
          main = "sI/I gal", ylim = c(0, max(h_freq) * 1.25), col = "red") # # #
  axis(1, at = seq_along(h_brks) - 1, labels = h_brks, las = 1)            # # #
  abline(v = med_pos, col = "blue", type = "h")                            # # #
                                                                           # # #
  # I / dI                                                                 # # #
  h_dat <- 1 / h_dat                                                       # # #
  h_dat[which(is.infinite(h_dat), arr.ind = T)] <- NA                      # # #
  last_brk <- ceiling(max(h_dat, na.rm = T))                               # # #
  hist_brks <- rev(1 / hist_brks)                                          # # #
  hist_brks[length(hist_brks)] <- last_brk                                 # # #
                                                                           # # #
  hist_dat <- hist(h_dat, breaks = hist_brks, plot = FALSE, right = FALSE) # # #
                                                                           # # #
  h_freq <- hist_dat$counts / sum(hist_dat$counts)                         # # #
  h_freq[2] <- h_freq[2] + h_freq[1]                                       # # #
  h_freq <- h_freq[-1]                                                     # # #
                                                                           # # #
  h_brks <- hist_dat$breaks[-1]                                            # # #
                                                                           # # #
  h_med <- median(h_dat, na.rm = T)                                        # # #
                                                                           # # #
  diffs <- abs(h_brks - h_med)                                             # # #
  diffs_min <- min(diffs)                                                  # # #
                                                                           # # #
  brk_i <- which(diffs == diffs_min, arr.ind = T)                          # # #
                                                                           # # #
  if(h_brks[brk_i] <= h_med){                                              # # #
    brk_w <- h_brks[brk_i + 1] - h_brks[brk_i]                             # # #
    med_pos <- (h_med - h_brks[brk_i]) / brk_w + brk_i - 1                 # # #
  }else{                                                                   # # #
    brk_w <- h_brks[brk_i] - h_brks[brk_i - 1]                             # # #
    med_pos <- brk_i - (h_brks[brk_i] - h_med) / brk_w - 1                 # # #
  }                                                                        # # #
                                                                           # # #
  barplot(h_freq, space = 0, xlab = "I/sI", ylab = "Frequency",            # # #
          main = "I/sI gal", ylim = c(0, max(h_freq) * 1.25), col = "red") # # #
  axis(1, at = seq_along(h_brks) - 1, labels = round(h_brks, 2), las = 1)  # # #
  abline(v = med_pos, col = "blue", type = "h")                            # # #
                                                                           # # #
  # dP / P                                                                 # # #
  h_temp <- abs(gP_bin[,,2] / gP_bin[,,1])                                 # # #
  h_temp[which(is.infinite(h_temp), arr.ind = T)] <- NA                    # # #
  h_dat <- as.vector(h_temp[which(!is.na(h_temp), arr.ind = T)])           # # #
  rm(h_temp)                                                               # # #
                                                                           # # #
  last_brk <- ceiling(max(h_dat, na.rm = T)) + 0.5                         # # #
  hist_brks <- c(0, seq(0.05, 1.05, 0.1), last_brk)                        # # #
                                                                           # # #
  hist_dat <- hist(h_dat, breaks = hist_brks, plot = FALSE, right = FALSE) # # #
                                                                           # # #
  h_freq <- hist_dat$counts / sum(hist_dat$counts)                         # # #
  h_freq[length(h_freq) - 1] <-                                            # # #
    h_freq[length(h_freq) - 1] + h_freq[length(h_freq)]                    # # #
  h_freq <- h_freq[-length(h_freq)]                                        # # #
                                                                           # # #
  h_brks <- hist_dat$breaks[-length(hist_brks)]                            # # #
                                                                           # # #
  h_med <- median(h_dat, na.rm = T)                                        # # #
                                                                           # # #
  diffs <- abs(h_brks - h_med)                                             # # #
  diffs_min <- min(diffs)                                                  # # #
                                                                           # # #
  brk_i <- which(diffs == diffs_min, arr.ind = T)                          # # #
                                                                           # # #
  if(h_brks[brk_i] <= h_med){                                              # # #
    brk_w <- h_brks[brk_i + 1] - h_brks[brk_i]                             # # #
    med_pos <- (h_med - h_brks[brk_i]) / brk_w + brk_i - 1                 # # #
  }else{                                                                   # # #
    brk_w <- h_brks[brk_i] - h_brks[brk_i - 1]                             # # #
    med_pos <- brk_i - (h_brks[brk_i] - h_med) / brk_w - 1                 # # #
  }                                                                        # # #
                                                                           # # #
  barplot(h_freq, space = 0, xlab = "sP/P", ylab = "Frequency",            # # #
          main = "sP/P gal", ylim = c(0, max(h_freq) * 1.25), col = "red") # # #
  axis(1, at = seq_along(h_brks) - 1, labels = h_brks, las = 1)            # # #
  abline(v = med_pos, col = "blue", type = "h")                            # # #
                                                                           # # #
  # P / dP                                                                 # # #
  h_dat <- 1 / h_dat                                                       # # #
  h_dat[which(is.infinite(h_dat), arr.ind = T)] <- NA                      # # #
  last_brk <- ceiling(max(h_dat, na.rm = T))                               # # #
  hist_brks <- rev(1 / hist_brks)                                          # # #
  hist_brks[length(hist_brks)] <- last_brk                                 # # #
                                                                           # # #
  hist_dat <- hist(h_dat, breaks = hist_brks, plot = FALSE, right = FALSE) # # #
                                                                           # # #
  h_freq <- hist_dat$counts / sum(hist_dat$counts)                         # # #
  h_freq[2] <- h_freq[2] + h_freq[1]                                       # # #
  h_freq <- h_freq[-1]                                                     # # #
                                                                           # # #
  h_brks <- hist_dat$breaks[-1]                                            # # #
                                                                           # # #
  h_med <- median(h_dat, na.rm = T)                                        # # #
                                                                           # # #
  diffs <- abs(h_brks - h_med)                                             # # #
  diffs_min <- min(diffs)                                                  # # #
                                                                           # # #
  brk_i <- which(diffs == diffs_min, arr.ind = T)                          # # #
                                                                           # # #
  if(h_brks[brk_i] <= h_med){                                              # # #
    brk_w <- h_brks[brk_i + 1] - h_brks[brk_i]                             # # #
    med_pos <- (h_med - h_brks[brk_i]) / brk_w + brk_i - 1                 # # #
  }else{                                                                   # # #
    brk_w <- h_brks[brk_i] - h_brks[brk_i - 1]                             # # #
    med_pos <- brk_i - (h_brks[brk_i] - h_med) / brk_w - 1                 # # #
  }                                                                        # # #
                                                                           # # #
  barplot(h_freq, space = 0, xlab = "P/sP", ylab = "Frequency",            # # #
          main = "P/sP gal", ylim = c(0, max(h_freq) * 1.25), col = "red") # # #
  axis(1, at = seq_along(h_brks) - 1, labels = round(h_brks, 2), las = 1)  # # #
  abline(v = med_pos, col = "blue", type = "h")                            # # #
  dev.off()                                                                # # #
  rm(gI_bin, sI_bin)                                                       # # #
  ##############################################################################
                                                                             # #
  #------------------------------------------------------------------------# # #
                                                                             # #
  ######################### CREATING P/uP vs uP PLOT ####################### # #
  ##############################################################################
  print(paste0("Performing statiscal analysis of P/σP for Gal ", b,        # # #
               " band..."))                                                # # #
                                                                           # # #
  temp_snr <- as.vector(gP_bin[,,1] / gP_bin[,,2])                         # # #
  temp_snr[which(is.infinite(temp_snr), arr.ind = T)] <- NA                # # #
  temp_NAs <- which(is.na(temp_snr), arr.ind = T)                          # # #
  temp_snr <- temp_snr[-temp_NAs]                                          # # #
  temp_n <- as.vector(gP_bin[,,2])[-temp_NAs]                              # # #
                                                                           # # #
  temp_g <- gP_bin[,,1]                                                    # # #
  non_na_len <- length(as.vector(temp_g[which(!is.na(temp_g))]))           # # #
  rm(temp_g)                                                               # # #
                                                                           # # #
  # Establishing a set of High quality bins                                # # #
  # This set size should be 1% of total amount of valid bins               # # #
  # or 20 bins, whichever is largest                                       # # #
  # High quality must not be smaller than bare minimum                     # # #
  min_snr <- 1.5                                                           # # #####
  snr_step <- 0.1                                                          # # #
  high_snr <- 10                                                           # # #
  temp_HIGHs <- which(temp_snr >= high_snr, arr.ind = T)                   # # #
  min_rep <- 0.05                                                          # # ###0.01
  min_num <- 35                                                            # # ###20
  while(length(temp_HIGHs) < (min_rep * non_na_len) &                      # # #
        (high_snr - min_snr) > snr_step){                                  # # #
    high_snr <- high_snr - snr_step                                        # # #
    temp_HIGHs <- which(temp_snr >= high_snr, arr.ind = T)                 # # #
  }                                                                        # # #
  # If 1% of total amount is less than 20 bins, then minimum size is 20    # # #
  # High quality must not be smaller than bare minimum                     # # #
  while(length(temp_HIGHs) < min_num & (high_snr - min_snr) > snr_step){   # # #
    high_snr <- high_snr - snr_step                                        # # #
    temp_HIGHs <- which(temp_snr >= high_snr, arr.ind = T)                 # # #
  }                                                                        # # #
                                                                           # # #
  mean_n_HIGHs <- median(temp_n[temp_HIGHs])                               # # #
  sd_n_HIGHs <- mad(temp_n[temp_HIGHs])                                    # # #
  rm(temp_HIGHs)                                                           # # #
                                                                           # # #
  print(paste0("Creating P/σP vs σP plot for Gal ", b," band..."))         # # #
                                                                           # # #
  pdf_name <- paste0(out_folder, g_name, "_PdσP_vs_σP-gal.pdf")            # # #
  pdf(pdf_name, width = 15, height = 15)                                   # # #
  par(cex.main = 3)                                                        # # #
  par(cex.axis = 2)                                                        # # #
  par(cex.lab = 2)                                                         # # #
                                                                           # # #
  y_max <- max(temp_n)                                                     # # #
  y_min <- min(temp_n)                                                     # # #
  x_max <- max(temp_snr)                                                   # # #
  x_min <- min(temp_snr)                                                   # # #
  h_max <- round(mean_n_HIGHs + 2 * sd_n_HIGHs, 4)                         # # #
  v_sel <- 1                                                               # # #
  v_min <- high_snr                                                        # # #
  rm(temp_HIGHs, mean_n_HIGHs, sd_n_HIGHs)                                 # # #
                                                                           # # #
  matplot(temp_snr, temp_n * 100, xlab = 'P/sP', ylab = 'sP (%)',          # # #
          xlim = c(x_min, x_max), ylim = c(y_min, y_max) * 100, cex = 2,   # # #
          main = paste0("P/sP vs sP (%) for Gal ", g_name), pch = '+')     # # #
  abline(h = h_max * 100, col = "grey", type = "dashed", lwd = 2)          # # #
  abline(v = v_min, col = "blue", type = "solid", lwd = 2)                 # # #
  abline(v = v_sel, col = "red", type = "solid", lwd = 2)                  # # #
  dev.off()                                                                # # #
  rm(temp_snr, temp_n, temp_NAs)                                           # # #
  ##############################################################################
                                                                             # #
  #------------------------------------------------------------------------# # #
                                                                             # #
  ############### FILTERING PIXELS WITH TOO HIGH UNCERTAINTY ############### # #
  ##############################################################################
  print(paste0("Filtering pixels which P/σP is lesser than ", v_min, " wi",# # #
               "th σP higher than ", h_max, " within Gal maps in ", b, " ",# # #
               "band."))                                                   # # #
  TN_not_NA <- length(which(!is.na(gP_bin[,,1]), arr.ind = T))/ 2          # # #
  rej_inds <- which((gP_bin[,,1] / gP_bin[,,2] < v_min &                   # # #
                      gP_bin[,,2] > h_max) | gP_bin[,,1] > 1, arr.ind = T) # # #
                                                                           # # #
  temp <- gP_bin[,,1]                                                      # # #
  temp[rej_inds] <- NA                                                     # # #
  gP_bin[,,1] <- temp                                                      # # #
  temp <- gP_bin[,,2]                                                      # # #
  temp[rej_inds] <- NA                                                     # # #
  gP_bin[,,2] <- temp                                                      # # #
  temp <- gIP_bin[,,1]                                                     # # #
  temp[rej_inds] <- NA                                                     # # #
  gIP_bin[,,1] <- temp                                                     # # #
  temp <- gIP_bin[,,2]                                                     # # #
  temp[rej_inds] <- NA                                                     # # #
  gIP_bin[,,2] <- temp                                                     # # #
  temp <- gX_bin[,,1]                                                      # # #
  temp[rej_inds] <- NA                                                     # # #
  gX_bin[,,1] <- temp                                                      # # #
  temp <- gX_bin[,,2]                                                      # # #
  temp[rej_inds] <- NA                                                     # # #
  gX_bin[,,2] <- temp                                                      # # #
  temp <- gQ_bin[,,1]                                                      # # #
  temp[rej_inds] <- NA                                                     # # #
  gQ_bin[,,1] <- temp                                                      # # #
  temp <- gQ_bin[,,2]                                                      # # #
  temp[rej_inds] <- NA                                                     # # #
  gQ_bin[,,2] <- temp                                                      # # #
  temp <- gU_bin[,,1]                                                      # # #
  temp[rej_inds] <- NA                                                     # # #
  gU_bin[,,1] <- temp                                                      # # #
  temp <- gU_bin[,,2]                                                      # # #
  temp[rej_inds] <- NA                                                     # # #
  gU_bin[,,2] <- temp                                                      # # #
  rm(temp)                                                                 # # #
                                                                           # # #
  N_not_NA <- length(which(!is.na(gP_bin[,,1]), arr.ind = T))/ 2           # # #
                                                                           # # #
  gP_crp <- bin_expand_map(gP_bin, binS, dims_crp, crp_tgt_xy)             # # #
  gIP_crp <- bin_expand_map(gIP_bin, binS, dims_crp, crp_tgt_xy)           # # #
  gX_crp <- bin_expand_map(gX_bin, binS, dims_crp, crp_tgt_xy)             # # #
  gQ_crp <- bin_expand_map(gQ_bin, binS, dims_crp, crp_tgt_xy)             # # #
  gU_crp <- bin_expand_map(gU_bin, binS, dims_crp, crp_tgt_xy)             # # #
  gP_uncut_crp <- bin_expand_map(gP_uncut_bin, binS, dims_crp, crp_tgt_xy) # # #
  gIP_uncut_crp <- bin_expand_map(gIP_uncut_bin, binS, dims_crp,crp_tgt_xy)# # #
  gX_uncut_crp <- bin_expand_map(gX_uncut_bin, binS, dims_crp, crp_tgt_xy) # # #
  gQ_uncut_crp <- bin_expand_map(gQ_uncut_bin, binS, dims_crp, crp_tgt_xy) # # #
  gU_uncut_crp <- bin_expand_map(gU_uncut_bin, binS, dims_crp, crp_tgt_xy) # # #
                                                                           # # #
  gP_crpS <- bin_expand_map(gP_bin, binS, dims_crp, crp_tgt_xy, dense = F) # # #
  gIP_crpS <- bin_expand_map(gIP_bin, binS, dims_crp, crp_tgt_xy, dense =F)# # #
  gX_crpS <- bin_expand_map(gX_bin, binS, dims_crp, crp_tgt_xy, dense = F) # # #
  gQ_crpS <- bin_expand_map(gQ_bin, binS, dims_crp, crp_tgt_xy, dense = F) # # #
  gU_crpS <- bin_expand_map(gU_bin, binS, dims_crp, crp_tgt_xy, dense = F) # # #
  gP_uncut_crpS <- bin_expand_map(gP_uncut_bin, binS, dims_crp,            # # #
                                  crp_tgt_xy, dense = F)                   # # #
  gIP_uncut_crpS <- bin_expand_map(gIP_uncut_bin, binS, dims_crp,          # # #
                                   crp_tgt_xy, dense = F)                  # # #
  gX_uncut_crpS <- bin_expand_map(gX_uncut_bin, binS, dims_crp,            # # #
                                  crp_tgt_xy, dense = F)                   # # #
  gQ_uncut_crpS <- bin_expand_map(gQ_uncut_bin, binS, dims_crp,            # # #
                                  crp_tgt_xy, dense = F)                   # # #
  gU_uncut_crpS <- bin_expand_map(gU_uncut_bin, binS, dims_crp,            # # #
                                  crp_tgt_xy, dense = F)                   # # #
                                                                           # # #
  cpix <- c(crp_tgt_xy, 1, 1)                                              # # #
  cval <- c(crp_tgt_val, valz, valz2)                                      # # #
                                                                           # # #
  ghold[,,,2] <- gP_crp                                                    # # #
  ghold[,,,3] <- gIP_crp                                                   # # #
  ghold[,,,4] <- gX_crp                                                    # # #
  ghold[,,,5] <- gQ_crp                                                    # # #
  ghold[,,,6] <- gU_crp                                                    # # #
  ghold_uncut[,,,2] <- gP_uncut_crp                                        # # #
  ghold_uncut[,,,3] <- gIP_uncut_crp                                       # # #
  ghold_uncut[,,,4] <- gX_uncut_crp                                        # # #
  ghold_uncut[,,,5] <- gQ_uncut_crp                                        # # #
  ghold_uncut[,,,6] <- gU_uncut_crp                                        # # #
                                                                           # # #
  ghold_s <- ghold                                                         # # #
  ghold_uncut_s <- ghold_uncut                                             # # #
  ghold_s[,,,2] <- gP_crpS                                                 # # #
  ghold_s[,,,3] <- gIP_crpS                                                # # #
  ghold_s[,,,4] <- gX_crpS                                                 # # #
  ghold_s[,,,5] <- gQ_crpS                                                 # # #
  ghold_s[,,,6] <- gU_crpS                                                 # # #
  ghold_uncut_s[,,,2] <- gP_uncut_crpS                                     # # #
  ghold_uncut_s[,,,3] <- gIP_uncut_crpS                                    # # #
  ghold_uncut_s[,,,4] <- gX_uncut_crpS                                     # # #
  ghold_uncut_s[,,,5] <- gQ_uncut_crpS                                     # # #
  ghold_uncut_s[,,,6] <- gU_uncut_crpS                                     # # #
                                                                           # # #
  writeFITSim(ghold, paste0(out_filt_folder, g_name, "_Merged_I-p-pI-X-q-",# # #
                            "u.fits"), crpixn = cpix, crvaln = cval,       # # #
              ctypen = ctyp, header = temp_head)                           # # #
  writeFITSim(ghold_s, paste0(out_filt_folder,g_name,"_Merged_I-p-pI-X-q-",# # #
                              "u_sparse.fits"),                            # # #
              crpixn = cpix, crvaln = cval, ctypen = ctyp,                 # # #
              header = temp_head)                                          # # #
  writeFITSim(ghold_uncut, paste0(out_unfilt_folder, g_name,"_Merged_unfi",# # #
                                  "ltered_I-p-pI-X-q-u.fits"),             # # #
              crpixn = cpix, crvaln = cval, ctypen = ctyp,                 # # #
              header = temp_head)                                          # # #
  writeFITSim(ghold_uncut_s,paste0(out_unfilt_folder,g_name,"_Merged_unfi",# # #
                                   "ltered_I-p-pI-X-q-u_sparse.fits"),     # # #
              crpixn = cpix, crvaln = cval, ctypen = ctyp,                 # # #
              header = temp_head)                                          # # #
  rm(ghold_s, ghold_uncut_s)                                               # # #
  ##############################################################################
                                                                             # #
  #------------------------------------------------------------------------# # #
                                                                             # #
  ######################### CALCULATING STATISTICS ######################### # #
  ##############################################################################
  print(paste0("Calculating stats for Gal ", b, " band..."))               # # #
  Gstats["<P>(%)", b] <- median(gP_bin[,,1], na.rm = T) * 100              # # #
  Gstats["unc_<P>(%)", b] <- unc_median(gP_bin[,,1], gP_bin[,,2], 0) * 100 # # #
  Gstats["<P/unc_P>", b] <- median(gP_bin[,,1] / gP_bin[,,2], na.rm = T)   # # #
  Gstats["<X>(deg)", b] <- median(gX_bin[,,1], na.rm = T)                  # # #
  Gstats["unc_<X>(deg)", b] <- unc_median(gX_bin[,,1], gX_bin[,,2], 0)     # # #
  Gstats["<IP>", b] <- median(gIP_bin[,,1], na.rm = T)                     # # #
  Gstats["unc_<IP>", b] <- unc_median(gIP_bin[,,1], gIP_bin[,,2], 0)       # # #
  Gstats["CF", b] <- round(N_not_NA / TN_not_NA, 3)                        # # #
                                                                           # # #
  G_uncut_stats["<P>(%)", b] <- median(gP_uncut_bin[,,1], na.rm = T) * 100 # # #
  G_uncut_stats["unc_<P>(%)", b] <- unc_median(gP_uncut_bin[,,1],          # # #
                                               gP_uncut_bin[,,2], 0) * 100 # # #
  G_uncut_stats["<P/unc_P>", b] <- median(gP_uncut_bin[,,1] /              # # #
                                            gP_uncut_bin[,,2], na.rm = T)  # # #
  G_uncut_stats["<X>(deg)", b] <- median(gX_uncut_bin[,,1], na.rm = T)     # # #
  G_uncut_stats["unc_<X>(deg)", b] <- unc_median(gX_uncut_bin[,,1],        # # #
                                                 gX_uncut_bin[,,2], 0)     # # #
  G_uncut_stats["<IP>", b] <- median(gIP_uncut_bin[,,1], na.rm = T)        # # #
  G_uncut_stats["unc_<IP>", b] <- unc_median(gIP_uncut_bin[,,1],           # # #
                                             gIP_uncut_bin[,,2], 0)        # # #
  G_uncut_stats["CF", b] <- NA                                             # # #
                                                                           # # #
  print(paste0("Calculating stats for Sky ", b, " band ..."))              # # #
  Sstats["<P>(%)", b] <- median(sP_bin[,,1], na.rm = T) * 100              # # #
  Sstats["unc_<P>(%)", b] <- unc_median(sP_bin[,,1], sP_bin[,,2], 0) * 100 # # #
  Sstats["<X>(deg)", b] <- median(sX_bin[,,1], na.rm = T)                  # # #
  Sstats["unc_<X>(deg)", b] <- unc_median(sX_bin[,,1], sX_bin[,,2], 0)     # # #
  Sstats["<IP>", b] <- median(sIP_bin[,,1], na.rm = T)                     # # #
  Sstats["unc_<IP>", b] <- unc_median(sIP_bin[,,1], sIP_bin[,,2], 0)       # # #
  rm(gP_bin, sP_bin, gIP_bin, sIP_bin, gX_bin, sX_bin)                     # # #
  ##############################################################################
                                                                             # #
  #------------------------------------------------------------------------# # #
                                                                             # #
  ########################## CREATING ARROW PLOTS ########################## # #
  ##############################################################################
  print(paste0("Removing outliers that may remain within Gal maps in ", b, # # #
               " band."))                                                  # # #
                                                                           # # #
  rej_inds <- which(gP_crpS[,,1] > (median(gP_crpS[,,1], na.rm = T) +      # # #
                                      5 * mad(gP_crpS[,,1], na.rm = T)),   # # #
                    arr.ind = T)                                           # # #
  temp <- gP_crpS[,,1]                                                     # # #
  temp[rej_inds] <- NA                                                     # # #
  gP_crpS[,,1] <- temp                                                     # # #
  temp <- gIP_crpS[,,1]                                                    # # #
  temp[rej_inds] <- NA                                                     # # #
  gIP_crpS[,,1] <- temp                                                    # # #
  temp <- gX_crpS[,,1]                                                     # # #
  temp[rej_inds] <- NA                                                     # # #
  gX_crpS[,,1] <- temp                                                     # # #
  rm(temp)                                                                 # # #
                                                                           # # #
  basetitle <- paste0("Polarization Map ", acq_name, " ", b)               # # #
  basetitle_F <- paste0("Polarized flux Map ", acq_name," ", b)            # # #
  basetitle_X <- paste0("Polarization Angle ", acq_name, " ", b)           # # #
                                                                           # # #
  dif_x <- crp_tgt_xy[1] - 0:dims_crp[1]                                   # # #
  dif_y <- 0:dims_crp[2] - crp_tgt_xy[2]                                   # # #
  raseq <- tgt_coords[1] + dif_x * scl_x / cosd                            # # #
  decseq <- tgt_coords[2] + dif_y * scl_y                                  # # #
                                                                           # # #
  ra_bot <- round(raseq[round(0.1 * dims_crp[1])], 2)                      # # #
  ra_top <- round(raseq[round(0.9 * dims_crp[1])], 2)                      # # #
  dec_bot <- round(decseq[round(0.235 * dims_crp[2])], 2)                  # # #
  dec_top <- round(decseq[round(0.935 * dims_crp[2])], 2)                  # # #
                                                                           # # #
  xticks <- seq(0.1, .9, length.out = 5)                                   # # #
  yticks <- seq(0.235, .935, length.out = 5)                               # # #
  leg_X_ticks <- seq(-90, 90, 30)                                          # # #
  leg_P_ticks <- round(seq(0, max(gP_crp[,,1], na.rm = T), length.out=7),4)# # #
  if(max(gP_crp[,,1], na.rm = T) < 10){                                    # # #
    r_m <- 0.95                                                            # # #
  }else{                                                                   # # #
    r_m <- 0.9                                                             # # #
  }                                                                        # # #
  xlabs <- round(seq(ra_bot, ra_top, length.out = 5), 2)                   # # #
  ylabs <- round(seq(dec_bot, dec_top, length.out = 5), 2)                 # # #
  leg_X_labs <- paste0(leg_X_ticks, "°")                                   # # #
  leg_P_labs <- paste0(leg_P_ticks * 100, "%")                             # # #
                                                                           # # #
  rat <- dims_crp[2] / dims_crp[1]                                         # # #
                                                                           # # #
  if(length(which(!is.na(gP_crp[,,1]),arr.ind = T)) != 0){                 # # #
    print(paste0("Creating polarization angle map for Gal in ", b,         # # #
                 " band..."))                                              # # #
                                                                           # # #
    pdf(paste0(out_filt_folder, g_name, "_Pol_Angle.pdf"),                 # # #
        width = round(30/rat), height = 30)                                # # #
    par(cex.main = 6)                                                      # # #
    par(cex.axis = 4)                                                      # # #
    par(cex.lab = 5)                                                       # # #
                                                                           # # #
    # Generate a circular color scale                                      # # #
    num_colors <- 181                                                      # # #
    hue_range <- c(0, 360)                                                 # # #
    hues <- seq(hue_range[1], hue_range[2], length.out = num_colors)       # # #
    colors <- c(hcl(h = hues, c = 180, l = 65, alpha = 1), "black")        # # #
    rm(num_colors, hue_range, hues)                                        # # #
                                                                           # # #
    delt <- .8                                                             # # #
                                                                           # # #
    t_X <- gX_crp[,,1]                                                     # # #
                                                                           # # #
    if(iso_flag){                                                          # # #
      ctr_arr[which(ctr_arr == 0, arr.ind = TRUE)] <- NA                   # # #
      t_X[ctr_arr == 1] <- 91                                              # # #
    }                                                                      # # #
                                                                           # # #
    # Main Plot                                                            # # #
    image.plot(t_X, xlab = ' ', ylab = ' ', xaxt = "n", yaxt = "n",        # # #
               col = colors, bigplot = c(.09, .06 + delt, .09, .95),       # # #
               smallplot = c(.07 + delt, .08 + delt, .09, .95),            # # #
               legend.lab = "Pol. Ang. (deg)", legend.cex = 6,             # # #
               legend.line = 15, axis.args = list(at = leg_X_ticks,        # # #
                                                  labels = leg_X_labs))    # # #
                                                                           # # #
    # Main title                                                           # # #
    mtext(basetitle_X, cex = 4, line=2)                                    # # #
    # X-axis                                                               # # #
    axis(1, at = xticks, labels = paste0(xlabs, "°"), line = 3, lwd = 0)   # # #
    axis(1, at = xticks, labels = rep("", length(xticks)), lwd = 0,        # # #
         lwd.ticks = 1)                                                    # # #
    mtext("Ra (deg)", side = 1, cex = 6, line = 11)                        # # #
    # Y-axis                                                               # # #
    axis(2, at = yticks, labels = paste0(ylabs, "°"))                      # # #
    mtext("Dec (deg)", side = 2, cex = 6, line = 8)                        # # #
    # N/E referencial                                                      # # #
    arrows(x0 = delt+.1, y0 = 0.05, x1 = delt+.1, y1 = 0.1, length =0.25,  # # #
           col = "red", lwd = 5)                                           # # #
    arrows(x0 = delt+.1, y0 = 0.05, x1 = delt+.1 - 0.05, y1 = 0.05,        # # #
           length = 0.25, col = "red", lwd = 5)                            # # #
    text(x = delt+.1, y = 0.105, labels = "N", pos = 3, col = "red", cex=3)# # #
    text(x = delt+.1-0.053, y = 0.05, labels = "E", pos = 2, col = "red",  # # #
         cex = 3)                                                          # # #
    # Target Center                                                        # # #
    points(crp_tgt_xy[1] / dims_crp[1], crp_tgt_xy[2] / dims_crp[2], pch=3,# # #
           col = "black", cex = 3)                                         # # #
    rm(delt)                                                               # # #
    dev.off()                                                              # # #
                                                                           # # #
    pdf(paste0(out_filt_folder, g_name, "_Pol_Degree.pdf"),                # # #
        width = round(30/rat), height = 30)                                # # #
    par(cex.main = 6)                                                      # # #
    par(cex.axis = 4)                                                      # # #
    par(cex.lab = 5)                                                       # # #
                                                                           # # #
    # Generate a circular color scale                                      # # #
    num_colors <- 181                                                      # # #
    hue_range <- c(0, 270)                                                 # # #
    hues <- seq(hue_range[1], hue_range[2], length.out = num_colors)       # # #
    colors <- c(hcl(h = hues, c = 180, l = 65, alpha = 1), "black")        # # #
    rm(hue_range, hues)                                                    # # #
                                                                           # # #
    delt <- .8                                                             # # #
                                                                           # # #
    t_P <- gP_crp[,,1]                                                     # # #
                                                                           # # #
    if(iso_flag){                                                          # # #
      ctr_arr[which(ctr_arr == 0, arr.ind = TRUE)] <- NA                   # # #
      t_P_seq <- seq(min(t_P, na.rm = T), max(t_P, na.rm = T),             # # #
                     length.out = num_colors)                              # # #
      t_P[ctr_arr == 1] <- diff(t_P_seq[1:2]) + t_P_seq[num_colors]        # # #
    }                                                                      # # #
                                                                           # # #
    # Main Plot                                                            # # #
    image.plot(t_P, xlab = ' ', ylab = ' ', xaxt = "n", yaxt = "n",        # # #
               col = colors, bigplot = c(.1, .1 + delt * r_m, .09, .95),   # # #
               smallplot = c(.11 + delt * r_m, .12 + delt * r_m, .09, .95),# # #
               legend.lab = "Pol. Degree (%)", legend.cex = 6,             # # #
               legend.line = 15, axis.args = list(at = leg_P_ticks,        # # #
                                                  labels = leg_P_labs))    # # #
                                                                           # # #
    # Main title                                                           # # #
    mtext(basetitle, cex = 4, line=2)                                      # # #
    # X-axis                                                               # # #
    axis(1, at = xticks, labels = paste0(xlabs, "°"), line = 3, lwd = 0)   # # #
    axis(1, at = xticks, labels = rep("", length(xticks)), lwd = 0,        # # #
         lwd.ticks = 1)                                                    # # #
    mtext("Ra (deg)", side = 1, cex = 6, line = 11)                        # # #
    # Y-axis                                                               # # #
    axis(2, at = yticks, labels = paste0(ylabs, "°"))                      # # #
    mtext("Dec (deg)", side = 2, cex = 6, line = 8)                        # # #
    # N/E referencial                                                      # # #
    arrows(x0 = delt+.1, y0 = 0.05, x1 = delt+.1, y1 = 0.1, length =0.25,  # # #
           col = "red", lwd = 5)                                           # # #
    arrows(x0 = delt+.1, y0 = 0.05, x1 = delt+.1 - 0.05, y1 = 0.05,        # # #
           length = 0.25, col = "red", lwd = 5)                            # # #
    text(x = delt+.1, y = 0.105, labels = "N", pos = 3, col = "red", cex=3)# # #
    text(x = delt+.1-0.053, y = 0.05, labels = "E", pos = 2, col = "red",  # # #
         cex = 3)                                                          # # #
    # Target Center                                                        # # #
    points(crp_tgt_xy[1] / dims_crp[1], crp_tgt_xy[2] / dims_crp[2], pch=3,# # #
           col = "black", cex = 3)                                         # # #
    rm(delt)                                                               # # #
    dev.off()                                                              # # #
                                                                           # # #
    pp2p <- matrix(crp_tgt_xy / dims_crp[1:2], ncol=2)                     # # #
                                                                           # # #
    print(paste0("Creating arrow plots for Gal polarization in ", b,       # # #
                 " band..."))                                              # # #
                                                                           # # #
    create_arrow_plot_pdf(mag = gP_crpS[,,1], u_mag = gP_crpS[,,1],        # # #
                          ang = gX_crpS[,,1], u_ang = gX_crpS[,,2],        # # #
                          bkg = bkgI_crp, outpath = out_filt_folder,       # # #
                          binsize = binS,                                  # # #
                          outname = paste0(g_name, "_Merged_Arrows_PX"),   # # #
                          mtitle = paste0(basetitle), a_col = "red",       # # #
                          points = pp2p, mask = msk_crp, x_ticks = xticks, # # #
                          y_ticks = yticks, x_labs = xlabs, y_labs = ylabs,# # #
                          ref_pix = crp_tgt_xy)                            # # #
                                                                           # # #
    print(paste0("Creating arrow plots for Gal polarized flux in ", b,     # # #
                 " band..."))                                              # # #
                                                                           # # #
    create_arrow_plot_pdf(mag = gIP_crpS[,,1], u_mag = gIP_crpS[,,1],      # # #
                          ang = gX_crpS[,,1], u_ang = gX_crpS[,,2],        # # #
                          bkg = bkgI_crp, outpath = out_filt_folder,       # # #
                          binsize = binS,                                  # # #
                          outname = paste0(g_name, "_Merged_Arrows_IPX"),  # # #
                          mtitle = paste0(basetitle_F), a_col = "red",     # # #
                          points = pp2p, mask = msk_crp, x_ticks = xticks, # # #
                          y_ticks = yticks, x_labs = xlabs, y_labs = ylabs,# # #
                          ref_pix = crp_tgt_xy, stat_flag = FALSE)         # # #
                                                                           # # #
    if(iso_flag){                                                          # # #
      print(paste0("Creating arrow plots over isophotes for Gal polarizat",# # #
                   "ion in ", b, " band..."))                              # # #
                                                                           # # #
      basetitle <- paste0("Polarized flux Map ", acq_name," ", b)          # # #
                                                                           # # #
      create_arrow_iso_plot_pdf(mag = gP_crpS[,,1], u_mag = gP_crpS[,,1],  # # #
                                ang = gX_crpS[,,1], u_ang = gX_crpS[,,2],  # # #
                                bkg = ctr_arr, outpath = out_filt_folder,  # # #
                                binsize = binS,                            # # #
                                outname = paste0(g_name, "_Iso_Arrows_PX"),# # #
                                mtitle = paste0(basetitle), points = pp2p, # # #
                                mask = msk_crp, x_ticks = xticks,          # # #
                                y_ticks = yticks, x_labs = xlabs,          # # #
                                y_labs = ylabs, ref_pix = crp_tgt_xy)      # # #
                                                                           # # #
      print(paste0("Creating arrow plots over isophotes for Gal polarized",# # #
                   " flux in ", b, " band..."))                            # # #
                                                                           # # #
      basetitle <- paste0("Polarized flux Map ", acq_name," ", b)          # # #
                                                                           # # #
      create_arrow_iso_plot_pdf(mag = gIP_crpS[,,1], u_mag = gIP_crpS[,,1],# # #
                                ang = gX_crpS[,,1], u_ang = gX_crpS[,,2],  # # #
                                bkg = ctr_arr, outpath = out_filt_folder,  # # #
                                binsize = binS, stat_flag = FALSE,         # # #
                                outname = paste0(g_name,"_Iso_Arrows_IPX"),# # #
                                mtitle = paste0(basetitle_F), points =pp2p,# # #
                                mask = msk_crp, x_ticks = xticks,          # # #
                                y_ticks = yticks, x_labs = xlabs,          # # #
                                y_labs = ylabs, ref_pix = crp_tgt_xy)      # # #
    }                                                                      # # #
  }else{                                                                   # # #
    print(paste0("There are no non-NA pixels available for Gal polarizati",# # #
                 "on in ", b, " band, arrow plot will not be made..."))    # # #
  }                                                                        # # #
                                                                           # # #
  basetitle <- paste0("Unfiltered Polarization Map ", acq_name, " ", b)    # # #
  basetitle_F <- paste0("Unfiltered Polarized flux Map ", acq_name," ", b) # # #
  basetitle_X <- paste0("Unfiltered Polarization Angle ", acq_name, " ", b)# # #
                                                                           # # #
  dif_x <- crp_tgt_xy[1] - 0:dims_crp[1]                                   # # #
  dif_y <- 0:dims_crp[2] - crp_tgt_xy[2]                                   # # #
  raseq <- tgt_coords[1] + dif_x * scl_x / cosd                            # # #
  decseq <- tgt_coords[2] + dif_y * scl_y                                  # # #
                                                                           # # #
  ra_bot <- round(raseq[round(0.1 * dims_crp[1])], 2)                      # # #
  ra_top <- round(raseq[round(0.9 * dims_crp[1])], 2)                      # # #
  dec_bot <- round(decseq[round(0.235 * dims_crp[2])], 2)                  # # #
  dec_top <- round(decseq[round(0.935 * dims_crp[2])], 2)                  # # #
                                                                           # # #
  xticks <- seq(0.1, .9, length.out = 5)                                   # # #
  yticks <- seq(0.235, .935, length.out = 5)                               # # #
  leg_X_ticks <- seq(-90, 90, 30)                                          # # #
  leg_P_ticks <- round(seq(0, max(gP_uncut_crp[,,1], na.rm=T),             # # #
                           length.out = 15), 4)                            # # #
  if(max(gP_uncut_crp[,,1], na.rm = T) < 10){                              # # #
    r_m <- 0.95                                                            # # #
  }else{                                                                   # # #
    r_m <- 0.9                                                             # # #
  }                                                                        # # #
  xlabs <- round(seq(ra_bot, ra_top, length.out = 5), 2)                   # # #
  ylabs <- round(seq(dec_bot, dec_top, length.out = 5), 2)                 # # #
  leg_X_labs <- paste0(leg_X_ticks, "°")                                   # # #
  leg_P_labs <- paste0(leg_P_ticks * 100, "%")                             # # #
                                                                           # # #
  rat <- dims_crp[2] / dims_crp[1]                                         # # #
                                                                           # # #
  if(length(which(!is.na(gP_uncut_crp[,,1]), arr.ind = T)) != 0){          # # #
    print(paste0("Creating unfiltered polarization angle map for Gal in ", # # #
                 b, " band..."))                                           # # #
                                                                           # # #
    pdf(paste0(out_unfilt_folder, g_name, "_unfiltered_Pol_Angle.pdf"),    # # #
        width = round(30/rat), height = 30)                                # # #
    par(cex.main = 6)                                                      # # #
    par(cex.axis = 4)                                                      # # #
    par(cex.lab = 5)                                                       # # #
                                                                           # # #
    # Generate a circular color scale                                      # # #
    num_colors <- 181                                                      # # #
    hue_range <- c(0, 360)                                                 # # #
    hues <- seq(hue_range[1], hue_range[2], length.out = num_colors)       # # #
    colors <- c(hcl(h = hues, c = 180, l = 65, alpha = 1), "black")        # # #
    rm(num_colors, hue_range, hues)                                        # # #
                                                                           # # #
    delt <- .8                                                             # # #
                                                                           # # #
    t_X <- gX_uncut_crp[,,1]                                               # # #
                                                                           # # #
    if(iso_flag){                                                          # # #
      ctr_arr[which(ctr_arr == 0, arr.ind = TRUE)] <- NA                   # # #
      t_X[ctr_arr == 1] <- 91                                              # # #
    }                                                                      # # #
                                                                           # # #
    # Main Plot                                                            # # #
    image.plot(t_X, xlab = ' ', ylab = ' ', xaxt = "n", yaxt = "n",        # # #
               col = colors, bigplot = c(.09, .06 + delt, .09, .95),       # # #
               smallplot = c(.07 + delt, .08 + delt, .09, .95),            # # #
               legend.lab = "Pol. Ang. (deg)", legend.cex = 6,             # # #
               legend.line = 15, axis.args = list(at = leg_X_ticks,        # # #
                                                  labels = leg_X_labs))    # # #
                                                                           # # #
    # Main title                                                           # # #
    mtext(basetitle_X, cex = 4, line=2)                                    # # #
    # X-axis                                                               # # #
    axis(1, at = xticks, labels = paste0(xlabs, "°"), line = 3, lwd = 0)   # # #
    axis(1, at = xticks, labels = rep("", length(xticks)), lwd = 0,        # # #
         lwd.ticks = 1)                                                    # # #
    mtext("Ra (deg)", side = 1, cex = 6, line = 11)                        # # #
    # Y-axis                                                               # # #
    axis(2, at = yticks, labels = paste0(ylabs, "°"))                      # # #
    mtext("Dec (deg)", side = 2, cex = 6, line = 8)                        # # #
    # N/E referencial                                                      # # #
    arrows(x0 = delt+.1, y0 = 0.05, x1 = delt+.1, y1 = 0.1, length =0.25,  # # #
           col = "red", lwd = 5)                                           # # #
    arrows(x0 = delt+.1, y0 = 0.05, x1 = delt+.1 - 0.05, y1 = 0.05,        # # #
           length = 0.25, col = "red", lwd = 5)                            # # #
    text(x = delt+.1, y = 0.105, labels = "N", pos = 3, col = "red", cex=3)# # #
    text(x = delt+.1-0.053, y = 0.05, labels = "E", pos = 2, col = "red",  # # #
         cex = 3)                                                          # # #
    # Target Center                                                        # # #
    points(crp_tgt_xy[1] / dims_crp[1], crp_tgt_xy[2] / dims_crp[2], pch=3,# # #
           col = "black", cex = 3)                                         # # #
    rm(delt)                                                               # # #
    dev.off()                                                              # # #
                                                                           # # #
    pdf(paste0(out_unfilt_folder, g_name, "_unfiltered_Pol_Degree.pdf"),   # # #
        width = round(30/rat), height = 30)                                # # #
    par(cex.main = 6)                                                      # # #
    par(cex.axis = 4)                                                      # # #
    par(cex.lab = 5)                                                       # # #
                                                                           # # #
    # Generate a circular color scale                                      # # #
    num_colors <- 181                                                      # # #
    hue_range <- c(0, 270)                                                 # # #
    hues <- seq(hue_range[1], hue_range[2], length.out = num_colors)       # # #
    colors <- c(hcl(h = hues, c = 180, l = 65, alpha = 1), "black")        # # #
    rm(hue_range, hues)                                                    # # #
                                                                           # # #
    delt <- .8                                                             # # #
                                                                           # # #
    t_P <- gP_uncut_crp[,,1]                                               # # #
                                                                           # # #
    if(iso_flag){                                                          # # #
      ctr_arr[which(ctr_arr == 0, arr.ind = TRUE)] <- NA                   # # #
      t_P_seq <- seq(min(t_P, na.rm = T), max(t_P, na.rm = T),             # # #
                     length.out = num_colors)                              # # #
      t_P[ctr_arr == 1] <- diff(t_P_seq[1:2]) + t_P_seq[num_colors]        # # #
    }                                                                      # # #
                                                                           # # #
    # Main Plot                                                            # # #
    image.plot(t_P, xlab = ' ', ylab = ' ', xaxt = "n", yaxt = "n",        # # #
               col = colors, bigplot = c(.1, .1 + delt * r_m, .09, .95),   # # #
               smallplot = c(.11 + delt * r_m, .12 + delt * r_m, .09, .95),# # #
               legend.lab = "Pol. Degree (%)", legend.cex = 6,             # # #
               legend.line = 15, axis.args = list(at = leg_P_ticks,        # # #
                                                  labels = leg_P_labs))    # # #
                                                                           # # #
    # Main title                                                           # # #
    mtext(basetitle, cex = 4, line=2)                                      # # #
    # X-axis                                                               # # #
    axis(1, at = xticks, labels = paste0(xlabs, "°"), line = 3, lwd = 0)   # # #
    axis(1, at = xticks, labels = rep("", length(xticks)), lwd = 0,        # # #
         lwd.ticks = 1)                                                    # # #
    mtext("Ra (deg)", side = 1, cex = 6, line = 11)                        # # #
    # Y-axis                                                               # # #
    axis(2, at = yticks, labels = paste0(ylabs, "°"))                      # # #
    mtext("Dec (deg)", side = 2, cex = 6, line = 8)                        # # #
    # N/E referencial                                                      # # #
    arrows(x0 = delt+.1, y0 = 0.05, x1 = delt+.1, y1 = 0.1, length =0.25,  # # #
           col = "red", lwd = 5)                                           # # #
    arrows(x0 = delt+.1, y0 = 0.05, x1 = delt+.1 - 0.05, y1 = 0.05,        # # #
           length = 0.25, col = "red", lwd = 5)                            # # #
    text(x = delt+.1, y = 0.105, labels = "N", pos = 3, col = "red", cex=3)# # #
    text(x = delt+.1-0.053, y = 0.05, labels = "E", pos = 2, col = "red",  # # #
         cex = 3)                                                          # # #
    # Target Center                                                        # # #
    points(crp_tgt_xy[1] / dims_crp[1], crp_tgt_xy[2] / dims_crp[2], pch=3,# # #
           col = "black", cex = 3)                                         # # #
    rm(delt)                                                               # # #
    dev.off()                                                              # # #
                                                                           # # #
    pp2p <- matrix(crp_tgt_xy / dims_crp[1:2], ncol=2)                     # # #
                                                                           # # #
    print(paste0("Creating arrow plots for Gal unfiltered polarization in",# # #
                 " ", b, " band..."))                                      # # #
                                                                           # # #
    create_arrow_plot_pdf(mag = gP_uncut_crpS[,,1],                        # # #
                          u_mag = gP_uncut_crpS[,,1],                      # # #
                          ang = gX_uncut_crpS[,,1],                        # # #
                          u_ang = gX_uncut_crpS[,,2],                      # # #
                          bkg = bkgI_crp, outpath = out_unfilt_folder,     # # #
                          binsize = binS, ref_pix = crp_tgt_xy,            # # #
                          outname = paste0(g_name, "_Merged_Arrows_unfilt",# # #
                                           "ered_PX"),                     # # #
                          mtitle = paste0(basetitle), a_col = "red",       # # #
                          points = pp2p, mask = msk_crp, x_ticks = xticks, # # #
                          y_ticks = yticks, x_labs = xlabs, y_labs = ylabs)# # #
                                                                           # # #
    print(paste0("Creating arrow plots for Gal untilred polarized flux in",# # #
                 " ", b, " band..."))                                      # # #
                                                                           # # #
    create_arrow_plot_pdf(mag = gIP_uncut_crpS[,,1],                       # # #
                          u_mag = gIP_uncut_crpS[,,1],                     # # #
                          ang = gX_uncut_crpS[,,1],                        # # #
                          u_ang = gX_uncut_crpS[,,2],                      # # #
                          bkg = bkgI_crp, outpath = out_unfilt_folder,     # # #
                          binsize = binS, stat_flag = FALSE,               # # #
                          outname = paste0(g_name, "_Merged_Arrows_unfilt",# # #
                                           "ered_IPX"),                    # # #
                          mtitle = paste0(basetitle_F), a_col = "red",     # # #
                          points = pp2p, mask = msk_crp, x_ticks = xticks, # # #
                          y_ticks = yticks, x_labs = xlabs, y_labs = ylabs,# # #
                          ref_pix = crp_tgt_xy)                            # # #
                                                                           # # #
    if(iso_flag){                                                          # # #
      print(paste0("Creating arrow plots over isophotes for Gal unfiltere",# # #
                   "d polarization in ", b, " band..."))                   # # #
                                                                           # # #
      create_arrow_iso_plot_pdf(mag = gP_uncut_crpS[,,1],                  # # #
                                u_mag = gP_uncut_crpS[,,1],                # # #
                                ang = gX_uncut_crpS[,,1],                  # # #
                                u_ang = gX_uncut_crpS[,,2],                # # #
                                bkg = ctr_arr, outpath = out_unfilt_folder,# # #
                                binsize = binS,                            # # #
                                outname = paste0(g_name, "_Iso_Arrows_unf",# # #
                                                 "iltered_PX"),            # # #
                                mtitle = paste0(basetitle), points = pp2p, # # #
                                mask = msk_crp, x_ticks = xticks,          # # #
                                y_ticks = yticks, x_labs = xlabs,          # # #
                                y_labs = ylabs, ref_pix = crp_tgt_xy)      # # #
                                                                           # # #
      print(paste0("Creating arrow plots over isophotes for Gal unfiltere",# # #
                   "d polarized flux in ", b, " band..."))                 # # #
                                                                           # # #
      create_arrow_iso_plot_pdf(mag = gIP_uncut_crpS[,,1],                 # # #
                                u_mag = gIP_uncut_crpS[,,1],               # # #
                                ang = gX_uncut_crpS[,,1],                  # # #
                                u_ang = gX_uncut_crpS[,,2],                # # #
                                bkg = ctr_arr, outpath = out_unfilt_folder,# # #
                                binsize = binS, stat_flag = FALSE,         # # #
                                outname = paste0(g_name, "_Iso_Arrows_unf",# # #
                                                 "iltered_IPX"),           # # #
                                mtitle = paste0(basetitle_F), points =pp2p,# # #
                                mask = msk_crp, x_ticks = xticks,          # # #
                                y_ticks = yticks, x_labs = xlabs,          # # #
                                y_labs = ylabs, ref_pix = crp_tgt_xy)      # # #
    }                                                                      # # #
  }else{                                                                   # # #
    print(paste0("There are no non-NA pixels available for Gal unfiltered",# # #
                 " polarization in ", b, " band, arrow plot will not be m",# # #
                 "ade..."))                                                # # #
  }                                                                        # # #
                                                                           # # #
  if(length(which(!is.na(sP_crp[,,1]),arr.ind = T)) != 0){                 # # #
    print(paste0("Creating arrow plots for Sky polarization in ", b,       # # #
                 " band..."))                                              # # #
    basetitle <- paste0("Polarization Map ", acq_name," ", b, " Sky")      # # #
                                                                           # # #
    create_arrow_plot_pdf(mag = sP_crpS[,,1], u_mag = sP_crpS[,,1],        # # #
                          ang = sX_crpS[,,1], u_ang = sX_crpS[,,2],        # # #
                          bkg = sI_crp[,,1], outpath = out_folder,         # # #
                          outname = paste0(s_name, "_Merged_Arrows_PX"),   # # #
                          mtitle = paste0(basetitle), a_col = "green",     # # #
                          points = pp2p, mask = msk_crp, x_ticks = xticks, # # #
                          y_ticks = yticks, x_labs = xlabs, y_labs = ylabs,# # #
                          ref_pix = crp_tgt_xy)                            # # #
                                                                           # # #
    print(paste0("Creating arrow plots for Sky polarized flux in ", b,     # # #
                 " band..."))                                              # # #
    basetitle <- paste0("Polarized flux Map ", acq_name," ", b, " Sky")    # # #
                                                                           # # #
    create_arrow_plot_pdf(mag = sIP_crpS[,,1], u_mag = sIP_crpS[,,1],      # # #
                          ang = sX_crpS[,,1], u_ang = sX_crpS[,,2],        # # #
                          bkg = sI_crp[,,1], outpath = out_folder,         # # #
                          outname = paste0(s_name, "_Merged_Arrows_IPX"),  # # #
                          mtitle = paste0(basetitle_F), a_col = "green",   # # #
                          points = pp2p, mask = msk_crp, x_ticks = xticks, # # #
                          y_ticks = yticks, x_labs = xlabs, y_labs = ylabs,# # #
                          ref_pix = crp_tgt_xy)                            # # #
  }else{                                                                   # # #
    print(paste0("There are no non-NA pixels available for Sky polarizati",# # #
                 "on in ", b, " band..."))                                 # # #
  }                                                                        # # #
}                                                                          # # #
################################################################################
                                                                               #
##################### CREATING P,X,Is/Ig vs BAND PLOTS ####################### #
################################################################################
print(paste0("Creating Pg, IPg and Xg vs band plots for Gal ", acq_name,     # #
             " ..."))                                                        # #
Gpdf <- paste0(out_filt_folder, basename, "_Merged_PX_vs_Band.pdf")          # #
G_uncut_pdf <- paste0(out_unfilt_folder, basename, "_Merged_PX_vs_Band.pdf") # #
                                                                             # #
pdf(Gpdf, width = 110, height = 35)                                          # #
par(mfrow = c(1, 3))                                                         # #
par(cex = 3)                                                                 # #
                                                                             # #
Pbind <- cbind(Gstats["<P>(%)",], Sstats["<P>(%)",])                         # #
Xbind <- cbind(Gstats["<X>(deg)",], Sstats["<X>(deg)",])                     # #
IPbind <- cbind(Gstats["<IP>",], Sstats["<IP>",])                            # #
                                                                             # #
y_max <- min(c(max(c(Gstats["<P>(%)",] + Gstats["unc_<P>(%)",],              # #
                     Sstats["<P>(%)",] + Sstats["unc_<P>(%)",]), na.rm = T), # #
               100))                                                         # #
y_min <- max(c(min(c(Gstats["<P>(%)",] - Gstats["unc_<P>(%)",],              # #
                     Sstats["<P>(%)",] - Sstats["unc_<P>(%)",]),             # #
                   na.rm = T) * .2, 0))                                      # #
y_lims <- c(y_min, y_max)                                                    # #
x_lims <- c(wl_seq_nm[1], 1000)                                              # #
                                                                             # #
matplot(wl_arr[,"filt_c"], Pbind, xlab = 'Wavelength (nm)',                  # #
        ylab = 'P(%)', xlim = x_lims, ylim = y_lims, cex = 2,                # #
        main = paste0("P(%) vs band for Gal ", acq_name),                    # #
        pch = c(0,1), col = c(1,2))                                          # #
arrows(wl_arr[,"filt_c"],                                                    # #
       y0 = Gstats["<P>(%)",] - Gstats["unc_<P>(%)",],                       # #
       y1 = Gstats["<P>(%)",] + Gstats["unc_<P>(%)",],                       # #
       code = 3, length = 0.02, angle = 90, col = 1)                         # #
arrows(wl_arr[,"filt_c"],                                                    # #
       y0 = Sstats["<P>(%)",] - Sstats["unc_<P>(%)",],                       # #
       y1 = Sstats["<P>(%)",] + Sstats["unc_<P>(%)",],                       # #
       code = 3, length = 0.02, angle = 90, col = 2)                         # #
legend(850, y_max, legend = c(acq_name, "Sky"), pch = c(0,1), col = c(1,2),  # #
       cex = 2)                                                              # #
                                                                             # #
y_max <- max(c(Gstats["<X>(deg)",] + Gstats["unc_<X>(deg)",],                # #
               Sstats["<X>(deg)",] + Sstats["unc_<X>(deg)",]), na.rm = T)    # #
y_min <- min(c(Gstats["<X>(deg)",] - Gstats["unc_<X>(deg)",],                # #
               Sstats["<X>(deg)",] - Sstats["unc_<X>(deg)",]), na.rm = T)    # #
y_lims <- c(y_min, y_max)                                                    # #
                                                                             # #
matplot(wl_arr[,"filt_c"], Xbind, xlab = 'Wavelength (nm)',                  # #
        ylab = 'X(deg)', xlim = x_lims, ylim = y_lims, cex = 2,              # #
        main = paste0("X(deg) vs band for Gal ", acq_name),                  # #
        pch = c(0,1), col = c(1,2))                                          # #
arrows(wl_arr[,"filt_c"],                                                    # #
       y0 = Gstats["<X>(deg)",] - Gstats["unc_<X>(deg)",],                   # #
       y1 = Gstats["<X>(deg)",] + Gstats["unc_<X>(deg)",],                   # #
       code = 3, length = 0.02, angle = 90, col = 1)                         # #
arrows(wl_arr[,"filt_c"],                                                    # #
       y0 = Sstats["<X>(deg)",] - Sstats["unc_<X>(deg)",],                   # #
       y1 = Sstats["<X>(deg)",] + Sstats["unc_<X>(deg)",],                   # #
       code = 3, length = 0.02, angle = 90, col = 2)                         # #
legend(850, y_max, legend = c(acq_name, "Sky"), pch = c(0,1), col = c(1,2),  # #
       cex = 2)                                                              # #
                                                                             # #
y_max <- max(c(Gstats["<IP>",] + Gstats["unc_<IP>",],                        # #
               Sstats["<IP>",] + Sstats["unc_<IP>",]), na.rm = T)            # #
y_min <- min(c(Gstats["<IP>",] - Gstats["unc_<IP>",],                        # #
               Sstats["<IP>",] - Sstats["unc_<IP>",]), na.rm=T) * .2         # #
y_lims <- c(y_min, y_max)                                                    # #
                                                                             # #
matplot(wl_arr[,"filt_c"], IPbind, xlab = 'Wavelength (nm)',                 # #
        ylab = 'IP', xlim = x_lims, ylim = y_lims, cex = 2,                  # #
        main = paste0("IP vs band for Gal ", acq_name),                      # #
        pch = c(0,1), col = c(1,2))                                          # #
arrows(wl_arr[,"filt_c"],                                                    # #
       y0 = Gstats["<IP>",] - Gstats["unc_<IP>",],                           # #
       y1 = Gstats["<IP>",] + Gstats["unc_<IP>",],                           # #
       code = 3, length = 0.02, angle = 90, col = 1)                         # #
arrows(wl_arr[,"filt_c"],                                                    # #
       y0 = Sstats["<IP>",] - Sstats["unc_<IP>",],                           # #
       y1 = Sstats["<IP>",] + Sstats["unc_<IP>",],                           # #
       code = 3, length = 0.02, angle = 90, col = 2)                         # #
legend(850, y_max, legend = c(acq_name, "Sky"), pch = c(0,1), col = c(1,2),  # #
       cex = 2)                                                              # #
dev.off()                                                                    # #
                                                                             # #
pdf(G_uncut_pdf, width = 110, height = 35)                                   # #
par(mfrow = c(1, 3))                                                         # #
par(cex = 3)                                                                 # #
                                                                             # #
Pbind <- cbind(G_uncut_stats["<P>(%)",], Sstats["<P>(%)",])                  # #
Xbind <- cbind(G_uncut_stats["<X>(deg)",], Sstats["<X>(deg)",])              # #
IPbind <- cbind(G_uncut_stats["<IP>",], Sstats["<IP>",])                     # #
                                                                             # #
y_max <- min(c(max(c(G_uncut_stats["<P>(%)",] + G_uncut_stats["unc_<P>(%)",],# #
                     Sstats["<P>(%)",] + Sstats["unc_<P>(%)",]), na.rm = T), # #
               100))                                                         # #
y_min <- max(c(min(c(G_uncut_stats["<P>(%)",] - G_uncut_stats["unc_<P>(%)",],# #
                     Sstats["<P>(%)",] - Sstats["unc_<P>(%)",]),             # #
                   na.rm = T) * .2, 0))                                      # #
y_lims <- c(y_min, y_max)                                                    # #
x_lims <- c(wl_seq_nm[1], 1000)                                              # #
                                                                             # #
matplot(wl_arr[,"filt_c"], Pbind, xlab = 'Wavelength (nm)',                  # #
        ylab = 'P(%)', xlim = x_lims, ylim = y_lims, cex = 2,                # #
        main = paste0("P(%) vs band for Gal ", acq_name),                    # #
        pch = c(0,1), col = c(1,2))                                          # #
arrows(wl_arr[,"filt_c"],                                                    # #
       y0 = G_uncut_stats["<P>(%)",] - G_uncut_stats["unc_<P>(%)",],         # #
       y1 = G_uncut_stats["<P>(%)",] + G_uncut_stats["unc_<P>(%)",],         # #
       code = 3, length = 0.02, angle = 90, col = 1)                         # #
arrows(wl_arr[,"filt_c"],                                                    # #
       y0 = Sstats["<P>(%)",] - Sstats["unc_<P>(%)",],                       # #
       y1 = Sstats["<P>(%)",] + Sstats["unc_<P>(%)",],                       # #
       code = 3, length = 0.02, angle = 90, col = 2)                         # #
legend(850, y_max, legend = c(acq_name, "Sky"), pch = c(0,1), col = c(1,2),  # #
       cex = 2)                                                              # #
                                                                             # #
y_max <- max(c(G_uncut_stats["<X>(deg)",] + G_uncut_stats["unc_<X>(deg)",],  # #
               Sstats["<X>(deg)",] + Sstats["unc_<X>(deg)",]), na.rm = T)    # #
y_min <- min(c(G_uncut_stats["<X>(deg)",] - G_uncut_stats["unc_<X>(deg)",],  # #
               Sstats["<X>(deg)",] - Sstats["unc_<X>(deg)",]), na.rm = T)    # #
y_lims <- c(y_min, y_max)                                                    # #
                                                                             # #
matplot(wl_arr[,"filt_c"], Xbind, xlab = 'Wavelength (nm)',                  # #
        ylab = 'X(deg)', xlim = x_lims, ylim = y_lims, cex = 2,              # #
        main = paste0("X(deg) vs band for Gal ", acq_name),                  # #
        pch = c(0,1), col = c(1,2))                                          # #
arrows(wl_arr[,"filt_c"],                                                    # #
       y0 = G_uncut_stats["<X>(deg)",] - G_uncut_stats["unc_<X>(deg)",],     # #
       y1 = G_uncut_stats["<X>(deg)",] + G_uncut_stats["unc_<X>(deg)",],     # #
       code = 3, length = 0.02, angle = 90, col = 1)                         # #
arrows(wl_arr[,"filt_c"],                                                    # #
       y0 = Sstats["<X>(deg)",] - Sstats["unc_<X>(deg)",],                   # #
       y1 = Sstats["<X>(deg)",] + Sstats["unc_<X>(deg)",],                   # #
       code = 3, length = 0.02, angle = 90, col = 2)                         # #
legend(850, y_max, legend = c(acq_name, "Sky"), pch = c(0,1), col = c(1,2),  # #
       cex = 2)                                                              # #
                                                                             # #
y_max <- max(c(G_uncut_stats["<IP>",] + G_uncut_stats["unc_<IP>",],          # #
               Sstats["<IP>",] + Sstats["unc_<IP>",]), na.rm = T)            # #
y_min <- min(c(G_uncut_stats["<IP>",] - G_uncut_stats["unc_<IP>",],          # #
               Sstats["<IP>",] - Sstats["unc_<IP>",]), na.rm=T) * .2         # #
y_lims <- c(y_min, y_max)                                                    # #
                                                                             # #
matplot(wl_arr[,"filt_c"], IPbind, xlab = 'Wavelength (nm)',                 # #
        ylab = 'IP', xlim = x_lims, ylim = y_lims, cex = 2,                  # #
        main = paste0("IP vs band for Gal ", acq_name),                      # #
        pch = c(0,1), col = c(1,2))                                          # #
arrows(wl_arr[,"filt_c"],                                                    # #
       y0 = G_uncut_stats["<IP>",] - G_uncut_stats["unc_<IP>",],             # #
       y1 = G_uncut_stats["<IP>",] + G_uncut_stats["unc_<IP>",],             # #
       code = 3, length = 0.02, angle = 90, col = 1)                         # #
arrows(wl_arr[,"filt_c"],                                                    # #
       y0 = Sstats["<IP>",] - Sstats["unc_<IP>",],                           # #
       y1 = Sstats["<IP>",] + Sstats["unc_<IP>",],                           # #
       code = 3, length = 0.02, angle = 90, col = 2)                         # #
legend(850, y_max, legend = c(acq_name, "Sky"), pch = c(0,1), col = c(1,2),  # #
       cex = 2)                                                              # #
dev.off()                                                                    # #
################################################################################
                                                                               #
#----------------------------------------------------------------------------# #
                                                                               #
############################# CREATING STATS CSV ############################# #
################################################################################
print(paste0("Creating stats files for Gal ", acq_name, " ..."))             # #
write.table(Gstats, paste0(out_filt_folder, basename, "_Merged_Stats.csv"),  # #
            sep = ";", dec = ".")                                            # #
                                                                             # #
print(paste0("Creating stats files for unfiltered Gal ", acq_name, " ..."))  # #
write.table(G_uncut_stats, paste0(out_unfilt_folder, basename, "_unfiltered",# #
                                  "_Merged_Stats.csv"), sep = ";", dec = ".")# #
                                                                             # #
print(paste0("Creating stats files for Sky ", acq_name, " ..."))             # #
write.table(Sstats, paste0(out_folder, basename, "_Merged_Sky_Stats.csv"),   # #
            sep = ";", dec = ".")                                            # #
################################################################################
print(paste0("Cropping of target ", og_obj, " with bin size ", binS,
             " has been finished."))