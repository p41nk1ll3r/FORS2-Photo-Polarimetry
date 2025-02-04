rm(list=ls())
require(FITSio)
require(stringr)
require(jjb)
require(reticulate)

home_folder <- "/media/joaomfras/GalaxyPol/Pol-Gal/"
libs_folder <- paste0(home_folder, "Polarimetric-Imaging-Reduction-Scripts-(FO",
                      "RS2)/Commit/")
lib_head_path <- paste0(libs_folder, "fetch_fits_header_info.R")
lib_pro_path <- paste0(libs_folder, "process_lib.R")
source(lib_head_path)
source(lib_pro_path)
rm(libs_folder, lib_head_path, lib_pro_path)

################################################################################
################################# User Prompts #################################
################################################################################
obj <- "NGC-7721"                                                              #
og_folder <- paste0(home_folder, obj, "/")                                     #
print("Provide the general input path: ")                                      #
og_folder <- edit(og_folder)                                                   #
rm(home_folder, obj)                                                           #
                                                                               #
input_folder <- paste0(og_folder)                                              #
print("Provide the sources input path: ")                                      #
input_folder <- edit(input_folder)                                             #
                                                                               #
acq_path <- paste0(input_folder, "ACQ/")                                       #
print("Provide a path to move the ACQ files to: ")                             #
acq_path <- edit(acq_path)                                                     #
if(!dir.exists(acq_path)){                                                     #
  mkdir(acq_path)                                                              #
}                                                                              #
                                                                               #
mv_flag <- 1                                                                   #
                                                                               #
if(input_folder == acq_path){                                                  #
  mv_flag <- 0                                                                 #
}                                                                              #
                                                                               #
prod_path <- paste0(og_folder, "Processed_ACQ/")                               #
print("Provide an output path for data products: ")                            #
prod_path <- edit(prod_path)                                                   #
if(!dir.exists(prod_path)){                                                    #
  mkdir(prod_path)                                                             #
}                                                                              #
rm(og_folder)                                                                  #
                                                                               #
mrg_path <- paste0(prod_path, "Merged_ACQ/")                                   #
if(!dir.exists(mrg_path)){                                                     #
  mkdir(mrg_path)                                                              #
}                                                                              #
                                                                               #
cr_path <- paste0(prod_path, "ACQ_Cosmic_Rays/")                               #
if(!dir.exists(cr_path)){                                                      #
  mkdir(cr_path)                                                               #
}                                                                              #
rm(prod_path)                                                                  #
                                                                               #
cr_t <- 0.4                                                                    #
print("Define the sigma-clip value for cosmic ray detection: ")                #
cr_t <- edit(cr_t)                                                             #
################################################################################

##_##################################_###_###################################_##
#/ \--------------------------------/ \|/ \---------------------------------/ \#
#\_/--------------------------------\_/|\_/---------------------------------\_/#
################################################################################

################################################################################
############################ Setting up parameters #############################
################################################################################
#Keywords of interest in fits files' headers                                   #
EXPTIME <- "EXPTIME"                                                           #
REFX <- "CRPIX1"                                                               #
REFY <- "CRPIX2"                                                               #
VALX <- "CRVAL1"                                                               #
VALY <- "CRVAL2"                                                               #
DATE <- "DATE-OBS"                                                             #
CHIP <- "EXTNAME"                                                              #
CATG <- "HIERARCH ESO DPR CATG"                                                #
TYPE <- "HIERARCH ESO DPR TYPE"                                                #
FILTER <- "HIERARCH ESO INS FILT1 NAME"                                        #
BINX <- "HIERARCH ESO DET WIN1 BINX"                                           #
BINY <- "HIERARCH ESO DET WIN1 BINY"                                           #
SIZEX <- "HIERARCH ESO DET CHIP1 PSZX"                                         #
SIZEY <- "HIERARCH ESO DET CHIP1 PSZY"                                         #
X_GAP <- "HIERARCH ESO DET CHIP1 XGAP"                                         #
Y_GAP <- "HIERARCH ESO DET CHIP1 YGAP"                                         #
ANG_GAP <- "HIERARCH ESO DET CHIP1 RGAP"                                       #
OVERSCAN <- "HIERARCH ESO DET OUT1 OVSCY"                                      #
PRESCAN <- "HIERARCH ESO DET OUT1 PRSCY"                                       #
GAIN <- "HIERARCH ESO DET OUT1 GAIN"                                           #
RON <- "HIERARCH ESO DET OUT1 RON"                                             #
################################################################################

##_##################################_###_###################################_##
#/ \--------------------------------/ \|/ \---------------------------------/ \#
#\_/--------------------------------\_/|\_/---------------------------------\_/#
################################################################################

################################################################################
############################## Loading BIAS files ##############################
################################################################################
bias_path <- paste0(input_folder, "MASTER_BIAS/")                              #
print("Provide a path to the folder holding the Master Bias files: ")          #
bias_path <- edit(bias_path)                                                   #
                                                                               #
bias_list <- list.files(bias_path, full.names = T)                             #
rm(bias_path)                                                                  #
                                                                               #
catg_list <- get_fits_header_list_str(paths = bias_list, KEY = CATG)           #
chip_list <- get_fits_header_list_str(paths = bias_list, KEY = CHIP)           #
                                                                               #
M_chip1_ind <- intersect(which(catg_list == "MASTER", arr.ind = TRUE),         #
                         which(chip_list == "CHIP1", arr.ind = TRUE))          #
M_chip2_ind <- intersect(which(catg_list == "MASTER", arr.ind = TRUE),         #
                         which(chip_list == "CHIP2", arr.ind = TRUE))          #
                                                                               #
error_no_bias <- paste0("ERROR: the number of Master Bias files is expected t",#
                        "o 1 per CHIP. There are ", length(M_chip1_ind), " fi",#
                        "les for CHIP1 and ", length(M_chip2_ind), " for CHIP",#
                        "2. Check the path towards the files or run 'process_",#
                        "BIAS.R'.")                                            #
stopifnot(error_no_bias = length(M_chip1_ind) == 1 &&                          #
            length(M_chip2_ind) == 1)                                          #
rm(error_no_bias)                                                              #
                                                                               #
bias1_M <- readFITS(bias_list[M_chip1_ind])$imDat                              #
bias2_M <- readFITS(bias_list[M_chip2_ind])$imDat                              #
rm(bias_list)                                                                  #
################################################################################

##_##################################_###_###################################_##
#/ \--------------------------------/ \|/ \---------------------------------/ \#
#\_/--------------------------------\_/|\_/---------------------------------\_/#
################################################################################

################################################################################
############################## Loading FLAT files ##############################
################################################################################
flat_flag <- 1                                                                 #
                                                                               #
if(flat_flag == 1){                                                            #
  flat_path <- paste0(input_folder, "MASTER_FLAT/")                            #
  print("Provide a path to the folder holding the Master FLAT files: ")        #
  flat_path <- edit(flat_path)                                                 #
                                                                               #
  flat_list <- list.files(flat_path, full.names = T)                           #
  rm(flat_path)                                                                #
                                                                               #
  catg_list <- get_fits_header_list_str(paths = flat_list, KEY = CATG)         #
  chip_list <- get_fits_header_list_str(paths = flat_list, KEY = CHIP)         #
                                                                               #
  M_chip1_ind <- intersect(which(catg_list == "MASTER", arr.ind = TRUE),       #
                           which(chip_list == "CHIP1", arr.ind = TRUE))        #
  M_chip2_ind <- intersect(which(catg_list == "MASTER", arr.ind = TRUE),       #
                           which(chip_list == "CHIP2", arr.ind = TRUE))        #
                                                                               #
  error_no_flat <- paste0("ERROR: the required Master and Uncertainty Flat fi",#
                          "les were not found, please either point towards th",#
                          "e files or run 'process_FLAT.R'.")                  #
  stopifnot(error_no_flat =  length(M_chip1_ind) != 0 &&                       #
              length(M_chip1_ind) == length(M_chip2_ind))                      #
  rm(error_no_flat)                                                            #
                                                                               #
  flat_ch1_list <- flat_list[M_chip1_ind]                                      #
  flat_ch2_list <- flat_list[M_chip2_ind]                                      #
  rm(flat_list)                                                                #
                                                                               #
  flat_band_N <- length(flat_ch1_list)                                         #
  f_bands_arr <- array(NA, dim = c(flat_band_N, 2),                            #
                       dimnames = list(NULL, c("CH1", "CH2")))                 #
  rm(flat_band_N)                                                              #
                                                                               #
  print("Listing Master Flat files...")                                        #
                                                                               #
  f_bands_arr[,"CH1"] <- get_fits_header_list_str(flat_ch1_list, FILTER)       #
  f_bands_arr[,"CH2"] <- get_fits_header_list_str(flat_ch2_list, FILTER)       #
                                                                               #
  f1_u_ind <- which(f_bands_arr[,"CH1"] == "u_HIGH", arr.ind = TRUE)           #
  f1_b_ind <- which(f_bands_arr[,"CH1"] == "b_HIGH", arr.ind = TRUE)           #
  f1_v_ind <- which(f_bands_arr[,"CH1"] == "v_HIGH", arr.ind = TRUE)           #
  f1_r_ind <- which(f_bands_arr[,"CH1"] == "R_SPECIAL", arr.ind = TRUE)        #
  f1_i_ind <- which(f_bands_arr[,"CH1"] == "I_BESS", arr.ind = TRUE)           #
  f2_u_ind <- which(f_bands_arr[,"CH2"] == "u_HIGH", arr.ind = TRUE)           #
  f2_b_ind <- which(f_bands_arr[,"CH2"] == "b_HIGH", arr.ind = TRUE)           #
  f2_v_ind <- which(f_bands_arr[,"CH2"] == "v_HIGH", arr.ind = TRUE)           #
  f2_r_ind <- which(f_bands_arr[,"CH2"] == "R_SPECIAL", arr.ind = TRUE)        #
  f2_i_ind <- which(f_bands_arr[,"CH2"] == "I_BESS", arr.ind = TRUE)           #
  rm(f_bands_arr)                                                              #
}                                                                              #
rm(catg_list, chip_list, M_chip1_ind, M_chip2_ind)                             #
################################################################################

##_##################################_###_###################################_##
#/ \--------------------------------/ \|/ \---------------------------------/ \#
#\_/--------------------------------\_/|\_/---------------------------------\_/#
################################################################################

################################################################################
########### Setting up files lists, confirming bands and HWP angles ############
################################################################################
#Listing of paths to, and names of, files in the specified folder              #
pathList <- list.files(input_folder, full.names = T)                           #
fileList <- list.files(input_folder, full.names = F)                           #
rm(input_folder)                                                               #
                                                                               #
#Exclusion of non-fits files from the previous lists                           #
is_fits <- grep(".fits", pathList)                                             #
                                                                               #
if(length(is_fits) != 0){                                                      #
  pathList <- pathList[is_fits]                                                #
  fileList <- fileList[is_fits]                                                #
}                                                                              #
rm(is_fits)                                                                    #
                                                                               #
#Determining category of fits files present                                    #
#Possible types of fits files:                                                 #
#"CALIB", "ACQUISITION", "SCIENCE"                                             #
catgs_arr <-  get_fits_header_list_str(pathList, CATG)                         #
rm(TYPE, CATG)                                                                 #
                                                                               #
#Standard chip dimensions for FORS2 CCDs                                       #
dimsCH <- c(2048, 1034)                                                        #
                                                                               #
#Create list of pathList indexes for ACQUISITION fits                          #
acq_ind <- get_fits_header_key_index(catgs_arr, pathList, dimsCH,              #
                                     "ACQUISITION")                            #
acqsCount <- length(acq_ind)                                                   #
rm(catgs_arr)                                                                  #
                                                                               #
stopifnot("No ACQ files in specified folder."= (acqsCount != 0))               #
                                                                               #
print("Checking if files are paired...")                                       #
stopifnot("At least one file is not paired."= (acqsCount %% 2) == 0)           #
print("All files are paired.")                                                 #
rm(acqsCount)                                                                  #
                                                                               #
#Determining bands in which acqs are available                                 #
#This section assumes the following filters were used:                         #
#"u_HIGH", "b_HIGH", "v_HIGH", "R_SPECIAL", "I_BESS"                           #
band_str <- c("u_HIGH", "b_HIGH", "v_HIGH", "R_SPECIAL", "I_BESS")             #
band_n <- length(band_str)                                                     #
bands <- 1:band_n                                                              #
bands_arr <- get_fits_header_list_str(pathList[acq_ind], FILTER)               #
                                                                               #
u_ind <- which(bands_arr == band_str[1], arr.ind = TRUE)                       #
b_ind <- which(bands_arr == band_str[2], arr.ind = TRUE)                       #
v_ind <- which(bands_arr == band_str[3], arr.ind = TRUE)                       #
r_ind <- which(bands_arr == band_str[4], arr.ind = TRUE)                       #
i_ind <- which(bands_arr == band_str[5], arr.ind = TRUE)                       #
rm(bands_arr)                                                                  #
                                                                               #
lengths_ind <- c(length(u_ind), length(b_ind), length(v_ind), length(r_ind),   #
                 length(i_ind))                                                #
                                                                               #
for(i in 1:band_n){                                                            #
  if(lengths_ind[i] == 0){                                                     #
    bands <- bands[-which(bands == i, arr.ind = T)]                            #
                                                                               #
    if(flat_flag){                                                             #
      switch(i,                                                                #
             rm(f1_u_ind, f2_u_ind), rm(f1_b_ind, f2_b_ind),                   #
             rm(f1_v_ind, f2_v_ind), rm(f1_r_ind, f2_r_ind),                   #
             rm(f1_i_ind, f2_i_ind))                                           #
    }                                                                          #
    switch(i, rm(u_ind), rm(b_ind), rm(v_ind), rm(r_ind), rm(i_ind))           #
  }                                                                            #
}                                                                              #
rm(band_n, i, lengths_ind)                                                     #
################################################################################

##_##################################_###_###################################_##
#/ \--------------------------------/ \|/ \---------------------------------/ \#
#\_/--------------------------------\_/|\_/---------------------------------\_/#
################################################################################

################################################################################
##################### Setting up Python required libraries #####################
################################################################################
py_run_string("import lacosmic")                                               #
                                                                               #
sat <- 60000                                                                   #
cont <- 1.5                                                                    #
cr_t_str <- paste0(cr_t)                                                       #
cont_str <- paste0(cont)                                                       #
n_t <- 15                                                                      #
itN <- as.integer(4)                                                           #
py_run_string("maxN_Py = r.itN")                                               #
py_run_string("cont_Py = r.cont")                                              #
py_run_string("cr_thrsh_Py = r.cr_t")                                          #
py_run_string("n_thrsh_Py = r.n_t")                                            #
rm(cr_t, cont, n_t, itN)                                                       #
################################################################################

##_##################################_###_###################################_##
#/ \--------------------------------/ \|/ \---------------------------------/ \#
#\_/--------------------------------\_/|\_/---------------------------------\_/#
################################################################################

################################################################################
################################ ACQ PROCESSOR #################################
################################################################################
run_date <-  paste0("_", gsub(" ", "_", Sys.time()))                           #
                                                                               #
for(k in bands){                                                               #
                                                                               #
  switch(k,                                                                    #
         c(band_ind <- u_ind, rm(u_ind)),                                      #
         c(band_ind <- b_ind, rm(b_ind)),                                      #
         c(band_ind <- v_ind, rm(v_ind)),                                      #
         c(band_ind <- r_ind, rm(r_ind)),                                      #
         c(band_ind <- i_ind, rm(i_ind)))                                      #
                                                                               #
  if(flat_flag == 1){                                                          #
    switch(k,                                                                  #
           c(flat1_ind <- f1_u_ind, flat2_ind <- f2_u_ind,                     #
             rm(f1_u_ind, f2_u_ind)),                                          #
           c(flat1_ind <- f1_b_ind, flat2_ind <- f2_b_ind,                     #
             rm(f1_b_ind, f2_b_ind)),                                          #
           c(flat1_ind <- f1_v_ind, flat2_ind <- f2_v_ind,                     #
             rm(f1_v_ind, f2_v_ind)),                                          #
           c(flat1_ind <- f1_r_ind, flat2_ind <- f2_r_ind,                     #
             rm(f1_r_ind, f2_r_ind)),                                          #
           c(flat1_ind <- f1_i_ind, flat2_ind <- f2_i_ind,                     #
             rm(f1_i_ind, f2_i_ind))                                           #
    )                                                                          #
  }                                                                            #
                                                                               #
  if(length(band_ind) == 0){                                                   #
    print(paste0("No reference files in band ", band_str[k], ". Skipping."))   #
    next                                                                       #
  }                                                                            #
  if(flat_flag == 1){                                                          #
    if(length(flat1_ind) == 0 || length(flat2_ind) == 0){                      #
      print(paste0("Missing flat field files in band ", band_str[k], ". Flat ",#
                   "for Chip 1 = ", length(flat1_ind), "; Flat for Chip 2 = ", #
                   length(flat2_ind), ". Skipping."))                          #
      next                                                                     #
    }                                                                          #
                                                                               #
    flat1_M <- readFITS(flat_ch1_list[flat1_ind])$imDat                        #
    flat2_M <- readFITS(flat_ch2_list[flat2_ind])$imDat                        #
    rm(flat1_ind, flat2_ind)                                                   #
  }                                                                            #
                                                                               #
  length_ind <- length(band_ind)                                               #
                                                                               #
  if(length_ind == 0 || length_ind %% 2 != 0){                                 #
    print(paste0("Either no ACQs were found or the ones found are not pa",     #
                 "ired for the band = ", band_str[k], "..."))                  #
    next                                                                       #
  }                                                                            #
                                                                               #
  for(o in seq(1, length_ind, 2)){                                             #
                                                                               #
    print(paste0("Loading files ", o, " and ", o + 1, " of ", length_ind,      #
                 " for band ", band_str[k], "..."))                            #
                                                                               #
    #In case there is more than one set of files                               #
    #for each ang & band combination                                           #
    number <- (o + 1) / 2                                                      #
                                                                               #
    #Checking if the consecutive files share the name structure                #
    name_ch1 <- strsplit(fileList[acq_ind[band_ind[o]]],                       #
                         ".[0-9]+.fits")[[1]][1]                               #
    name_ch2 <- strsplit(fileList[acq_ind[band_ind[o + 1]]],                   #
                         ".[0-9]+.fits")[[1]][1]                               #
                                                                               #
    #Loading one chip of a pair to figure to which CCD it corresponds          #
    temp <- readFITS(pathList[acq_ind[band_ind[o]]])                           #
                                                                               #
    chip <- get_fits_header_str(temp$header, CHIP)                             #
                                                                               #
    stopifnot("ERROR: no valid chip info." = (chip == "CHIP1" ||               #
                                                chip == "CHIP2"))              #
    #Loading the files                                                         #
    if(chip == "CHIP1"){                                                       #
      ch1 <- temp                                                              #
      ch2 <- readFITS(pathList[acq_ind[band_ind[o + 1]]])                      #
      ch1_src <- fileList[acq_ind[band_ind[o]]]                                #
      ch2_src <- fileList[acq_ind[band_ind[o + 1]]]                            #
    }else{                                                                     #
      ch1 <- readFITS(pathList[acq_ind[band_ind[o + 1]]])                      #
      ch2 <- temp                                                              #
      ch1_src <- fileList[acq_ind[band_ind[o + 1]]]                            #
      ch2_src <- fileList[acq_ind[band_ind[o]]]                                #
    }                                                                          #
    rm(temp, chip)                                                             #
                                                                               #
    #Moving source files to a new folder                                       #
    if(mv_flag){                                                               #
      file.copy(pathList[acq_ind[band_ind[o]]],                                #
                paste0(acq_path, fileList[acq_ind[band_ind[o]]]))              #
      file.remove(pathList[acq_ind[band_ind[o]]])                              #
      file.copy(pathList[acq_ind[band_ind[o + 1]]],                            #
                paste0(acq_path, fileList[acq_ind[band_ind[o + 1]]]))          #
      file.remove(pathList[acq_ind[band_ind[o + 1]]])                          #
    }                                                                          #
                                                                               #
    print(paste0("ACQ set #", number, " for filter ", band_str[k],             #
                 " has been loaded: ", o + 1, "/", length_ind))                #
    print("Checking if both files match...")                                   #
                                                                               #
    #And if there has been a fluke creating or applying "lengths_ind"          #
    band_ch1 <- get_fits_header_str(ch1$header, FILTER)                        #
    band_ch2 <- get_fits_header_str(ch2$header, FILTER)                        #
                                                                               #
    stopifnot("The files are not compatible" = name_ch1 == name_ch2,           #
              band_ch1 == band_str[k],                                         #
              band_ch1 == band_ch2                                             #
    )                                                                          #
    rm(band_ch1, band_ch2, name_ch1, name_ch2)                                 #
                                                                               #
    print("Checking header info...")                                           #
                                                                               #
    ref_ch <- array(0, dim = c(1, 4),                                          #
                    dimnames = list("CH1", c("x", "y", "valx", "valy")))       #
    info_ch <- array(0, dim = c(2, 7),                                         #
                     dimnames = list(c("CH1", "CH2"),                          #
                                     c("exptime", "xSize", "ySize", "xGap",    #
                                       "yGap", "xBin", "yBin")))               #
    #Checking details                                                          #
    ref_ch["CH1", "x"] <- get_fits_header_num(ch2$header, REFX)                #
    ref_ch["CH1", "y"] <- get_fits_header_num(ch2$header, REFY)                #
    ref_ch["CH1", "valx"] <- get_fits_header_num(ch2$header, VALX)             #
    ref_ch["CH1", "valy"] <- get_fits_header_num(ch2$header, VALY)             #
                                                                               #
    info_ch["CH1", "exptime"] <- get_fits_header_num(ch1$header, EXPTIME)      #
    info_ch["CH2", "exptime"] <- get_fits_header_num(ch2$header, EXPTIME)      #
    info_ch["CH1", "xSize"] <- get_fits_header_num(ch1$header, SIZEX)          #
    info_ch["CH2", "xSize"] <- get_fits_header_num(ch2$header, SIZEX)          #
    info_ch["CH1", "ySize"] <- get_fits_header_num(ch1$header, SIZEY)          #
    info_ch["CH2", "ySize"] <- get_fits_header_num(ch2$header, SIZEY)          #
    info_ch["CH1", "xGap"] <- get_fits_header_num(ch1$header, X_GAP)           #
    info_ch["CH2", "xGap"] <- get_fits_header_num(ch2$header, X_GAP)           #
    info_ch["CH1", "yGap"] <- get_fits_header_num(ch1$header, Y_GAP)           #
    info_ch["CH2", "yGap"] <- get_fits_header_num(ch2$header, Y_GAP)           #
    info_ch["CH1", "xBin"] <- get_fits_header_num(ch1$header, BINX)            #
    info_ch["CH2", "xBin"] <- get_fits_header_num(ch2$header, BINX)            #
    info_ch["CH1", "yBin"] <- get_fits_header_num(ch1$header, BINY)            #
    info_ch["CH2", "yBin"] <- get_fits_header_num(ch2$header, BINY)            #
                                                                               #
    ang_gap <- get_fits_header_num(ch1$header, ANG_GAP) * pi / 180             #
    prescan <- get_fits_header_num(ch1$header, PRESCAN)                        #
    overscan <- get_fits_header_num(ch1$header, OVERSCAN)                      #
    valid_ch1_ys <- (prescan + 1):dimsCH[2]                                    #
    valid_ch2_ys <- 1:(dimsCH[2] - overscan)                                   #
    rm(prescan, overscan)                                                      #
                                                                               #
    #Dimensions of the pixel map of each CHIP                                  #
    ch1_dat <- ch1$imDat[,valid_ch1_ys]                                        #
    ch2_dat <- ch2$imDat[,valid_ch2_ys]                                        #
    ch1_B <- bias1_M[,valid_ch1_ys,]                                           #
    ch2_B <- bias2_M[,valid_ch2_ys,]                                           #
                                                                               #
    if(flat_flag == 1){                                                        #
      ch1_F <- flat1_M[,valid_ch1_ys,]                                         #
      ch2_F <- flat2_M[,valid_ch2_ys,]                                         #
    }                                                                          #
    rm(valid_ch1_ys, valid_ch2_ys)                                             #
                                                                               #
    dim_ch1 <- dim(ch1_dat)                                                    #
    dim_ch2 <- dim(ch2_dat)                                                    #
    N <- dim_ch2[1] * dim_ch2[2]                                               #
                                                                               #
    #Checking if all parameters checkout between CHIPs                         #
    stopifnot(info_ch["CH1", "exptime"] == info_ch["CH2", "exptime"],          #
              info_ch["CH1", "xSize"] == info_ch["CH2", "xSize"],              #
              info_ch["CH1", "ySize"] == info_ch["CH2", "ySize"],              #
              info_ch["CH1", "xGap"] == info_ch["CH2", "xGap"],                #
              info_ch["CH1", "yGap"] == info_ch["CH2", "yGap"],                #
              info_ch["CH1", "xBin"] == info_ch["CH2", "xBin"],                #
              info_ch["CH1", "yBin"] == info_ch["CH2", "yBin"])                #
                                                                               #
    #Corrections in spacing between chips                                      #
    #to take into account when merging                                         #
    xBinSize <- info_ch["CH1", "xSize"] * info_ch["CH1", "xBin"]               #
    yBinSize <- info_ch["CH1", "ySize"] * info_ch["CH1", "yBin"]               #
    xShift <- info_ch["CH1", "xGap"] / xBinSize                                #
    yShift <- info_ch["CH1", "yGap"] / yBinSize                                #
    rm(info_ch, xBinSize, yBinSize)                                            #
                                                                               #
    print("Headers' info match.")                                              #
    print("Preparing output files' header...")                                 #
                                                                               #
    #Preparing attributes for output files' header                             #
    refx_2 <- get_fits_header_num(ch2$header, REFX)                            #
    refy_2 <- get_fits_header_num(ch2$header, REFY)                            #
    refz <- 1                                                                  #
    valx_2 <- get_fits_header_num(ch2$header, VALX)                            #
    valy_2 <- get_fits_header_num(ch2$header, VALY)                            #
    valz <- NA                                                                 #
    typex_2 <- "  'PIXEL     '              / Coordinate system of x-axis "    #
    typey_2 <- "  'PIXEL     '              / Coordinate system of y-axis "    #
    typez_2 <- "  'DATA & UNC'              / Coordinate system of z-axis "    #
                                                                               #
    crpix <- c(refx_2, refy_2, refz)                                           #
    crval <- c(valx_2, valy_2, valz)                                           #
    ctype <- c(typex_2, typey_2, typez_2)                                      #
    rm(refx_2, refy_2, refz, valx_2, valy_2, valz, typex_2, typey_2, typez_2)  #
                                                                               #
    date_name <- get_fits_header_str(ch1$header, DATE)                         #
    fits_name <- paste0("FORS2_", date_name, "_", band_str[k], "_#", number,   #
                        run_date)                                              #
    rm(date_name)                                                              #
                                                                               #
    #Preparing merged chips ACQ files' header                                  #
    temp_head <- ch2$header                                                    #
    temp_head <- modVal('EXTNAME', 'MERGED', "Extension name", temp_head)      #
    temp_head <- addKwv('ORIGIN_1', ch1_src, "CHIP1 source file", temp_head)   #
    temp_head <- addKwv('ORIGIN_2', ch2_src, "CHIP2 source file", temp_head)   #
    temp_head <- addKwv('CAL-BIAS', 'YES', "Was BIAS subtracted?", temp_head)  #
    rm(ch1_src, ch2_src)                                                       #
                                                                               #
    if(flat_flag == 1){                                                        #
      temp_head <- addKwv('CAL-FLAT', 'YES', "Was FLAT divided?", temp_head)   #
    }else{                                                                     #
      temp_head <- addKwv('CAL-FLAT', 'NO', "Was FLAT divided?", temp_head)    #
    }                                                                          #
                                                                               #
    temp_head <- addKwv('CAL-COSM', 'YES', "Were cosmic-rays (CR) removed?",   #
                        temp_head)                                             #
    temp_head <- addKwv('COSM-SIG', cr_t_str, "Sigma used for CR removal",     #
                        temp_head)                                             #
    temp_head <- addKwv('COSM-CTR', cont_str, "Contrast used for CR removal",  #
                        temp_head)                                             #
    temp_head <- addKwv('FILETYP1', 'ACQUISITION', "File type", temp_head)     #
    temp_head <- addKwv('FILETYP2', 'MIXED BEAMS', "File sub-type", temp_head) #
    temp_head <- addKwv('SOURCE', 'process_ACQ.R',                             #
                        "Script used to generate this file", temp_head)        #
    temp_head <- addKwv('AUTHOR', 'JRS2023', "Script creator", temp_head)      #
                                                                               #
    ### Compensate for translational and rotational gaps and merge chips ##### #
    ########################################################################## #
    print("Preparing chips for merger...")                                   # #
    rot_info <- array(0, dim = c(N, 11),                                     # #
                      dimnames = list(NULL, c("x1", "y1", "x2", "y2", "r",   # #
                                              "ang", "ang'", "x1'", "y1'",   # #
                                              "x2'", "y2'")))                # #
    x_c <- dim_ch2[1] / 2                                                    # #
    y_c <- dim_ch2[2] / 2                                                    # #
                                                                             # #
    print("Checking coordinates...")                                         # #
    #Entering present pixel cartesian coordinates                            # #
    rot_info[,"x2"] <- rep(1:dim_ch2[1], times = dim_ch2[2])                 # #
    rot_info[,"y2"] <- rep(1:dim_ch2[2], each = dim_ch2[1])                  # #
    rot_info[,"x1"] <- rep(1:dim_ch1[1], times = dim_ch1[2])                 # #
    rot_info[,"y1"] <- rep(1:dim_ch1[2], each = dim_ch1[1])                  # #
    rm(dim_ch1)                                                              # #
                                                                             # #
    #Calculating pixel polar coordinates                                     # #
    rot_info[,"r"] <- sqrt((rot_info[,"x2"] - x_c)^2 +                       # #
                             (rot_info[,"y2"] - y_c)^2)                      # #
                                                                             # #
    rot_info[,"ang"] <- atan2(rot_info[,"y2"] - y_c, rot_info[,"x2"] - x_c)  # #
                                                                             # #
    #Correcting angle for the angular gap between CHIP1 and CHIP2            # #
    rot_info[,"ang'"] <- rot_info[,"ang"] + ang_gap                          # #
    rm(ang_gap)                                                              # #
                                                                             # #
    print("Calculating new coordinates...")                                  # #
    #Calculating new pixel cartesian coordinates                             # #
    rot_info[,"x2'"] <- rot_info[,"r"] * cos(rot_info[,"ang'"]) + x_c        # #
    rot_info[,"y2'"] <- rot_info[,"r"] * sin(rot_info[,"ang'"]) + y_c        # #
    rm(x_c, y_c)                                                             # #
                                                                             # #
    #Calculating pixel coordinate adjustment to fit with merger              # #
    rot_min <- c(min(floor(rot_info[,"x2'"]), na.rm = TRUE),                 # #
                 min(floor(rot_info[,"y2'"]), na.rm = TRUE))                 # #
    move_rot <- 1 - rot_min                                                  # #
    rm(rot_min)                                                              # #
                                                                             # #
    print("Compensating for off-grid coordinates...")                        # #
    #Adjusting pixels coordinates so that they're all positive               # #
    #(The rotation may have transported some pixels to other quadrants)      # #
    rot_info[,"x2'"] <- rot_info[,"x2'"] + move_rot[1]                       # #
    rot_info[,"y2'"] <- rot_info[,"y2'"] + move_rot[2]                       # #
    rot_info[,"x1'"] <- rot_info[,"x1"] + move_rot[1] + xShift               # #
    rot_info[,"y1'"] <- rot_info[,"y1"] + move_rot[2] + dim_ch2[2] + yShift  # #
    rm(dim_ch2, yShift)                                                      # #
                                                                             # #
    mrg_dim <- c(max(c(ceiling(rot_info[,"x2'"]),                            # #
                       ceiling(rot_info[,"x1'"])), na.rm = TRUE),            # #
                 max(ceiling(rot_info[,"y1'"]), na.rm = TRUE))               # #
    mrg_data <- array(0, dim = c(mrg_dim, 2),                                # #
                      dimnames = list(NULL, NULL, c("Data", "Unc")))         # #
    ref_mrg <- c(ref_ch["CH1", "x"] + move_rot[1] + xShift,                  # #
                 ref_ch["CH1", "y"] + move_rot[2])                           # #
    rm(ref_ch, move_rot, xShift)                                             # #
                                                                             # #
    print("Calculating uncertainties...")                                    # #
    #De-bias and load chips onto a temporary holder                          # #
    temp_ch1 <- array(0, dim = dim(ch1_dat))                                 # #
    temp_ch2 <- array(0, dim = dim(ch2_dat))                                 # #
                                                                             # #
    gain_1 <- get_fits_header_num(ch1$header, GAIN)                          # #
    ron_1 <- get_fits_header_num(ch1$header, RON)                            # #
    gain_2 <- get_fits_header_num(ch2$header, GAIN)                          # #
    ron_2 <- get_fits_header_num(ch2$header, RON)                            # #
    rm(ch1, ch2)                                                             # #
                                                                             # #
    temp_ch1 <- ch1_dat - ch1_B[,,1]                                         # #
    temp_ch2 <- ch2_dat - ch2_B[,,1]                                         # #
    temp_U1 <- unc_add(sqrt(ch1_dat * gain_1 + (ron_1 * gain_1)^2),          # #
                       ch1_B[,,2])                                           # #
    temp_U2 <- unc_add(sqrt(ch2_dat * gain_2 + (ron_2 * gain_2)^2),          # #
                       ch2_B[,,2])                                           # #
    tch1_neg <- which(temp_ch1 <= 0, arr.ind = T)                            # #
    tch2_neg <- which(temp_ch2 <= 0, arr.ind = T)                            # #
    temp_ch1[tch1_neg] <- NA                                                 # #
    temp_U1[tch1_neg] <- NA                                                  # #
    temp_ch2[tch2_neg] <- NA                                                 # #
    temp_U2[tch2_neg] <- NA                                                  # #
    rm(ch1_B, ch2_B, ch1_dat, ch2_dat, tch1_neg, tch2_neg, gain_1, gain_2,   # #
       ron_1, ron_2)                                                         # #
                                                                             # #
    if(flat_flag == 1){                                                      # #
      t_ch1 <- temp_ch1 / ch1_F[,,1]                                         # #
      t_ch2 <- temp_ch2 / ch2_F[,,1]                                         # #
      temp_U1 <- unc_div(temp_ch1, ch1_F[,,1], temp_U1, ch1_F[,,2], t_ch1)   # #
      temp_U2 <- unc_div(temp_ch2, ch2_F[,,1], temp_U2, ch2_F[,,2], t_ch2)   # #
      temp_ch1 <- t_ch1                                                      # #
      temp_ch2 <- t_ch2                                                      # #
      rm(ch1_F, ch2_F, t_ch1, t_ch2)                                         # #
    }                                                                        # #
                                                                             # #
    print("Merging CHIPs...")                                                # #
    #Loading CHIPS data into new array                                       # #
    for(n in 1:N){                                                           # #
      mrg_data[rot_info[n, "x1'"], rot_info[n, "y1'"], "Data"] <-            # #
        temp_ch1[rot_info[n, "x1"], rot_info[n, "y1"]]                       # #
      mrg_data[rot_info[n, "x1'"], rot_info[n, "y1'"], "Unc"] <-             # #
        temp_U1[rot_info[n, "x1"], rot_info[n, "y1"]]                        # #
      mrg_data[rot_info[n, "x2'"], rot_info[n, "y2'"], "Data"] <-            # #
        temp_ch2[rot_info[n, "x2"], rot_info[n, "y2"]]                       # #
      mrg_data[rot_info[n, "x2'"], rot_info[n, "y2'"], "Unc"] <-             # #
        temp_U2[rot_info[n, "x2"], rot_info[n, "y2"]]                        # #
    }                                                                        # #
                                                                             # #
    mrg_data[which(is.na(mrg_data), arr.ind = TRUE)] <- 0                    # #
    rm(temp_ch1, temp_ch2, temp_U1, temp_U2, n, rot_info)                    # #
                                                                             # #
    #Must update crpix after aligning chips since chip 2 had to be moved     # #
    crpix[1:2] <- ref_mrg[1:2]                                               # #
    rm(ref_mrg)                                                              # #
                                                                             # #
    #Creating Merged Acquisition file                                        # #
    writeFITSim(mrg_data, file = paste0(mrg_path, fits_name, "-Merged-Acq",  # #
                                        ".fits"), crpixn = crpix,            # #
                crvaln = crval, ctypen = ctype, header = temp_head)          # #
                                                                             # #
    print(paste0("Merged chips, debiased ACQ file number ", number, " ",     # #
                 "for filter ", band_str[k], " has been created."))          # #
    ############################################################################
                                                                               #
    #------------------------------------------------------------------------# #
                                                                               #
    ############################# Removing Cosmic Rays ####################### #
    ########################################################################## #
    print("Starting cosmic ray removal preprocessing for lacosmic pkg...")   # #
                                                                             # #
    #Masking NAs, Saturated pixels and pixels with flux <= 0                 # #
    print("Identifying and masking NA, saturated and negative pixels...")    # #
    na_ind <- which(is.na(mrg_data[,,"Data"]), arr.ind = T)                  # #
    sat_ind <- which(mrg_data[,,"Data"] == sat, arr.ind = T)                 # #
    neg_ind <- which(mrg_data[,,"Data"] <= 0, arr.ind = T)                   # #
                                                                             # #
    mask <- array(FALSE, dim = mrg_dim)                                      # #
    mask[na_ind] <- TRUE                                                     # #
    mask[sat_ind] <- TRUE                                                    # #
    mask[neg_ind] <- TRUE                                                    # #
    rm(na_ind, sat_ind, neg_ind, mrg_dim)                                    # #
                                                                             # #
    #Loading data to Python environment                                      # #
    print("Loading masked data to Python env...")                            # #
    temp <- mrg_data[,,"Data"]                                               # #
    py_run_string("data_Py = r.temp")                                        # #
    py_run_string("mask_Py = r.mask")                                        # #
    temp <- mrg_data[,,"Unc"]                                                # #
    py_run_string("data_Err_Py = r.temp")                                    # #
    rm(mask, temp)                                                           # #
                                                                             # #
    #Setting up lacosmic param. and loading to Python environment            # #
    print("Setting up lacosmic param. and loading them to Python env...")    # #
                                                                             # #
    lacosm_cmd <- paste0("noCR_data_Py, cr_map = lacosmic.lacosmic(data =",  # #
                         " data_Py, contrast = cont_Py, cr_threshold = cr",  # #
                         "_thrsh_Py, neighbor_threshold = n_thrsh_Py, err",  # #
                         "or = data_Err_Py, mask = mask_Py, maxiter = max",  # #
                         "N_Py, border_mode = 'constant')")                  # #
                                                                             # #
    #Running lacosmic                                                        # #
    print("Running lacosmic...")                                             # #
    py_run_string(lacosm_cmd)                                                # #
    mrg_data[,,"Data"] <- py$noCR_data_Py                                    # #
    rm(lacosm_cmd)                                                           # #
                                                                             # #
    #Creating CR removed Acquisition file                                    # #
    writeFITSim(mrg_data, file = paste0(mrg_path, fits_name,                 # #
                                        "-Merged-Acq-no-CRs.fits"),          # #
                crpixn = crpix, crvaln = crval, ctypen = ctype,              # #
                header = temp_head)                                          # #
    writeFITSim(py$cr_map, paste0(cr_path, fits_name, "-Merged-Acq-CR.fits"))# #
    rm(mrg_data, temp_head)                                                  # #
                                                                             # #
    print(paste0("Merged chips, debiased, Cosmic Ray removed Acquisition fi",# #
                 "le number ", number, " for filter ", band_str[k],          # #
                 " has been created."))                                      # #
    ############################################################################
                                                                               #
  }                                                                            #
  rm(o, length_ind, band_ind, crpix, crval, ctype, fits_name, number)          #
}                                                                              #
rm(acq_path, mrg_path, cr_path, k, EXPTIME, REFX, REFY, VALX, CHIP, VALY, BINX,#
   BINY, SIZEX, SIZEY, DATE, X_GAP, Y_GAP, FILTER, ANG_GAP, sat, PRESCAN, N,   #
   OVERSCAN, band_str, pathList, cr_t_str, mv_flag, bias1_M, bias2_M, dimsCH,  #
   cont_str, acq_ind, bands, fileList, run_date, RON, GAIN)                    #
                                                                               #
if(flat_flag == 1){                                                            #
  rm(flat1_M, flat2_M,flat_ch1_list, flat_ch2_list)                            #
}                                                                              #
rm(flat_flag)                                                                  #
################################################################################
################################################################################