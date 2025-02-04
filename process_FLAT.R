rm(list=ls())
require(FITSio)
require(stringr)
require(jjb)

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
input_folder <- paste0(home_folder, obj, "/")                                  #
print("Provide an input path: ")                                               #
input_folder <- edit(input_folder)                                             #
rm(home_folder, obj)                                                           #
                                                                               #
flat_path <- paste0(input_folder, "FLAT/")                                     #
print("Provide a path to move the FLAT files to: ")                            #
flat_path <- edit(flat_path)                                                   #
                                                                               #
if(!dir.exists(flat_path)){                                                    #
  mkdir(flat_path)                                                             #
}                                                                              #
                                                                               #
mv_flag <- 1                                                                   #
                                                                               #
if(input_folder == flat_path){                                                 #
  mv_flag <- 0                                                                 #
}                                                                              #
                                                                               #
Mflat_path <- paste0(input_folder, "MASTER_FLAT/")                             #
print("Provide an output path: ")                                              #
Mflat_path <- edit(Mflat_path)                                                 #
                                                                               #
if(!dir.exists(Mflat_path)){                                                   #
  mkdir(Mflat_path)                                                            #
}                                                                              #
################################################################################

##_##################################_###_###################################_##
#/ \--------------------------------/ \|/ \---------------------------------/ \#
#\_/--------------------------------\_/|\_/---------------------------------\_/#
################################################################################

################################################################################
############################ Setting up parameters #############################
################################################################################
#Keywords of interest in fits files' headers                                   #
REFX <- "CRPIX1"                                                               #
REFY <- "CRPIX2"                                                               #
VALX <- "CRVAL1"                                                               #
VALY <- "CRVAL2"                                                               #
DATE <- "DATE-OBS"                                                             #
CHIP <- "EXTNAME"                                                              #
CATG <- "HIERARCH ESO DPR CATG"                                                #
TYPE <- "HIERARCH ESO DPR TYPE"                                                #
FILTER <- "HIERARCH ESO INS FILT1 NAME"                                        #
GAIN <- "HIERARCH ESO DET OUT1 GAIN"                                           #
RON <- "HIERARCH ESO DET OUT1 RON"                                             #
################################################################################

##_##################################_###_###################################_##
#/ \--------------------------------/ \|/ \---------------------------------/ \#
#\_/--------------------------------\_/|\_/---------------------------------\_/#
################################################################################

################################################################################
#################### Loading the required calibration files ####################
################################################################################
bias_path <- paste0(input_folder, "MASTER_BIAS/")                              #
print("Provide a path to the folder holding the Master Bias files: ")          #
bias_path <- edit(bias_path)                                                   #
                                                                               #
bias_list <- list.files(bias_path, full.names = T)                             #
                                                                               #
catg_list <- get_fits_header_list_str(paths = bias_list, KEY = CATG)           #
chip_list <- get_fits_header_list_str(paths = bias_list, KEY = CHIP)           #
                                                                               #
M_chip1_ind <- intersect(which(catg_list == "MASTER", arr.ind = TRUE),         #
                         which(chip_list == "CHIP1", arr.ind = TRUE))          #
M_chip2_ind <- intersect(which(catg_list == "MASTER", arr.ind = TRUE),         #
                         which(chip_list == "CHIP2", arr.ind = TRUE))          #
rm(catg_list, chip_list)                                                       #
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
bias1_M <- readFITS(bias_list[M_chip1_ind])                                    #
bias2_M <- readFITS(bias_list[M_chip2_ind])                                    #
rm(bias_list, M_chip1_ind, M_chip2_ind)                                        #
################################################################################

##_##################################_###_###################################_##
#/ \--------------------------------/ \|/ \---------------------------------/ \#
#\_/--------------------------------\_/|\_/---------------------------------\_/#
################################################################################

################################################################################
############################# Setting up files list ############################
################################################################################
pathList <- list.files(input_folder, full.names = TRUE)                        #
fileList <- list.files(input_folder, full.names = FALSE)                       #
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
#Determining type of fits files present                                        #
#Possible types of fits files: "BIAS", "DARK", "FLAT,SKY", "OBJECT", "STD"     #
types_arr <- get_fits_header_list_str(pathList, TYPE)                          #
rm(TYPE)                                                                       #
                                                                               #
#Standard chip dimensions for FORS2 CCDs                                       #
dimsCH <- c(2048, 1034)                                                        #
################################################################################

##_##################################_###_###################################_##
#/ \--------------------------------/ \|/ \---------------------------------/ \#
#\_/--------------------------------\_/|\_/---------------------------------\_/#
################################################################################

################################################################################
############ Find Flats and Create List of Indexes for each Filter #############
################################################################################
flat_ind <- get_fits_header_key_index(types_arr, pathList, dimsCH, "FLAT,SKY") #
rm(types_arr)                                                                  #
                                                                               #
#Determining bands in which flats are available                                #
#This section assumes the following filters were used:                         #
#"u_HIGH", "b_HIGH", "v_HIGH", "R_SPECIAL", "I_BESS"                           #
bands_str <- c("u_HIGH", "b_HIGH", "v_HIGH", "R_SPECIAL", "I_BESS")            #
band_n <- length(bands_str)                                                    #
bands <- 1:band_n                                                              #
bands_array <- get_fits_header_list_str(pathList[flat_ind], FILTER)            #
rm(FILTER)                                                                     #
                                                                               #
u_ind <- which(bands_array == bands_str[1], arr.ind = TRUE)                    #
b_ind <- which(bands_array == bands_str[2], arr.ind = TRUE)                    #
v_ind <- which(bands_array == bands_str[3], arr.ind = TRUE)                    #
r_ind <- which(bands_array == bands_str[4], arr.ind = TRUE)                    #
i_ind <- which(bands_array == bands_str[5], arr.ind = TRUE)                    #
rm(bands_array)                                                                #
                                                                               #
lengths_ind <- c(length(u_ind), length(b_ind), length(v_ind), length(r_ind),   #
                 length(i_ind))                                                #
                                                                               #
loadFlat_u_ch1 <- array(NA, dim = c(dimsCH, lengths_ind[1] / 2))               #
loadFlat_u_ch2 <- array(NA, dim = c(dimsCH, lengths_ind[1] / 2))               #
loadFlat_b_ch1 <- array(NA, dim = c(dimsCH, lengths_ind[2] / 2))               #
loadFlat_b_ch2 <- array(NA, dim = c(dimsCH, lengths_ind[2] / 2))               #
loadFlat_v_ch1 <- array(NA, dim = c(dimsCH, lengths_ind[3] / 2))               #
loadFlat_v_ch2 <- array(NA, dim = c(dimsCH, lengths_ind[3] / 2))               #
loadFlat_r_ch1 <- array(NA, dim = c(dimsCH, lengths_ind[4] / 2))               #
loadFlat_r_ch2 <- array(NA, dim = c(dimsCH, lengths_ind[4] / 2))               #
loadFlat_i_ch1 <- array(NA, dim = c(dimsCH, lengths_ind[5] / 2))               #
loadFlat_i_ch2 <- array(NA, dim = c(dimsCH, lengths_ind[5] / 2))               #
                                                                               #
for(i in 1:band_n){                                                            #
  if(lengths_ind[i] == 0){                                                     #
    bands <- bands[-which(bands == i, arr.ind = T)]                            #
                                                                               #
    switch(i,                                                                  #
           rm(u_ind, loadFlat_u_ch1, loadFlat_u_ch2),                          #
           rm(b_ind, loadFlat_b_ch1, loadFlat_b_ch2),                          #
           rm(v_ind, loadFlat_v_ch1, loadFlat_v_ch2),                          #
           rm(r_ind, loadFlat_r_ch1, loadFlat_r_ch2),                          #
           rm(i_ind, loadFlat_i_ch1, loadFlat_i_ch2))                          #
    print(paste0("No flats were found or the ones found are not paired for fi",#
                 "lter: ", bands_str[i]))                                      #
  }                                                                            #
}                                                                              #
rm(band_n, i)                                                                  #
################################################################################

##_##################################_###_###################################_##
#/ \--------------------------------/ \|/ \---------------------------------/ \#
#\_/--------------------------------\_/|\_/---------------------------------\_/#
################################################################################

################################################################################
######### Creating Master FLATs for CHIP 1 and CHIP 2 for each Filter ##########
############## (this step assumes no FLAT has been pre-processed) ##############
################################################################################
for(k in bands){                                                               #
                                                                               #
  switch(k,                                                                    #
    c(t_ld1 <- loadFlat_u_ch1, t_ld2 <- loadFlat_u_ch2, temp_ind <- u_ind,     #
      temp_U1 <- loadFlat_u_ch1, temp_U2 <- loadFlat_u_ch2,                    #
      rm(loadFlat_u_ch1, loadFlat_u_ch2, u_ind)),                              #
    c(t_ld1 <- loadFlat_b_ch1, t_ld2 <- loadFlat_b_ch2, temp_ind <- b_ind,     #
      temp_U1 <- loadFlat_b_ch1, temp_U2 <- loadFlat_b_ch2,                    #
      rm(loadFlat_b_ch1, loadFlat_b_ch2, b_ind)),                              #
    c(t_ld1 <- loadFlat_v_ch1, t_ld2 <- loadFlat_v_ch2, temp_ind <- v_ind,     #
      temp_U1 <- loadFlat_v_ch1, temp_U2 <- loadFlat_v_ch2,                    #
      rm(loadFlat_v_ch1, loadFlat_v_ch2, v_ind)),                              #
    c(t_ld1 <- loadFlat_r_ch1, t_ld2 <- loadFlat_r_ch2, temp_ind <- r_ind,     #
      temp_U1 <- loadFlat_r_ch1, temp_U2 <- loadFlat_r_ch2,                    #
      rm(loadFlat_r_ch1, loadFlat_r_ch2, r_ind)),                              #
    c(t_ld1 <- loadFlat_i_ch1, t_ld2 <- loadFlat_i_ch2, temp_ind <- i_ind,     #
      temp_U1 <- loadFlat_i_ch1, temp_U2 <- loadFlat_i_ch2,                    #
      rm(loadFlat_i_ch1, loadFlat_i_ch2, i_ind))                               #
  )                                                                            #
                                                                               #
  for(f in seq(1, lengths_ind[k], 2)){                                         #
                                                                               #
    fi <- (f + 1) / 2                                                          #
                                                                               #
    ############################ Loading the files ########################### #
    ############################################################################
    tempA <- readFITS(pathList[flat_ind[temp_ind[f]]])                       # #
    tempB <- readFITS(pathList[flat_ind[temp_ind[f + 1]]])                   # #
                                                                             # #
    chipA <- get_fits_header_str(tempA$header, CHIP)                         # #
                                                                             # #
    if(chipA == "CHIP1"){                                                    # #
      ch1 <- tempA                                                           # #
      ch2 <- tempB                                                           # #
    }else{                                                                   # #
      ch1 <- tempB                                                           # #
      ch2 <- tempA                                                           # #
    }                                                                        # #
    rm(tempA, tempB, chipA)                                                  # #
                                                                             # #
    #Cutting Margins                                                         # #
    ch1$imDat[,964:dimsCH[2]] <- NA                                          # #
    ch1$imDat[,1:7] <- NA                                                    # #
    ch1$imDat[1:184,] <- NA                                                  # #
    ch1$imDat[1860:dimsCH[1],] <- NA                                         # #
    ch2$imDat[,1028:dimsCH[2]] <- NA                                         # #
    ch2$imDat[,1:317] <- NA                                                  # #
    ch2$imDat[1:192,] <- NA                                                  # #
    ch2$imDat[1862:dimsCH[1],] <- NA                                         # #
                                                                             # #
    #If assumptions are correct, we can load each chip into the corresponding# #
    #loader array after subtracting the BIAS                                 # #
    gain_1 <- get_fits_header_num(ch1$header, GAIN)                          # #
    ron_1 <- get_fits_header_num(ch1$header, RON)                            # #
    gain_2 <- get_fits_header_num(ch2$header, GAIN)                          # #
    ron_2 <- get_fits_header_num(ch2$header, RON)                            # #
                                                                             # #
    u_ch1 <- sqrt(ch1$imDat * gain_1 + (ron_1 * gain_1)^2)                   # #
    u_ch2 <- sqrt(ch2$imDat * gain_2 + (ron_2 * gain_2)^2)                   # #
    rm(gain_1, gain_2, ron_1, ron_2)                                         # #
                                                                             # #
    ch1$imDat <- ch1$imDat - bias1_M$imDat[,,1]                              # #
    ch2$imDat <- ch2$imDat - bias2_M$imDat[,,1]                              # #
                                                                             # #
    uf_ch1 <- unc_add(u_ch1, bias1_M$imDat[,,2])                             # #
    uf_ch2 <- unc_add(u_ch2, bias2_M$imDat[,,2])                             # #
                                                                             # #
    medF1 <- median(ch1$imDat, na.rm = TRUE)                                 # #
    medFU1 <- unc_median(ch1$imDat, uf_ch1, 0)                               # #
    medF2 <- median(ch2$imDat, na.rm = TRUE)                                 # #
    medFU2 <- unc_median(ch2$imDat, uf_ch2, 0)                               # #
                                                                             # #
    t_ld1[,,fi] <- ch1$imDat / medF1                                         # #
    t_ld2[,,fi] <- ch2$imDat / medF2                                         # #
    temp_U1[,,fi] <- unc_div(ch1$imDat, medF1, uf_ch1, medFU2, t_ld1[,,fi])  # #
    temp_U2[,,fi] <- unc_div(ch2$imDat, medF2, uf_ch2, medFU2, t_ld2[,,fi])  # #
    rm(u_ch1, uf_ch1, medF1, medFU1, u_ch2, uf_ch2, medF2, medFU2)           # #
                                                                             # #
    #Handling header information for master flat files                       # #
    if(f == (lengths_ind[k] - 1)){                                           # #
      refx_1 <- get_fits_header_num(ch1$header, REFX)                        # #
      refy_1 <- get_fits_header_num(ch1$header, REFY)                        # #
      refx_2 <- get_fits_header_num(ch2$header, REFX)                        # #
      refy_2 <- get_fits_header_num(ch2$header, REFY)                        # #
      refz <- 1                                                              # #
      crpix1 <- c(refx_1, refy_1, refz)                                      # #
      crpix2 <- c(refx_2, refy_2, refz)                                      # #
      rm(refx_1, refx_2, refy_1, refy_2, refz)                               # #
                                                                             # #
      valx_1 <- get_fits_header_num(ch1$header, VALX)                        # #
      valy_1 <- get_fits_header_num(ch1$header, VALY)                        # #
      valx_2 <- get_fits_header_num(ch2$header, VALX)                        # #
      valy_2 <- get_fits_header_num(ch2$header, VALY)                        # #
      valz <- NA                                                             # #
      crval1 <- c(valx_1, valy_1, valz)                                      # #
      crval2 <- c(valx_2, valy_2, valz)                                      # #
      rm(valy_1, valy_2, valx_1, valx_2, valz)                               # #
                                                                             # #
      typex <- "  'PIXEL     '              / Coordinate system of x-axis "  # #
      typey <- "  'PIXEL     '              / Coordinate system of y-axis "  # #
      typez <- "  'DATA & UNC'              / Coordinate system of z-axis "  # #
      ctype <- c(typex, typey, typez)                                        # #
      rm(typex, typey, typez)                                                # #
                                                                             # #
      catg <- get_fits_header_str(ch1$header, CATG)                          # #
                                                                             # #
      catg_repl_M <- str_replace(grep(CATG, ch1$header, value = TRUE),       # #
                                 paste0(catg, " "), "MASTER")                # #
      rm(catg)                                                               # #
                                                                             # #
      date_name <- get_fits_header_str(ch2$header, DATE)                     # #
                                                                             # #
      M_ch1_header <- ch1$header                                             # #
      M_ch2_header <- ch2$header                                             # #
                                                                             # #
      M_ch1_header[match(grep(CATG, M_ch1_header, value = TRUE),             # #
                         M_ch1_header)] <- catg_repl_M                       # #
      M_ch2_header[match(grep(CATG, M_ch2_header, value = TRUE),             # #
                         M_ch2_header)] <- catg_repl_M                       # #
      rm(catg_repl_M)                                                        # #
    }                                                                        # #
    rm(ch1, ch2)                                                             # #
                                                                             # #
    if(mv_flag){                                                             # #
      file.copy(pathList[flat_ind[temp_ind[f]]],                             # #
                paste0(flat_path, fileList[flat_ind[temp_ind[f]]]))          # #
      file.remove(pathList[flat_ind[temp_ind[f]]])                           # #
      file.copy(pathList[flat_ind[temp_ind[f + 1]]],                         # #
                paste0(flat_path, fileList[flat_ind[temp_ind[f + 1]]]))      # #
      file.remove(pathList[flat_ind[temp_ind[f + 1]]])                       # #
    }                                                                        # #
                                                                             # #
    print(paste0("Flats loaded for filter ", bands_str[k], ": ", f + 1, "/", # #
                 lengths_ind[k]))                                            # #
  }                                                                          # #
  rm(temp_ind, f, fi)                                                        # #
  ##############################################################################
                                                                               #
  #--------------------------------------------------------------------------# #
                                                                               #
  ######################## Creating MASTER FLAT arrays ####################### #
  ##############################################################################
  print(paste0("Combining loaded flats at ", bands_str[k], " for CH1..."))   # #
  band_masterFlat_ch1 <- array(NA, dim = c(dimsCH, 2),                       # #
                               dimnames = list(NULL, NULL, c("Data", "Unc")))# #
  band_masterFlat_ch1[,,"Data"] <- apply(t_ld1, 1:2, median, na.rm = TRUE)   # #
  band_masterFlat_ch1[,,"Unc"] <- unc_median(t_ld1, temp_U1, 1:2)            # #
  rm(temp_U1, t_ld1)                                                         # #
                                                                             # #
  print(paste0("Combining loaded flats at ", bands_str[k], " for CH2..."))   # #
  band_masterFlat_ch2 <- array(NA, dim = c(dimsCH, 2),                       # #
                               dimnames = list(NULL, NULL, c("Data", "Unc")))# #
  band_masterFlat_ch2[,,"Data"] <- apply(t_ld2, 1:2, median, na.rm = TRUE)   # #
  band_masterFlat_ch2[,,"Unc"] <- unc_median(t_ld2, temp_U2, 1:2)            # #
  rm(temp_U2, t_ld2)                                                         # #
  ##############################################################################
                                                                               #
  #--------------------------------------------------------------------------# #
                                                                               #
  ######################### Creating MASTER FLAT files ####################### #
  ##############################################################################
  writeFITSim(band_masterFlat_ch1,                                           # #
              file = paste0(Mflat_path,"FORS2-", date_name, "_Master-Flat-", # #
                            bands_str[k], "-CH1.fits"), crpixn = crpix1,     # #
              crvaln = crval1, ctypen = ctype, header = M_ch1_header)        # #
  writeFITSim(band_masterFlat_ch2,                                           # #
              file = paste0(Mflat_path,"FORS2-", date_name, "_Master-Flat-", # #
                            bands_str[k], "-CH2.fits"), crpixn = crpix2,     # #
              crvaln = crval2, ctypen = ctype, header = M_ch2_header)        # #
                                                                             # #
  print(paste0("Master FLAT chips for ", bands_str[k],": Created."))         # #
  rm(band_masterFlat_ch1, band_masterFlat_ch2, M_ch1_header, M_ch2_header,   # #
     crpix1, crpix2, crval1, crval2, ctype, date_name)                       # #
  ##############################################################################
}                                                                              #
rm(pathList, fileList, bands_str, Mflat_path, flat_ind, flat_path, bias_path,  #
   bias1_M, bias2_M, lengths_ind, CATG, CHIP, REFX, REFY, VALY, VALX, DATE, k, #
   GAIN, RON, bands, dimsCH, mv_flag)                                          #
################################################################################
################################################################################