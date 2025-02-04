rm(list=ls())
require(FITSio)
require(stringr)
require(jjb)
require(profmem)

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
bias_path <- paste0(input_folder, "BIAS/")                                     #
print("Provide a path to move the BIAS files to: ")                            #
bias_path <- edit(bias_path)                                                   #
                                                                               #
if(!dir.exists(bias_path)){                                                    #
  mkdir(bias_path)                                                             #
}                                                                              #
                                                                               #
mv_flag <- 1                                                                   #
                                                                               #
if(input_folder == bias_path){                                                 #
  mv_flag <- 0                                                                 #
}                                                                              #
                                                                               #
Mbias_path <- paste0(input_folder, "MASTER_BIAS/")                             #
print("Provide an output path: ")                                              #
Mbias_path <- edit(Mbias_path)                                                 #
                                                                               #
if(!dir.exists(Mbias_path)){                                                   #
  mkdir(Mbias_path)                                                            #
}                                                                              #
################################################################################

##_##################################_###_###################################_##
#/ \--------------------------------/ \|/ \---------------------------------/ \#
#\_/--------------------------------\_/|\_/---------------------------------\_/#
################################################################################

################################################################################
#################### Setting up files lists and parameters #####################
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
GAIN <- "HIERARCH ESO DET OUT1 GAIN"                                           #
RON <- "HIERARCH ESO DET OUT1 RON"                                             #
                                                                               #
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
#Possible types of fits files: "BIAS", "DARK", "FLAT", "OBJECT", "SKY", "STD"  #
types_arr <- get_fits_header_list_str(pathList, TYPE)                          #
rm(TYPE)                                                                       #
                                                                               #
#Standard chip dimensions for FORS2 CCDs                                       #
dimsChip <- c(2048, 1034)                                                      #
                                                                               #
#Create list of pathList indexes for BIASes fits                               #
bias_ind <- get_fits_header_key_index(types_arr, pathList, dimsChip, "BIAS")   #
rm(types_arr)                                                                  #
                                                                               #
bias_n <- length(bias_ind) / 2                                                 #
                                                                               #
dummy_arr <- 1:bias_n                                                          #
ram_req <- 1.1 * sum(profmem(median(dummy_arr))$bytes) * dimsChip[1] *         #
  dimsChip[2]                                                                  #
rm(dummy_arr)                                                                  #
################################################################################

##_##################################_###_###################################_##
#/ \--------------------------------/ \|/ \---------------------------------/ \#
#\_/--------------------------------\_/|\_/---------------------------------\_/#
################################################################################

################################################################################
##### Loading all bias files to an array to create master bias for CHIP 1 ######
################################################################################
print("Creating Master Bias arrays for CH1 and CH2")                           #
masterBias_ch1 <- array(NA, dim = c(dimsChip, 2),                              #
                        dimnames = list(NULL, NULL, c("Data", "Unc")))         #
masterBias_ch2 <- array(NA, dim = c(dimsChip, 2),                              #
                        dimnames = list(NULL, NULL, c("Data", "Unc")))         #
                                                                               #
# Making sure that there's enough ram to process each iteration                #
abort_prevention <- 5                                                          #
rf <- c(0, round(seq(dimsChip[1] / abort_prevention, dimsChip[1],              #
                 length.out = abort_prevention)))                              #
cf <- c(0, round(seq(dimsChip[2] / abort_prevention, dimsChip[2],              #
                 length.out = abort_prevention)))                              #
                                                                               #
for(r in 2:length(rf)){                                                        #
  rs <- (rf[r - 1] + 1):rf[r]                                                  #
  for(c in 2:length(cf)){                                                      #
    cs <- (cf[c - 1] + 1):cf[c]                                                #
                                                                               #
    lr <- length(rs)                                                           #
    lc <- length(cs)                                                           #
                                                                               #
    loadBias_px1 <- array(NA, dim = c(lr, lc, bias_n, 2),                      #
                          dimnames = list(NULL, NULL, NULL, c("Data","Unc")))  #
    loadBias_px2 <- array(NA, dim = c(lr, lc, bias_n, 2),                      #
                          dimnames = list(NULL, NULL, NULL, c("Data","Unc")))  #
                                                                               #
    for(b in seq(1, length(bias_ind), 2)){                                     #
      #Loading the files                                                       #
      tempA <- readFITS(pathList[bias_ind[b]])                                 #
      tempB <- readFITS(pathList[bias_ind[b + 1]])                             #
                                                                               #
      chipA <- get_fits_header_str(tempA$header, CHIP)                         #
                                                                               #
      if(chipA == "CHIP1"){                                                    #
        ch1 <- tempA                                                           #
        ch2 <- tempB                                                           #
      }else{                                                                   #
        ch1 <- tempB                                                           #
        ch2 <- tempA                                                           #
      }                                                                        #
      rm(chipA, tempA, tempB)                                                  #
                                                                               #
      #Loading each pixel in the corresponding load array.                     #
      print(paste0("-> Loading chip data of pixels in rows: ", rs[1], ":",     #
                   rs[lr], "/", dimsChip[1], " and cols: ", cs[1], ":", cs[lc],#
                   "/", dimsChip[2], " of bias nº ", (b + 1) / 2, " to temp a",#
                   "rrays..."))                                                #
      loadBias_px1[,, (b + 1) / 2, "Data"] <- ch1$imDat[rs, cs]                #
      loadBias_px2[,, (b + 1) / 2, "Data"] <- ch2$imDat[rs, cs]                #
                                                                               #
      #Handling header information for master bias files                       #
      if(rf[r] == dimsChip[1] & cf[c] == dimsChip[2] &                         #
         b == (length(bias_ind) - 1)){                                         #
        print("-> Creating headers...")                                        #
                                                                               #
        refx_1 <- get_fits_header_num(ch1$header, REFX)                        #
        refy_1 <- get_fits_header_num(ch1$header, REFY)                        #
        refz <- 1                                                              #
        crpix1 <- c(refx_1, refy_1, refz)                                      #
        rm(refx_1, refy_1)                                                     #
                                                                               #
        refx_2 <- get_fits_header_num(ch2$header, REFX)                        #
        refy_2 <- get_fits_header_num(ch2$header, REFY)                        #
        crpix2 <- c(refx_2, refy_2, refz)                                      #
        rm(refx_2, refy_2)                                                     #
                                                                               #
        valx_1 <- get_fits_header_num(ch1$header, VALX)                        #
        valy_1 <- get_fits_header_num(ch1$header, VALY)                        #
        valz <- NA                                                             #
        crval1 <- c(valx_1, valy_1, valz)                                      #
        rm(valy_1, valx_1)                                                     #
                                                                               #
        valx_2 <- get_fits_header_num(ch2$header, VALX)                        #
        valy_2 <- get_fits_header_num(ch2$header, VALY)                        #
        crval2 <- c(valx_2, valy_2, valz)                                      #
        rm(valy_2, valx_2)                                                     #
                                                                               #
        typex <- "  'PIXEL     '              / Coordinate system of x-axis "  #
        typey <- "  'PIXEL     '              / Coordinate system of y-axis "  #
        typez <- "  'DATA & UNC'              / Coordinate system of z-axis "  #
        ctype <- c(typex, typey, typez)                                        #
        rm(typex, typey)                                                       #
                                                                               #
        catg_1 <- get_fits_header_str(ch1$header, CATG)                        #
        catg_2 <- get_fits_header_str(ch2$header, CATG)                        #
                                                                               #
        catg_repl_M1 <- str_replace(grep(CATG, ch1$header, value = TRUE),      #
                                    paste0(catg_1, " "), "MASTER")             #
        catg_repl_M2 <- str_replace(grep(CATG, ch2$header, value = TRUE),      #
                                    paste0(catg_2, " "), "MASTER")             #
        rm(catg_1, catg_2)                                                     #
                                                                               #
        date_name <- get_fits_header_str(ch2$header, DATE)                     #
        rm(DATE)                                                               #
                                                                               #
        M_ch1_header <- ch1$header                                             #
        M_ch2_header <- ch2$header                                             #
                                                                               #
        M_ch1_header[match(grep(CATG, M_ch1_header, value = TRUE),             #
                           M_ch1_header)] <- catg_repl_M1                      #
        M_ch2_header[match(grep(CATG, M_ch2_header, value = TRUE),             #
                           M_ch2_header)] <- catg_repl_M2                      #
        rm(catg_repl_M1, catg_repl_M2)                                         #
      }                                                                        #
    }                                                                          #
                                                                               #
    gain1 <- get_fits_header_num(ch1$header, GAIN)                             #
    gain2 <- get_fits_header_num(ch2$header, GAIN)                             #
    ron1 <- get_fits_header_num(ch1$header, RON)                               #
    ron2 <- get_fits_header_num(ch2$header, RON)                               #
    rm(ch1, ch2)                                                               #
                                                                               #
    loadBias_px1[,,,"Unc"] <- sqrt(loadBias_px1[,,,"Data"] * gain1 +           #
                                     (ron1 * gain1)^2)                         #
    loadBias_px2[,,,"Unc"] <- sqrt(loadBias_px2[,,,"Data"] * gain2 +           #
                                     (ron2 * gain2)^2)                         #
    rm(gain1, ron1, gain2, ron2)                                               #
                                                                               #
    print(paste0("Calculating master bias pixels in rows: ", rs[1], ":",       #
                 rs[lr], "/", dimsChip[1], " and cols: ", cs[1], ":", cs[lc],  #
                 "/", dimsChip[2], "..."))                                     #
    masterBias_ch1[rs, cs, "Data"] <- apply(loadBias_px1[,,,"Data"], 1:2,      #
                                            median, na.rm = T)                 #
    masterBias_ch1[rs, cs, "Unc"] <- unc_median(loadBias_px1[,,,"Data"],       #
                                                loadBias_px1[,,,"Unc"])        #
    masterBias_ch2[rs, cs, "Data"] <- apply(loadBias_px2[,,, "Data"], 1:2,     #
                                            median, na.rm = T)                 #
    masterBias_ch2[rs, cs, "Unc"] <- unc_median(loadBias_px2[,,,"Data"],       #
                                                loadBias_px2[,,,"Unc"])        #
  }                                                                            #
}                                                                              #
rm(r, rs, rf, c, cs, cf, abort_prevention, ram_req, lr, lc,                    #
   loadBias_px1, loadBias_px2)                                                 #
                                                                               #
if(mv_flag){                                                                   #
  for(b in 1:length(bias_ind)){                                                #
    file.copy(pathList[bias_ind[b]], paste0(bias_path, fileList[bias_ind[b]])) #
    file.remove(pathList[bias_ind[b]])                                         #
    pathList[bias_ind[b]] <- paste0(bias_path, fileList[bias_ind[b]])          #
  }                                                                            #
}                                                                              #
rm(bias_ind, pathList, fileList, b, mv_flag)                                   #
                                                                               #
#Cutting Margins                                                               #
masterBias_ch1[,964:dimsChip[2],] <- NA                                        #
masterBias_ch1[,1:7,] <- NA                                                    #
masterBias_ch1[1:184,,] <- NA                                                  #
masterBias_ch1[1860:dimsChip[1],,] <- NA                                       #
masterBias_ch2[,1028:dimsChip[2],] <- NA                                       #
masterBias_ch2[,1:317,] <- NA                                                  #
masterBias_ch2[1:192,,] <- NA                                                  #
masterBias_ch2[1862:dimsChip[1],,] <- NA                                       #
rm(dimsChip, CATG, CHIP, GAIN, REFX, REFY, refz, RON, typez, VALX, VALY, valz) #
################################################################################

##_##################################_###_###################################_##
#/ \--------------------------------/ \|/ \---------------------------------/ \#
#\_/--------------------------------\_/|\_/---------------------------------\_/#
################################################################################

################################################################################
########################## Creating MASTER BIAS files ##########################
################################################################################
writeFITSim(masterBias_ch1, file = paste0(Mbias_path,"FORS2-", date_name,      #
                                          "_Master-Bias-CH1.fits"),            #
            crvaln = crval1, crpixn = crpix1, ctypen = ctype,                  #
            header = M_ch1_header)                                             #
writeFITSim(masterBias_ch2, file = paste0(Mbias_path,"FORS2-", date_name,      #
                                          "_Master-Bias-CH2.fits"),            #
            crvaln = crval2, crpixn = crpix2, ctypen = ctype,                  #
            header = M_ch2_header)                                             #
                                                                               #
print("Master BIAS chips: Created.")                                           #
rm(masterBias_ch1, masterBias_ch2, M_ch1_header, M_ch2_header, Mbias_path,     #
   crpix1, crpix2, crval1, crval2, ctype, date_name)                           #
################################################################################
################################################################################