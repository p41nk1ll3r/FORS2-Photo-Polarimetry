require(stringr)

###################################################################
####### To retrieve a numeric parameter from a FITS header ########
###################################################################
get_fits_header_num <- function(header, KEY = "CRPIX1"){          #
  num_val <- grep(KEY, header, value = TRUE)[1]                   #
  num_val <- strsplit(num_val,"=")[[1]][2]                        #
  num_val <- strsplit(num_val,"/")[[1]][1]                        #
  return(as.numeric(num_val))                                     #
}                                                                 #
###################################################################

########################################################################
########## To retrieve a string parameter from a FITS header ###########
########################################################################
get_fits_header_str <- function(header, KEY = "HIERARCH ESO DPR CATG"){#
  str_val <- strsplit(grep(KEY, header, value = TRUE)[1],"'")[[1]][2]  #
                                                                       #
  f_i <- unlist(gregexpr('[a-z+A-Z+0-9]', str_val))[1]                 #
  l_i <- tail(unlist(gregexpr('[a-z+A-Z+0-9]', str_val)), n=1)         #
                                                                       #
  str_val <- str_sub(str_val, f_i, l_i)                                #
  return(str_val)                                                      #
}                                                                      #
########################################################################

########################################################################
##################### To retrieve a parameter and ######################
################# its description from a FITS header ###################
########################################################################
get_fits_header_val <- function(header, KEY = "HIERARCH ESO DPR CATG"){#
  val <- sub(".*?=","", grep(KEY, header, value = TRUE)[1])            #
  val <- str_remove_all(substr(val, 1, min(c(68, nchar(val)))), "'")   #
  return(val)                                                          #
}                                                                      #
########################################################################

##############################################################
########### To retrieve a list of the same string ############
######## parameter from different FITS files headers #########
##############################################################
get_fits_header_list_str <- function(paths, KEY = "EXTNAME"){#
                                                             #
  count <- length(paths)                                     #
                                                             #
  types <- array(NaN, dim = count)                           #
                                                             #
  print(paste0("Listing ", KEY, "..."))                      #
                                                             #
  for(fx in 1:count){                                        #
                                                             #
    temp_fits <- readFITS(paths[fx])                         #
    temp_type <- get_fits_header_str(temp_fits$header, KEY)  #
                                                             #
    types[fx] <- temp_type                                   #
                                                             #
    rm(temp_fits, temp_type)                                 #
  }                                                          #
  return(types)                                              #
}                                                            #
##############################################################

##############################################################
########### To retrieve a list of the same numeric ###########
######### parameter from different FITS files headers ########
##############################################################
get_fits_header_list_num <- function(paths, KEY = "CRPIX1"){ #
                                                             #
  count <- length(paths)                                     #
                                                             #
  values <- array(NaN, dim = count)                          #
                                                             #
  print(paste0("Listing ", KEY, "..."))                      #
                                                             #
  for(fx in 1:count){                                        #
                                                             #
    temp_fits <- readFITS(paths[fx])                         #
    temp_val <- get_fits_header_num(temp_fits$header, KEY)   #
                                                             #
    values[fx] <- temp_val                                   #
                                                             #
    rm(temp_fits, temp_val)                                  #
  }                                                          #
  return(values)                                             #
}                                                            #
##############################################################

################################################################################
####### To retrieve a list of indexes for FITS files with a given value ########
###### on a specific parameter. Assumes that every FITS is part of a duo #######
######## and rejects FITS files that do not comply with that assumption ########
################################################################################
get_fits_header_key_index <- function(types, paths, dimExp, KEY = "BIAS"){     #
                                                                               #
  dimExp <- as.integer(dimExp)                                                 #
                                                                               #
  if(length(KEY) != 1){                                                        #
    print(paste0("ERROR: you have inserted ", length(KEY), "keys. Expected 1.",#
                 " Returning NULL."))                                          #
    return(NULL)                                                               #
  }                                                                            #
  if(length(types) != length(paths)){                                          #
    print(paste0("ERROR: Number of files is ", length(paths), " which is not ",#
                 "compatible with the number type flags provided (",           #
                 length(types), "). Returning NULL."))                         #
    return(NULL)                                                               #
  }                                                                            #
  if(length(dimExp) != 2){                                                     #
    print(paste0("ERROR: dimExp has length ", length(dimExp), " instead of 2.",#
                 " Returning NULL."))                                          #
    return(NULL)                                                               #
  }                                                                            #
                                                                               #
  CHIP <- "EXTNAME"                                                            #
                                                                               #
  prelim_inds <- which(types == KEY, arr.ind = TRUE)                           #
  prelength <- length(prelim_inds)                                             #
  inds <- prelim_inds                                                          #
                                                                               #
  to_rm <- NULL                                                                #
                                                                               #
  p <- 1                                                                       #
                                                                               #
  while(p <= prelength){                                                       #
                                                                               #
    if(p == prelength){                                                        #
      print("No pair available for the last file. Skiping file.")              #
      to_rm <- append(to_rm, p)                                                #
      p <- p + 1                                                               #
      next                                                                     #
    }                                                                          #
                                                                               #
    #Checking if the consecutive files share the name structure                #
    nameA <- strsplit(fileList[prelim_inds[p]], ".[0-9]+.fits")[[1]][1]        #
    nameB <- strsplit(fileList[prelim_inds[p + 1]], ".[0-9]+.fits")[[1]][1]    #
                                                                               #
    if(nameA != nameB){                                                        #
      print(paste0("Consecutive files ", nameA, " and ", nameB, " are not ",   #
                   "from the same exposure. Skiping files."))                  #
      to_rm <- append(to_rm, p)                                                #
      p <- p + 1                                                               #
      next                                                                     #
    }                                                                          #
                                                                               #
    #Loading the files                                                         #
    chA <- readFITS(paths[prelim_inds[p]])                                     #
    chB <- readFITS(paths[prelim_inds[p + 1]])                                 #
                                                                               #
    chipA <- get_fits_header_str(chA$header, CHIP)                             #
    chipB <- get_fits_header_str(chB$header, CHIP)                             #
                                                                               #
    if(!((chipA == "CHIP1" && chipB == "CHIP2") ||                             #
         (chipB == "CHIP1" && chipA == "CHIP2"))){                             #
      print(paste0("An anomaly was detected in ", nameA, " and ", nameB, ". ", #
                   "CHIP data is not the expected from files of the same ",    #
                   "exposure: CHIP_A = ", chipA, " and CHIP_B = ", chipB, ", ",#
                   "instead of CHIP_A, CHIP_B = {1,2} & CHIP_A != CHIP_B. ",   #
                   "Skipping Files."))                                         #
      to_rm <- append(to_rm, c(p, p + 1))                                      #
      p <- p + 2                                                               #
      next                                                                     #
    }                                                                          #
                                                                               #
    #Regarding the dimensions of the pixel map of each CHIP                    #
    #And the CHIP file order                                                   #
    dimA <- dim(chA$imDat)                                                     #
    dimB <- dim(chB$imDat)                                                     #
                                                                               #
    #Additionally, if both chip share the acquisition procedure                #
    if(!identical(dimA, dimB)){                                                #
      print(paste0("An anomaly was detected in ", nameA, " and ", nameB, ". ", #
                   "CHIP data have different dimensions: (", dimA[1], ", ).",  #
                   dimA[2], ") and (", dimB[1], ", ", dimB[2], ". They should",#
                   " be the same. Skipping files."))                           #
      to_rm <- append(to_rm, c(p, p + 1))                                      #
    }                                                                          #
    if(!identical(dimA, dimExp)){                                              #
      print(paste0("An anomaly was detected, ", nameA, " and ", nameB, " data",#
                   " have dimensions: (", dimA[1], ", ", dimA[2], ") which is",#
                   " not the expected (", dimExp[1], ", ", dimExp[2], "). ",   #
                   "Skipping files."))                                         #
      to_rm <- append(to_rm, c(p, p + 1))                                      #
    }                                                                          #
                                                                               #
    p <- p + 2                                                                 #
    rm(chA, chB)                                                               #
  }                                                                            #
  rm(p, prelength)                                                             #
                                                                               #
  if(!is.null(to_rm)){                                                         #
    inds <- prelim_inds[-to_rm]                                                #
  }else{                                                                       #
    inds <- prelim_inds                                                        #
  }                                                                            #
                                                                               #
  return(inds)                                                                 #
}                                                                              #
################################################################################