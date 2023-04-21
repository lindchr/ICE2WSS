#####################################
#Extract slope from ATL13 data
#Load .dat files produced by sript Extract_ATL13
#Choose data representing the river, and calculate the slope
###################################


#' Process all files specified by use in parallel.
#'
#' @param Paths List of 3 paths: c(Outdir, ATL13_data_path, SWORD_dir)
#' @param Files List of files to process.
#' @param SWORD list. List containing the SWORD version and SWORD area
#' @param Max_reg_dist A number. Maximum distance along river accepted for slope calculation
#' @param Min_reg_dist A number. Minimum distance along river accepted for slope calculation
#' @param Min_reg_p A number. Minimum number of water surface elevation points needed for slope calculation
#' @param Occ_thr A number. Minimum water surface occurrence accepted. See https://global-surface-water.appspot.com/
#' @return Nothing. File is produced in output path
ICE2WSS <- function(Paths, Files, SWORD, Max_reg_dist, Min_reg_dist,
                     Min_reg_p, Occ_thr){

   Max_reg_dist <<- Max_reg_dist
   Min_reg_dist <<- Min_reg_dist
   Min_reg_p <<- Min_reg_p
   Occ_thr <<- Occ_thr

   # Handle needed packages
   my_packages <- c("sp", "rgdal", "hdf5r")                          # Specify your packages
   not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])]    # Extract not installed packages
   if(length(not_installed)) install.packages(not_installed)

   library(sp)
   library(rgdal)
   options(scipen=999)
   library(hdf5r)


   closeAllConnections()
   ifelse(!dir.exists(Paths[1]), dir.create(Paths[1]), FALSE)

   if(is.na(Occ_thr)){
      Colnames <- c("DecYear","Lon [deg]","Lat [deg]","WSS [cm/km]","StdError [cm/km]","IntersectAngle [deg]","Nobs",
            "Resolution [m]","pVal","R2","WSE","WaterID","Beam","ReachID","NodeID","ReachDEMslope [cm/km]")
   } else {
      Colnames <- data.frame(matrix(c("DecYear","Lon [deg]","Lat [deg]","WSS [cm/km]","StdError [cm/km]","IntersectAngle [deg]","Nobs",
            "Resolution [m]","pVal","R2","WSE","WaterID","Beam","ReachID","NodeID","ReachDEMSlope [cm/km]")   ,nrow=1))
   }
   file_timestamp <<- format(Sys.time(), "%d%m%Y_%H%M")
   output_file <<- paste(Paths[1],"WSS_",file_timestamp,".txt",sep="")
   write.table(Colnames, file = output_file, row.names = FALSE, append = FALSE, col.names = FALSE, sep = ", ",quote = FALSE)

   handle_SWORD(Paths[3], SWORD[1], SWORD[2])

   #### Process files
   filelist <<- list.files(Paths[2],full.names = TRUE)
   file <- c(Files[1]:Files[2])

   start_time <- Sys.time()

   #output <- parallel::mclapply(file, function(file)Find_slope(file))
   output <- lapply(file, function(file)Find_slope(file))

   end_time <- Sys.time()
   cat(end_time - start_time)

}

