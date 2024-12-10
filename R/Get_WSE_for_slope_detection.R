# This file is part of ICE2WSS
# ICE2WSS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# [Full details at https://www.gnu.org/licenses/gpl-3.0.en.html]


#' Calculate River water surface slopes using satellite data and river centre lines
#'
#' Process all files specified by user in parallel. Takes water surface elevation data and assign to SWORD centre lines and nodes.
#'  Apply outlier filtering and divides data in o sections of the river depending on user settings.
#'  Computes water surface slope based on a regression of the distance along the river vs the water surface elevation.
#'  The water surface slope is the slope of the regression. Data is outputted in a text file with time, position, slope, and statistical metrics.
#'
#' @param Paths List of 3 paths: c(Outdir, ATL13_data_path, SWORD_dir)
#' @param Files List of files to process.
#'  Input data must contain csv files with comma as delimiter with the following columns:
#'  DecYear, Latitude, Longitude, Water_height, Water_ID, Beam.
#'  An additional column containing surface water occurrence in whole numbers between 0 and 100 is optional.
#' @param SWORD list. List containing the SWORD version and SWORD area
#' @param Max_reg_dist A number. Maximum distance along river accepted for slope calculation
#' @param Min_reg_dist A number. Minimum distance along river accepted for slope calculation
#' @param Min_reg_p A number. Minimum number of water surface elevation points needed for slope calculation
#' @param Occ_thr A number. Minimum water surface occurrence accepted. See https://global-surface-water.appspot.com/
#' @importFrom foreach %dopar%
#' @return Nothing. File is produced in output path
#' @export
ICE2WSS <- function(Paths, Files, SWORD, Max_reg_dist=8000, Min_reg_dist = 400,
                     Min_reg_p=10 , Occ_thr=NA){
   SWORD_dat <- NA
   i <- NA
   output_file <- NA

   if (missing(Paths) == TRUE) {stop(paste(as.character(Sys.time()),"Error: Required paths are not provided. See documentation."))}
   if (missing(Files)== TRUE) {stop(paste(as.character(Sys.time()),"Error: List of file to process is not provided."))}
   if (missing(SWORD)== TRUE) {stop(paste(as.character(Sys.time()),"Error: No specifications for SWORD area and version is specified"))}
   if (missing(Max_reg_dist)== TRUE) {warning(paste(as.character(Sys.time()),"Warning: No maximum distance for slope estimate is provided. Using default value Max_reg_dist = 8000m."))}
   if (missing(Min_reg_dist)== TRUE) {warning(paste(as.character(Sys.time()),"Warning: No minimum distance for slope estimate is provided. Using default value Min_reg_dist = 400m."))}
   if (missing(Min_reg_p)== TRUE) {warning(paste(as.character(Sys.time()),"Warning: No minimum number of data points for slope estimate is provided. Using default value Min_reg_p = 10."))}
   if (missing(Occ_thr)== TRUE) {warning(paste(as.character(Sys.time()),"Warning: No value for water occurrence is provided. Using default value Occ_thr = NA"))}

   if (!is.numeric(Max_reg_dist)) {stop(paste(as.character(Sys.time()),"Error: Maximum distance for slope estimate is not a numeric value."))}
   if (!is.numeric(Min_reg_dist)) {stop(paste(as.character(Sys.time()),"Error: Minimum distance for slope estimate is not a numeric value."))}
   if (!Min_reg_p == round(Min_reg_p)) {stop(paste(as.character(Sys.time()),"Error: Minimum number of points for slope estimate is not an integer."))}

   if (!is.character(Paths) | !length(Paths)==3) {stop(paste(as.character(Sys.time()),"Error: Paths is not character of required folder paths or does not contain 3 elements (Outdir, ATL13dir,SWORDdir)."))}
   if (!is.character(SWORD) | !length(SWORD)==2) {stop(paste(as.character(Sys.time()),"Error: SWORD is not character of required SWORD info or does not contain 2 elements (version, area)."))}
   if (!is.numeric(Files)) {stop(paste(as.character(Sys.time()),"Error: Files does not contains indices for files to process."))}
   if (! SWORD[2] %in% c("as", "oc", "na", "sa", "af", "eu")) {stop(paste(as.character(Sys.time()),"SWORD area invalid."))}
   options(scipen=999)

   #Check paths and create output file
   closeAllConnections()
   if (!dir.exists(Paths[2])) {stop(paste(as.character(Sys.time()),"Error: Path to input ATL13 data does not exist."))}
   if (!dir.exists(Paths[3])) {stop(paste(as.character(Sys.time()),"Error: Path to input SWORD data does not exist."))}
   ifelse(!dir.exists(Paths[1]), dir.create(Paths[1]), FALSE)

   #### Process files
   filelist <- list.files(Paths[2],full.names = TRUE)
   if (max(Files) > length(filelist)) {stop(paste(as.character(Sys.time()),"Error: Files to process exceeds maximum files available."))}


   Colnames <- data.frame(matrix(c("DecYear","Lon [deg]","Lat [deg]","WSS [cm/km]","StdError [cm/km]","IntersectAngle [deg]","Nobs",
                                 "Resolution [m]","pVal","R2","WSE","WaterID","Beam","ReachID","NodeID","ReachDEMSlope [cm/km]")   ,nrow=1))

   file_timestamp <- format(Sys.time(), "%d%m%Y_%H%M")
   output_file <<- paste(Paths[1],"/WSS_",file_timestamp,".txt",sep="")
   cat(substr(as.character(Sys.time()),1,16),"Results written to output file: ", output_file,"\n")
   utils::write.table(Colnames, file = output_file, row.names = FALSE, append = FALSE, col.names = FALSE, sep = ",",quote = FALSE)

   Log_file <- paste(Paths[1],"/Log_",file_timestamp,".txt",sep="")
   cat(substr(as.character(Sys.time()),1,16),"Progress written to log file: ", Log_file,"\n")
   utils::write.table(paste0("Current settings: Paths=",Paths[1],", ",Paths[2],", ",Paths[3], ", SWORD=",SWORD[1]," ",SWORD[2], ", Max_reg_dist=",Max_reg_dist,", Min_reg_dist=",Min_reg_dist,", Min_reg_p",Min_reg_p,", Occ_thr=",Occ_thr), file = Log_file, row.names = FALSE, append = FALSE, col.names = FALSE, sep = ",",quote = FALSE)
   utils::write.table(" ", file = Log_file, row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",",quote = FALSE)

   SWORD_dat <- handle_SWORD(Paths[3], SWORD[1], SWORD[2],Log_file)

   #Cores to be used: All-1
   C <- parallel::detectCores(all.tests = FALSE, logical = TRUE) -1
   cl <- parallel::makeCluster(C,outfile = Log_file)
   utils::write.table(paste(substr(as.character(Sys.time()),1,16),"Running program with", length(cl),"threads \n"), file = Log_file, row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",",quote = FALSE)
   doParallel::registerDoParallel(cl)
   parallel::clusterExport(cl, "SWORD_dat")
   parallel::clusterExport(cl, "output_file")
   #parallel::clusterExport(cl, "Log_file")
   foreach::foreach(i=Files, .packages=c("ICE2WSS")) %dopar% {
      ICE2WSS::Find_slope(Paths, i, SWORD_dat, Max_reg_dist, Min_reg_dist,Min_reg_p, Occ_thr,filelist,output_file,SWORD_dat)
   }
   parallel::stopCluster
   parallel::stopCluster(cl)

   # lapply(Files, function(Files)Find_slope(Paths, Files, SWORD, Max_reg_dist, Min_reg_dist,
   #                                                    Min_reg_p, Occ_thr,filelist,output_file))

   utils::write.table(paste(substr(as.character(Sys.time()),1,16),"Computations completed"), file = Log_file, row.names = FALSE, append = TRUE, col.names = FALSE, sep = ",",quote = FALSE)
   cat(substr(as.character(Sys.time()),1,16),"Computations completed \n")

}

