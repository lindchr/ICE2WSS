# This file is part of ICE2WSS
# ICE2WSS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# [Full details at https://www.gnu.org/licenses/gpl-3.0.en.html]


##############################################3
#Extract data from SWORD
##############################################

#' Process all files specified by use in parallel.
#'
#' @param SWORD_dir Path pointing to folder with SWORD data.
#' @param Area String. Contains area of interest of SWORD data
#' @param Version String. SWORD data
#' @return Nothing. File is produced in SWORD subfolder
Produce_SWORD_data <- function(SWORD_dir,Area,Version){

   cat(substr(as.character(Sys.time()),1,16),"Loading SWORD data from hdf5 file.\n")

   SWORD_nc_dir <- paste(SWORD_dir,"/SWORD_",Version,"_nc/netcdf/",Area,"_sword_",Version,".nc",sep="")
   file <- hdf5r::h5file(SWORD_nc_dir)

   Lon_node <- file[["nodes/x"]][]
   Lat_node <- file[["nodes/y"]][]
   Node_id <- file[["nodes/node_id"]][]
   Reach_id <- file[["nodes/reach_id"]][]
   Wse <- file[["nodes/wse"]][]
   Node_dat <- cbind(Lon_node,Lat_node,Wse, Node_id,Reach_id)#,Node_length)
   Node_dat <- as.data.frame(Node_dat,row.names = FALSE,col.names = names(Node_dat))

   Rch_id_down <- file[["reaches/rch_id_dn"]][,1:4]
   Rch_id_up <- file[["reaches/rch_id_up"]][,1:4]
   Reach_slope <- file[["reaches/slope"]][]
   Reach_id <- file[["reaches/reach_id"]][]

   Reach_dat <- cbind(Reach_id, Reach_slope, Rch_id_up,Rch_id_down)
   Names <- c("Reach_id", "Reach_slope","Rch_id_up1","Rch_id_up2","Rch_id_up3","Rch_id_up4","Rch_id_down1","Rch_id_down2","Rch_id_down3","Rch_id_down4")
   colnames(Reach_dat) <- Names
   Reach_dat <- as.data.frame(Reach_dat,row.names = FALSE,col.names = Names)

   df_new <- merge(Node_dat, Reach_dat,all.x=TRUE)
   df_new$Area <- Area

   utils::write.table(df_new, file = paste(SWORD_dir,"/Processed_SWORD_",Version,"_",Area,".txt",sep="" ), row.names = FALSE, append = FALSE, col.names = TRUE, sep = ",",quote = FALSE)

}
