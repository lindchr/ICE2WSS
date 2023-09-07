##############################################3
#Extract data from SWORD
##############################################

#' Process all files specified by use in parallel.
#'
#' @param SWORD_dir Path pointing to folder with SWORD data.
#' @param Area String. Contains area of interest of SWORD data
#' @param Version String. SWORD data
#' @return Nothing. File is produced in SWORD subfolder
#' @export
Produce_SWORD_data <- function(SWORD_dir,Area,Version){

   cat(as.character(Sys.time()),"Reformatting SWORD data. This can take a while \n")
   #existing_dirs <- list.dirs(SWORD_dir,recursive=FALSE,full.names = FALSE)
   #path <- intersect(grep(Version, existing_dirs, value = TRUE),grep("nc", existing_dirs, value = TRUE))
   SWORD_nc_dir <- paste(SWORD_dir,"/SWORD_",Version,"_nc/netcdf/",Area,"_sword_",Version,".nc",sep="")
   file <- hdf5r::h5file(SWORD_nc_dir)

   Lon_node <- file[["nodes/x"]][]
   Lat_node <- file[["nodes/y"]][]
   Node_id <- file[["nodes/node_id"]][]
   Reach_id <- file[["nodes/reach_id"]][]
   Wse <- file[["nodes/wse"]][]
   #Node_length <- file[["nodes/node_length"]][] #Removed may 23
   Node_dat <- cbind(Lon_node,Lat_node,Wse, Node_id,Reach_id)#,Node_length)
   Node_dat <- as.data.frame(Node_dat,row.names = FALSE,col.names = names(Node_dat))

   #N_rch_down <- file[["reaches/n_rch_down"]][] #Removed may 23
   #N_rch_up <- file[["reaches/n_rch_up"]][]
   Rch_id_down <- file[["reaches/rch_id_dn"]][,1:4]
   Rch_id_up <- file[["reaches/rch_id_up"]][,1:4]
   Reach_slope <- file[["reaches/slope"]][]
   Reach_id <- file[["reaches/reach_id"]][]

   Reach_dat <- cbind(Reach_id, Reach_slope, Rch_id_up,Rch_id_down)#N_rch_up,N_rch_down
   #Names <- c("Reach_id", "Reach_slope", "N_rch_up","N_rch_down","Rch_id_up1","Rch_id_up2","Rch_id_up3","Rch_id_up4","Rch_id_down1","Rch_id_down2","Rch_id_down3","Rch_id_down4")
   Names <- c("Reach_id", "Reach_slope","Rch_id_up1","Rch_id_up2","Rch_id_up3","Rch_id_up4","Rch_id_down1","Rch_id_down2","Rch_id_down3","Rch_id_down4")
   colnames(Reach_dat) <- Names
   Reach_dat <- as.data.frame(Reach_dat,row.names = FALSE,col.names = Names)

   df_new <- merge(Node_dat, Reach_dat,all.x=TRUE)

   utils::write.table(df_new, file = paste(SWORD_dir,"/Processed_SWORD_",Version,"_",Area,".txt",sep="" ), row.names = FALSE, append = FALSE, col.names = TRUE, sep = ",",quote = FALSE)

}
