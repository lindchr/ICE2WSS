############################
#Functions to extract river slope from icesat ATL13
# Find Slope loads data and extract occurrence > 85
# Assign_SWORD Project all measurements onto centerline and assignmen SWORD data to satellite data
# The end of Find slope calculates the river slope for WSE measured at the same time using linear regression
# within 2 and 10 km.
###########################

#' Find_slope
#' Checks input data, call all other functions and saves output data to file
#'
#' @param file2 Number. File number in satellite input data folder.
#' @return Nothing. File is produced in output path
#' Find_slope()
Find_slope <- function(file2){

   if(file2 %% 10 == 0){
      cat(as.character(Sys.time()),"Progress: Running file",file2,"\n") # in parallel with other files. \n")
   }

   #Current_file <<- file2 #Used for extracting method data

   test <- try(read.table(filelist[file2],sep=",",header=FALSE),silent=TRUE)

   if(class(test) == "try-error"){
      cat(substr(as.character(Sys.time()),1,16), "Load input data failed for file",filelist[file2],".\n",test[1],"\n")
      return()
   } else if (file.info(filelist[file2])$size == 0){
      cat(substr(as.character(Sys.time()),1,16), "Load input data failed for file",filelist[file2],". File is empty. \n")
      return()
   } else {
      dat <- read.table(filelist[file2],sep=",",header=FALSE)
   }
   dat <- read.table(filelist[file2],sep=",",header=FALSE)

   ##### Check data formats

   test <- try((dim(dat)[2] == 6 & is.na(Occ_thr))| (dim(dat)[2] == 7 & !is.na(Occ_thr)),silent=TRUE)

   if(test ==FALSE){
      cat(substr(as.character(Sys.time()),1,16), "Load input data failed for file",filelist[file2],". Wrong number of columns:",dim(dat)[2]," \n")
      return()
   }

   if(dim(dat)[2] == 7 & !is.na(Occ_thr)){ #Check occurrence values
      dat <- dat[dat[,7] > Occ_thr,1:7]
      colnames(dat) <- c("DecYear", "Lat", "Lon","H_ortho","WaterID","Beam","Occurrence")
      if(any(dat$Occurrence < 0)) {warning(paste(as.character(Sys.time()),"Warning: Water occurrence is less than 0. Check input data."))}
      if(any(dat$Occurrence > 100)) {warning(paste(as.character(Sys.time()),"Warning: Water occurrence is larger than 100. Check input data."))}

      if(dim(dat)[1] < Min_reg_p){
         return()
      }

      if(dim(dat)[2] < 5 & dim(dat)[1] == 1){
         return()
      }

   } else { #If occurrence is not included and we are happy
      colnames(dat) <- c("DecYear", "Lat", "Lon","H_ortho","WaterID","Beam")
   }

   #if(max(dat$DecYear)- min(dat$DecYear) > 0.00025){
   #   warning(paste(as.character(Sys.time()),"Warning: Data may not contain multiple passes. Time duration covered in file is more than 2 hours."))
   #}

   if(any(dat$DecYear < 0)) {warning(paste(as.character(Sys.time()),"Warning: Decimal year contains negative values. Ensure that input data are correct."))}
   if(any(abs(dat$Lat) > 90)) {warning(paste(as.character(Sys.time()),"Warning: Latitude contains values larger than +/-90. Ensure that input data are correct."))}
   if(any(abs(dat$Lon) > 360)) {warning(paste(as.character(Sys.time()),"Warning: Longitude contains values larger than 360. Ensure that input data are correct."))}
   if(any(dat$Beam > 6 | dat$Beam < 1)) {warning(paste(as.character(Sys.time()),"Warning: Beam number is not value between 1 and 6. Ensure that input data are correct."))}

   i <- c(1:6)
   beams_out <- lapply(i,function(i)Assign_SWORD(dat,i))
   beams_out <- as.data.frame(sapply(do.call(rbind, beams_out),as.numeric))

   if(dim(beams_out)[1] == 35 & dim(beams_out)[2] == 1 | dim(beams_out)[1] == 36 & dim(beams_out)[2] == 1){
      return()
   }



   if(dim(beams_out)[1] < Min_reg_p){
      return()
   }

   if(dim(beams_out)[2] == 1 & dim(beams_out)[1] == 32){
      return()
   }


   beams_out <- beams_out[order(beams_out$Node_id),]

   #Calculate UTM corr to get distance along centerline in meters
   UTM_zone <- round((beams_out$centerline_lon[1]+180)/6)
   myProj<-as.character(paste("+proj=utm +zone=",UTM_zone," +north ellps=WGS84",sep=""))

   xydat<-project(as.matrix(cbind(beams_out$centerline_lon,beams_out$centerline_lat)), as.character(myProj))
   beams_out$UTM_E <- xydat[,1]
   beams_out$UTM_N <- xydat[,2]

   beams_out <- Filter_data(beams_out)
   if(is.null(beams_out)){
      return()
   }

   slope_output <- Calc_slope(beams_out,myProj)

   if(is.null(slope_output)){
      return()
   }

   if(any(is.na(slope_output$WSS))){
      id <- which(is.na(slope_output$WSS))
      output <- slope_output[-c(id),]
   } else {
      output <- slope_output
   }
   output$Reach_slope <- format(round(output$Reach_slope*100,digits=4),nsmall = 4)

   output <- output[,c("DecYear","Slope_lon","Slope_lat","WSS","StdError","IntersectAngle","Nobs",
                  "Resolution","pVal","R2","Wse","WaterID","Beam","Reach_id","Node_id","Reach_slope")]

   output[,c(1,2,3)] <- format(round(output[,c(1,2,3)],digits=6),nsmall = 6)
   output[,6] <- format(round(as.numeric(output[,6]),digits=1),nsmall = 1)
   output[,c(4,5,8,9,10,11)] <- format(round(output[,c(4,5,8,9,10,11)], digits=4), nsmall = 4)

   write.table(output, file = output_file, row.names = FALSE, append = TRUE, col.names = FALSE, sep = ", ",quote = FALSE)
   return(output)
}


#' Process data loaded in Find_slope. Assigns SWORD data to each satellite measurement.
#'
#' @param dat Main data frame containing satellite data loaded in Find_slope
#' @param i List of beams (1:6).
#' @return Data frame containing original data and additional columns of SWORD data.
Assign_SWORD <- function(dat,i){
   dat_beam <- dat[dat$Beam == i,]

   if(dim(dat_beam)[1] < 1){
      return()
   }

   knn <- as.data.frame(RANN::nn2(SWORD_dat[,c(2,3)],dat_beam[,c(3,2)],k=1))
   knn <- knn[!(knn$nn.idx=0), ] #Remove end nodes - rest of script does not work for these.

   dat_beam <- cbind(dat_beam,SWORD_dat[knn$nn.idx,])
   dat_beam$type <- as.numeric(substr(SWORD_dat$Node_id[knn$nn.idx],14,15))

   #next_VS_idx <- knn$nn.idx+1 #Removed may 3rd
   #pre_VS_idx <- knn$nn.idx-1

   VS_up_idx <- ifelse(SWORD_dat$Node_id[knn$nn.idx+1] > SWORD_dat$Node_id[knn$nn.idx-1],knn$nn.idx+1,knn$nn.idx-1 )
   VS_dn_idx <- ifelse(SWORD_dat$Node_id[knn$nn.idx+1] > SWORD_dat$Node_id[knn$nn.idx-1],knn$nn.idx-1,knn$nn.idx+1 )
   dat_beam$VS_lon_dn <- SWORD_dat$Lon[VS_up_idx]
   dat_beam$VS_lat_dn <- SWORD_dat$Lat[VS_up_idx]
   dat_beam$VS_id_dn <- SWORD_dat$Node_id[VS_up_idx]

   dat_beam$VS_lon_up <- SWORD_dat$Lon[VS_dn_idx]
   dat_beam$VS_lat_up <- SWORD_dat$Lat[VS_dn_idx]
   dat_beam$VS_id_up <- SWORD_dat$Node_id[VS_dn_idx]
   dat_beam$dist_to_VS <- knn$nn.dists

   dat_beam <- dat_beam[dat_beam$type == 1 | dat_beam$type == 3 ,]
   if(dim(dat_beam)[1] < 1){
      return()
   }

   projected_dat <- project_to_centerline(dat_beam)
   dat_beam$centerline_lon <- projected_dat[,1]
   dat_beam$centerline_lat <- projected_dat[,2]
   dat_beam$IntersectAngle <- projected_dat[,3]

   return(dat_beam)
}


#' Outlier filtering for each SWORD node. Removes global and contextual outliers using quantiles and Gaussian kernel distribution.
#'
#' @param dat Main data frame containing satellite and SWORD data
#' @param myProj String containing UTM projection. Is determined based on lat and lon.
#' @return Main data frame containing outlier filtered data.
Filter_data <- function(dat,myProj){
  return_dat <- c()

  for(oo in unique(dat$Node_id)){
    dat_id <- dat[dat$Node_id == oo,]
    dat_reach <- dat[dat$Reach_id == 42264500051,]
    if(length(unique(dat_id$WaterID)) > 1){
      dist1 <- dat_id$dist_to_VS[dat_id$WaterID == unique(dat_id$WaterID)[1]]
      dist2 <- dat_id$dist_to_VS[dat_id$WaterID == unique(dat_id$WaterID)[2]]
      if(mean(dist1) > mean(dist2)){
        id <- which(dat_id$WaterID == unique(dat_id$WaterID)[2])
      } else {
        id <- which(dat_id$WaterID == unique(dat_id$WaterID)[1])
      }
      dat_id <- dat_id[id,]
    }

    if(dim(dat_id)[1] > 2 ){

      if(sd(dat_id$H_ortho) > 0.2){
        # Detect global point outliers
        GKD <- density(dat_id$H_ortho) #Gaussian kernel distribution
        WSE <- GKD$x[which.max(GKD$y)]

        quartiles <- quantile(dat_id$H_ortho, probs=c(.25, .75), na.rm = FALSE)
        IQR <- IQR(dat_id$H_ortho)

        Lower <- WSE - 1.7*IQR
        Upper <- WSE + 1.7*IQR

        dat_id <- dat_id[dat_id$H_ortho < Upper & dat_id$H_ortho > Lower,]

        if(dim(dat_id)[1] < 1){
          next()
        }
      }

      if((max(dat_id$H_ortho)-min(dat_id$H_ortho)) >  1.3 | sd(dat_id$H_ortho > 0.5)){
        next()
      }
    }#If dim > 2
    return_dat <- rbind(return_dat,dat_id)
  }#For each node id
  return(return_dat)
}



#' Iterative process that splits the data frame into sections along the river and calculates the WSS and corresponding statistics.
#'
#' @param dat Main data frame containing satellite and SWORD data
#' @param myProj String containing UTM projection. Is determined based on lat and lon.
#' @return Data frame of slopes with additional information
Calc_slope <- function(dat,myProj){
   output <- dat[1,]
   output$Slope_lon <- NA
   output$Slope_lat <- NA
   output$WSS <- NA
   output$R2 <- NA
   output$StdError <- NA
   output$pVal <- NA
   output$Nobs <- NA
   output$Resolution <- NA
   output$Beam <- NA
   output$IntersectAngle <- NA

   dat <- dat[order(dat$Reach_id),]

   XX <- 1
   Maks <- dim(dat)[1]
   count <- 0
   id_reuse <- c(1)
   use_id <- c(1)
   total_used <- c()

   while (XX <= (Maks-Min_reg_p)){
      use_id <- c(1:length(dat[,1]))

      if(length(total_used) > 0){
         use_id <- use_id[-c(total_used)]
      }

      ID_list <- c(unique(dat[XX,]$Reach_id),dat[XX,]$Rch_id_up1, dat[XX,]$Rch_id_up2, dat[XX,]$Rch_id_up3, dat[XX,]$Rch_id_up4, dat[XX,]$Rch_id_down1, dat[XX,]$Rch_id_down2, dat[XX,]$Rch_id_down3, dat[XX,]$Rch_id_down4)
      approved_ID <- which(dat$Reach_id[use_id] %in% ID_list)
      use_id <- use_id[approved_ID]

      if(length(use_id) < Min_reg_p){                          #cat("Not enough regression_dat with this reach id, XX jumps ",2*Min_reg_p,"\n")
         XX <- XX + 2*Min_reg_p
         next()
      }

      regression_dat <- dat[use_id,]

      if(length(unique(na.omit(regression_dat$WaterID))) > 1){  #cat("Contains multiple water bodies. This can remove the reach id assigned to XX \n")
         dist1 <- regression_dat$dist_to_VS[regression_dat$WaterID == unique(regression_dat$WaterID)[1]]
         dist2 <- regression_dat$dist_to_VS[regression_dat$WaterID == unique(regression_dat$WaterID)[2]]
         if(mean(dist1) > mean(dist2)){
            id <- which(regression_dat$WaterID == unique(regression_dat$WaterID)[2])
            use_id <- use_id[id]
            regression_dat <- regression_dat[id,]

            id_reject <- which(regression_dat$WaterID == unique(regression_dat$WaterID)[1])
            total_used <- c(total_used,id_reject)              #This rejects the points with waterID far from centerline
                                                               #If this is not done, we can get a slope using only a lake.
                                                               #This will be assigned to the reach id and give the impression that the slope changes
            if(length(use_id) < Min_reg_p){                    #cat("Too few points within the 8 kilometerw, XX jumps",2*Min_reg_p,"\n")
               XX <- XX+2*Min_reg_p
               next()
            }
         } else {
            id <- which(regression_dat$WaterID == unique(regression_dat$WaterID)[1])
            use_id <- use_id[id]
            regression_dat <- regression_dat[id,]

            id_reject <- which(regression_dat$WaterID == unique(regression_dat$WaterID)[2])
            total_used <- c(total_used,id_reject)

            if(length(use_id) < Min_reg_p){                    #cat("Too few points within the 8 kilometerw, XX jumps",2*Min_reg_p,"\n")
               XX <- XX+2*Min_reg_p
               next()
            }
         }
      }

      if(min(regression_dat$dist_to_VS) > 0.01){               #In case we only have accepted points in a lake
         XX <- XX+round(length(regression_dat[,1])/2)          #cat("too large dist to centerline (maybe lake ID), XX jumps 30 \n")
         next()
      }

      pos_UTM <- as.matrix(cbind(regression_dat$UTM_E,regression_dat$UTM_N))
      current_reaches <- SWORD_dat[(SWORD_dat$Reach_id < max(regression_dat$Reach_id)+20 & SWORD_dat$Reach_id>min(regression_dat$Reach_id)-20 ),]
      current_reaches <- rbind(current_reaches,SWORD_dat[which(SWORD_dat$Reach_id %in% ID_list),])
      current_reaches <- unique(current_reaches)
      dist <- CalcDist(pos_UTM,as.matrix(cbind(current_reaches$Lon_node,current_reaches$Lat_node)),as.character(myProj))
      regression_dat$river_dist <- dist[,3]-dist[1,3]

      ID <- which(regression_dat$river_dist > -Max_reg_dist)   #These are the points that I must reject as they are too far from other points
      if(length(ID) < Min_reg_p){                              #cat("Too few points within the 8 kilometerw, XX jumps",2*Min_reg_p,"\n")
         XX <- XX+2*Min_reg_p
         next()
      }

      regression_dat <- regression_dat[ID,]

      ID <- which(regression_dat$river_dist < (Max_reg_dist-abs(min(regression_dat$river_dist))) )#Accept points that cover less than Max_reg_dist
      if(length(ID) < Min_reg_p){                              #cat("Too few points within acceptable regression dist, XX jumps",2*Min_reg_p,"\n")
         XX <- XX+2*Min_reg_p
         next()
      }
      regression_dat <- regression_dat[ID,]

      if(max(regression_dat$river_dist)-min(regression_dat$river_dist) < Min_reg_dist ){ #If all dists within Max_reg_dist are less than Min_reg_dist, then all dists are less than Min_reg_dist meters apart
         XX <- XX + 20                                         #cat("Regression_dist:", max(regression_dat$river_dist)-min(regression_dat$river_dist), "XX jumps",20,\n")
         next()
      }

      LSR <- lm(regression_dat$H_ortho ~ regression_dat$river_dist)
      WSS <- coef(LSR)[2]
      LSR_Rsquared <- summary(LSR)$r.squared
      StdError <- coef(summary(LSR))[2,2]
      pVal <- coef(summary(LSR))[2,4]

      this_reach <- current_reaches[current_reaches$Reach_id == ID_list[1],]
      this_reach_min <- which(this_reach$Node_id == min(this_reach$Node_id))
      up_reach <- current_reaches[current_reaches$Reach_id == ID_list[2],]
      dn_reach <- current_reaches[current_reaches$Reach_id == ID_list[6],]

      if(dim(dn_reach)[1] > 0){
         dn_reach_max <- which(dn_reach$Node_id == max(dn_reach$Node_id))
         if(sqrt( (this_reach$Lon_node[this_reach_min] - dn_reach$Lon_node[dn_reach_max])^2 +   (this_reach$Lat_node[this_reach_min] - dn_reach$Lat_node[dn_reach_max])^2   ) > 0.05 & WSS < 0 & StdError*100*1000 < 3){
            regression_dat$river_dist <- -1*regression_dat$river_dist

            LSR <- lm(regression_dat$H_ortho ~ regression_dat$river_dist)
            WSS <- coef(LSR)[2]
            LSR_Rsquared <- summary(LSR)$r.squared
            StdError <- coef(summary(LSR))[2,2]
            pVal <- coef(summary(LSR))[2,4]

         }
      } else if (dim(up_reach)[1] > 0){
         up_reach_max <- which(up_reach$Node_id == max(up_reach$Node_id))
         if(sqrt( (this_reach$Lon_node[this_reach_min] - up_reach$Lon_node[up_reach_max])^2 +   (this_reach$Lat_node[this_reach_min] - up_reach$Lat_node[up_reach_max])^2   ) > 0.05 & WSS < 0){
            #print("   Inverting the distance")
            regression_dat$river_dist <- -1*regression_dat$river_dist

            LSR <- lm(regression_dat$H_ortho ~ regression_dat$river_dist)
            WSS <- coef(LSR)[2]
            LSR_Rsquared <- summary(LSR)$r.squared
            StdError <- coef(summary(LSR))[2,2]
            pVal <- coef(summary(LSR))[2,4]
         }
      }

      count=count+1
      output[count,] <- dat[XX,]
      output$R2[count] <- LSR_Rsquared
      output$WSS[count] <- WSS*100*1000                        #cm/km
      output$StdError[count] <- StdError*100*1000              #cm/km
      output$pVal[count] <- pVal
      output$Nobs[count] <- dim(regression_dat)[1]
      output$Beam[count] <- paste0(sort(unique(regression_dat$Beam)),collapse='')
      output$IntersectAngle[count] <- median(regression_dat$IntersectAngle)

      mean_dist <- (abs(max(regression_dat$river_dist))+abs(min(regression_dat$river_dist)))/2
      center_id <- regression_dat$Node_id[1]+ mean_dist/200*10
      IDD <- which.min(abs(current_reaches$Node_id - center_id))

      output$Slope_lon[count] <- current_reaches$Lon_node[IDD]
      output$Slope_lat[count] <- current_reaches$Lat_node[IDD]
      output$Node_id[count] <- current_reaches$Node_id[IDD]
      output$Reach_id[count] <- current_reaches$Reach_id[IDD]
      output$Resolution[count] <- max(regression_dat$river_dist)-min(regression_dat$river_dist)
      output$Reach_slope[count] <- current_reaches$Reach_slope[IDD]

      XX <- XX+length(use_id)+1
      total_used <- c(total_used,use_id)

   }
   if(!(count == 0)){
      return(output)
   } else {
      return()
   }
}



#' Projects satellite data onto river centrline based on nearby SWORD nodes.
#'
#' @param dat Main data frame containing satellite and SWORD data
#' @return Main data frame including coordinates for the projected position and the intersection angle between river and satellite
project_to_centerline <- function(dat){

   Project_angle <- c()
   Projected_lon <- c()
   Projected_lat <- c()
   for(node in unique(dat$Node_id)){
      newid <- which(dat$Node_id == node)

      newdat <- dat[dat$Node_id == node,]
      borrow_point=FALSE
      borrow_point2=FALSE
      if(dim(newdat)[1] < 2 & dim(dat)[1] > 1){ #If there is only 1 point, use another icesat point to get direction
         borrow_point = TRUE
         newdat <- rbind(newdat,newdat)
         borrow_id <- ifelse(newid < dim(dat)[1],newid+1,newid-1)

         newdat$Lon[2] <- dat$Lon[borrow_id]
         newdat$Lat[2] <- dat$Lat[borrow_id]
      } else if( dim(dat)[1] < 1){
         borrow_point2 = TRUE
         newdat <- rbind(newdat,newdat)
      }

      #Find correction to apply to the icesat points to move intersection to node
      correction <- c(newdat$Lon[1]-unique(newdat$Lon_node),newdat$Lat[1]-unique(newdat$Lat_node))
      newdat$Lon_corr <- newdat$Lon-correction[1]
      newdat$Lat_corr <- newdat$Lat-correction[2]

      #The node is now origo of the intersection
      dist_up <- sqrt( (mean(newdat$Lon)-unique(newdat$VS_lon_up))^2 + (mean(newdat$Lat)-unique(newdat$VS_lat_up))^2)
      dist_dn <- sqrt( (mean(newdat$Lon)-unique(newdat$VS_lon_dn))^2 + (mean(newdat$Lat)-unique(newdat$VS_lat_dn))^2)
      if(dist_up < dist_dn){
         b <- cbind(newdat$VS_lon_up-newdat$Lon_node,newdat$VS_lat_up-newdat$Lat_node)
      } else {
        b <- cbind(newdat$VS_lon_dn-newdat$Lon_node,newdat$VS_lat_dn-newdat$Lat_node)
      }

      Vector_dat1 <-cbind(newdat$Lon_corr[dim(newdat)[1]]-newdat$Lon_node[1],newdat$Lat_corr[dim(newdat)[1]]-newdat$Lat_node[1])
      Angle <- Project_angle(c(Vector_dat1,b[1,]))

      Vector_dat <- cbind(newdat$Lon-newdat$Lon_node,newdat$Lat-newdat$Lat_node)
      Project_coor <- apply(cbind(Vector_dat,b),1,Project)

      if(borrow_point == TRUE){
         Project_coor <- Project_coor[,1]
         Projected_lon <- c(Projected_lon,newdat$Lon_node[1] + Project_coor[1])
         Projected_lat <- c(Projected_lat,newdat$Lat_node[1] + Project_coor[2])
      } else if (borrow_point2 == TRUE){
         Angle <- NA
         Project_coor <- Project_coor[,1]
         Projected_lon <- c(Projected_lon,newdat$Lon_node[1] + Project_coor[1])
         Projected_lat <- c(Projected_lat,newdat$Lat_node[1] + Project_coor[2])
      } else {
         Projected_lon <- c(Projected_lon,newdat$Lon_node + Project_coor[1,])
         Projected_lat <- c(Projected_lat,newdat$Lat_node + Project_coor[2,])
      }

      Project_angle[newid] <- Angle
   }

   return(cbind(centerline_lon = Projected_lon, centerline_lat=Projected_lat,IntersectAngle = Project_angle))
}

#' Calculates the distance of all points along the river to get the x-dimension of the regression.
#'
#' @param dat Data coordinates in UTM
#' @param centerLineLL Coordinates of SWORD centerline in degrees (lon,lat)
#' @param myProj String containing UTM projection. Is determined based on lat and lon.
#' @return Standard data frame now including the position along the river relative to arbitrary point
CalcDist<-function(dat,centerLineLL,myProj){

   xy<-project(as.matrix(centerLineLL[,1:2]), myProj)
   xydat <- dat[,1:2]

   # calculate length of line segments of centerLine
   NP<-nrow(xy)
   BeginPx<-xy[1:(NP-1),1]
   BeginPy<-xy[1:(NP-1),2]
   begin<-cbind(BeginPx,BeginPy)
   EndPx<-xy[2:NP,1]
   EndPy<-xy[2:NP,2]
   end<-cbind(EndPx,EndPy)

   lineSegLength<-sqrt((EndPx-BeginPx)^2+(EndPy-BeginPy)^2)
   alpha<-(EndPy-BeginPy)/(EndPx-BeginPx)

   #identify line segment and get dist
   MyPos<-matrix(NA,nrow=nrow(xydat),ncol=3)
   sat <- xydat[,1:2] ################ This is my line
   MyPos[,1:2]<-sat[,1:2]

   for(i in 1:nrow(MyPos)){
      id<-which.min(sqrt((sat[i,1]-begin[,1])^2+(sat[i,2]-begin[,2])^2))

      if(abs(sat[i,1]-begin[id,1])<0.00001){
         alphap<-(sat[i,2]-end[id,2])/(sat[i,1]-end[id,1])
      } else{alphap<-(sat[i,2]-begin[id,2])/(sat[i,1]-begin[id,1])}

      test<-abs(alphap-alpha[id])<0.0000000001

      if(test){
         if (id==1) distTOT<-getdist(sat[i,],begin[id,])
         else distTOT<-sum(lineSegLength[1:(id-1)])+getdist(sat[i,],begin[id,])
      } else {
         if(id==2){
            distTOT<-getdist(sat[i,],begin[id-1,])
         } else if(id==1){
            distTOT<-getdist(sat[i,],begin[id,])
         } else {
            distTOT<-sum(lineSegLength[1:(id-1-1)])+getdist(sat[i,],begin[id-1,])
         }
      }
      MyPos[i,3]<-distTOT
   }

   return(MyPos)
}



#' Calculate distance between points
#'
#' @param v Point coordinates
#' @param w input coordinate
#' @return Distance
getdist<-function(v, w)sqrt((v[1] - w[1])^2 + (v[2] - w[2])^2)

#' Projects satellite data onto the centerline
#'
#' @param dat Main data frame
#' @return Reprojected data
Project <- function(dat){
   as.vector(dat[1:2]%*%dat[3:4]) / (as.vector(dat[3:4] %*% dat[3:4])) *as.vector(dat[3:4])

}


#' Calculates angle between satellite beam and local river reach
#'
#' @param dat Main data frame
#' @return Intersection angle
Project_angle <- function(dat){

   v5 <- atan2(dat[2],dat[1])*180/pi
   v6 <- atan2(dat[4],dat[3])*180/pi

   v7 <- v5-v6
   v7 <- abs(v7)

   v7 <- ifelse(v7 > 180,360-v7,v7)
   v7 <- ifelse(v7 > 90,180-v7,v7)

   return(v7)
}


#' Loads SWORD centerlines and nodes as a text file. If this doesn't exist, the file is created based on the nc SWORD file
#'
#' @param SWORD_dir String. Path to folder with SWORD data either as text file or as nc
#' @param Version String. String containing SWORD version to use
#' @param Area String. String containing SWORD area to use
#' @return Nothing. Data is loaded
handle_SWORD <- function(SWORD_dir, Version, Area){

   if(exists("SWORD_dat")){#If data is already loaded into R

   } else { #If txt file already exist, create path name
      ready_SWORD_data <- paste(SWORD_dir,"/Processed_SWORD_",Version,"_",Area,".txt",sep="" )

      if(!file.exists(ready_SWORD_data)){ #If file does not exist, create
         start_time <- Sys.time()
         Produce_SWORD_data(SWORD_dir,Area,Version)
         end_time <- Sys.time()
         print(end_time - start_time)

      }
      #Load SWORD text file
      start_time <- Sys.time()
      SWORD_dat <<- read.table(ready_SWORD_data,sep=",",header=TRUE)
      end_time <- Sys.time()
      print(end_time - start_time)
   }
}

#'
#' #' Checks data formats and give warnings if data is not as expected. Assigns column names to data frame
#' #'
#' #' @param dat Main data frame containing satellite and SWORD data
#' #' @return Data frame with columns names
#' Check_data_format <- function(dat,Filename){
#'
#'    ##### Check data formats
#'    browser()
#'    test <- try((dim(dat)[2] == 6 & is.na(Occ_thr))| (dim(dat)[2] == 7 & !is.na(Occ_thr)),silent=TRUE)
#'
#'    if(test ==FALSE){
#'
#'    }
#'    if(class(test) == "try-error"){
#'       cat(substr(as.character(Sys.time()),1,16), "Load input data failed for file",filelist[file2],".\n",test[1],"\n")
#'       return()
#'    }
#'
#'    if(dim(dat)[2] < 6 | (dim(dat)[2] < 7 & !is.na(Occ_thr))){
#'       cat(substr(as.character(Sys.time()),1,16),"Load input data failed for file",Filename,". Wrong number of columns \n")
#'       browser()
#'       return()
#'    }
#'
#'    if(dim(dat)[2] == 7 & !is.na(Occ_thr)){ #If occurrence is included and we are happy
#'       dat <- dat[dat[,7] > Occ_thr,1:7]
#'       colnames(dat) <- c("DecYear", "Lat", "Lon","H_ortho","WaterID","Beam","Occurrence")
#'       if(any(dat$Occurrence < 0)) {warning(paste(as.character(Sys.time()),"Warning: Water occurrence is less than 0. Check input data."))}
#'       if(any(dat$Occurrence > 100)) {warning(paste(as.character(Sys.time()),"Warning: Water occurrence is larger than 100. Check input data."))}
#'
#'    } else if(dim(dat)[2] == 7 & is.na(Occ_thr)) { #If dimensions of data does not correspond to 'no occurrence'
#'       stop(paste("Error: Data contains 7 columns but no occurrence threshold has been defined."))
#'
#'    } else if (!(dim(dat)[2] == 6) & !(dim(dat)[2] == 7) ){ #If dimensions of data is completely wrong
#'       warning(paste("Warning: Data must contain 6 or 7 columns (water occurrence optional). Data has",dim(dat)[2],"columns."))
#'
#'    } else { #If occurrence is not included and we are happy
#'       colnames(dat) <- c("DecYear", "Lat", "Lon","H_ortho","WaterID","Beam") #Removed "Occurrence" in 2nd of may.
#'    }
#'
#'
#'    if(max(dat$DecYear)- min(dat$DecYear) > 0.00025){
#'       warning(paste(as.character(Sys.time()),"Warning: Data may not contain multiple passes. Time duration covered in file is more than 2 hours."))
#'    }
#'
#'    if(any(dat$DecYear < 0)) {warning(paste(as.character(Sys.time()),"Warning: Decimal year contains negative values. Ensure that input data are correct."))}
#'    if(any(abs(dat$Lat) > 90)) {warning(paste(as.character(Sys.time()),"Warning: Latitude contains values larger than +/-90. Ensure that input data are correct."))}
#'    if(any(abs(dat$Lon) > 360)) {warning(paste(as.character(Sys.time()),"Warning: Longitude contains values larger than 360. Ensure that input data are correct."))}
#'    #if(any(!(nchar(dat$WaterID) == 6 |nchar(dat$WaterID) == 7) )) {warning(paste(as.character(Sys.time()),"Warning: WaterID does not contain 6 or 7 digits. Ensure that input data are correct."))}
#'    if(any(dat$Beam > 6 | dat$Beam < 1)) {warning(paste(as.character(Sys.time()),"Warning: Beam number is not value between 1 and 6. Ensure that input data are correct."))}
#'
#'    return(dat)
#' }
