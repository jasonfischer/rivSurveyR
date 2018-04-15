#' xSec.secondary
#' 
#' Compiles ADCP cell data exported from SonTek RiverSurveyor as MATLAB file into an object of class "adcp.secondary" which inherits from data.table
#' @param data	A list of MATLAB files exported from RiverSurveyor to be compiled
#' @param transectNames	A list of transect names corresponding to the transects represented by the MATLAB files. See 'details' for additional information.
#' @param depthReference	Defines the depth measurements to be used. The default ("unit") uses the method defined in RiverSurveyor, but \code{depthReference} can be set to use measurements from the vertical beam ("VB"), bottom-track beams ("BT"), or both ("composite").
#' @param project	Logical.  If TRUE (default), transect replicates are projected to a mean transect line using an orthogonal projection of the x and y coordinates. If FALSE transects are not projected to a mean transect line.
#' @details	When assigning transectNames, order matters, the transect names must be in the same order of the MATLAB files assigned to the data argument. Ex., if the first two files in the list are Transect1a and Transect1b, which are replicates of Transect1, t first two files in the name list are "Transect1" and "Transect1".  The processing functions will spatially average and compile all data from the same transects, this tells the processing functions to average and compile these two files. Listing the names as "Transect1a" and "Transect1b" results in the transects being processed separately.  \code{transectNames} will only expect characters and names must be enclosed in "".
#' @return	A \link[data.table]{data.table} with the transect name (transectName), cellID, depth averaged velocity to the east (Mean.Vel.E), depth averaged velocity to the north (Mean.Vel.N), depth, east velocity within a cell (Vel.E), north velocity within a cell (Vel.N), vertical velocity within a cell (Vel.up), error velocity within a cell (Vel.error), blanking distance or depth where measurements first begin (StartDepth), cell height, cell number within an ensemble (cell one is closest to the surface; cellNumber), cell depth, UTM x coordinates (UTM_X), UTM y coordinates (UTM_Y), longitude, latitude, altitude of water surface, and distance from starting point. If project = TRUE, orthogonally projected longitude, latitude, and UTM coordinates (Longitude_Proj, Latitude_Proj, UTM_X_Proj, UTM_Y_Proj) are also returned.
#' @export
#' @examples
#' data(mNine)
#' #mNine is a list of MATLAB files
#' names(mNine)
#' #Drop the last letter (l or r) from the MATLAB file names
#' tNames <- substr(names(mNine), 0, nchar(names(mNine))-1)
#' xSec.secondary(mNine, tNames)
#' @seealso	\link{project.transect}
 
xSec.secondary <- function(data, depthReference = "unit", transectNames, project = TRUE) {
  if(missing(transectNames) == TRUE){
    transectNames <- as.list(as.character(seq_along(data)))
  } else {
    transectNames <- as.list(transectNames)
  }
  stopifnot(length(data) == length(transectNames), is.list(data) == TRUE)
  #Create empty list to populate with velocity and coordinate data  
  secondaryFlow.1 <- list()
  #iteratively populate list with data.tables of data from each transect replicate
  for (i in seq_along(data)){
    #RiverSurveyor version 3.9.50 includes site info. in the MATLAB export, to make exports from this version compatiable with older versions, remove the site info.
    if("SiteInfo" %in% names(data[[i]])){
      data[[i]] <- data[[i]][-2]
    } else {
    }
    #Check that coordinate system of each transect is ENU
    if(!data[[i]][[1]][4,,1][[1]]==2){
      stop("Coordinate System of transect ", i, " in data list is not ENU. Coordinate system can be changed to ENU in RiverSurveyor Live.")
    } else {      
    }
    
    transect <- data[[i]]
    meanVel <- data.frame(transect[[6]][10,,1][[1]])
    names(meanVel) <- c("Mean.Vel.1", "Mean.Vel.2")
    #pull out depth data measured with the vertical beam, bottom track mode, or both, based on user defined depthReference
    if(depthReference == "unit"){
      if(transect[[1]][7,,1][[1]] == 0){
        depth <- transect[[7]][1,,1][[1]]
      } else if(transect[[1]][7,,1][[1]] == 1){
        depth <- transect[[7]][2,,1][[1]]
      } else {
        depth <- transect[[6]][7,,1][[1]]
      }
    } else if(depthReference == "VB"){
      depth <- transect[[7]][1,,1][[1]]
    } else if(depthReference == "BT"){
      depth <- transect[[7]][2,,1][[1]]
    } else if(depthReference == "composite"){
      depth <- transect[[6]][7,,1][[1]]
    } else {
      stop("depthReference must either be defined as unit, VB, BT, or composite")
    }
    
    Vel.E <- transect[[9]][1,,1][[1]][,1,]
    Vel.N <- transect[[9]][1,,1][[1]][,2,]
    Vel.up <- transect[[9]][1,,1][[1]][,3,]
    Vel.error <- transect[[9]][1,,1][[1]][,4,]
    
    startDepth <- transect[[5]][7,,1][[1]]
    startDepth <- t(startDepth)
    cellHeight <- transect[[5]][8,,1][[1]]
    cellHeight <- t(cellHeight)
    distance <- sqrt(transect[[6]][[9]][,1]^2 + transect[[6]][[9]][,2]^2)
    distance <- t(distance)
    sampleNum <- transect[[5]][[5]]
    sampleNum <- t(sampleNum)
    
    #create a matrix containing the number of each cell, each column represents an ensemble and each row a cell
    cellNumber <- matrix(rep(seq_along(transect[[9]][1,,1][[1]][,1,1]), times = length(transect[[9]][1,,1][[1]][1,1,])), nrow = length(transect[[9]][1,,1][[1]][,1,1]), ncol = length(transect[[9]][1,,1][[1]][1,1,]))
    #create a matrix to populate with depths of each measured cell
    cellDepth <- matrix(0, nrow = length(transect[[9]][1,,1][[1]][,1,1]), ncol = length(transect[[9]][1,,1][[1]][1,1,]))
    for (w in seq_along(cellNumber[,1])){
      cellDepth[w,] <- (cellNumber[w,] * cellHeight) + startDepth
    }
    
    cellNumber <- as.vector(cellNumber)
    cellDepth <- as.vector(cellDepth)
    
    cells <- t(transect[[6]][8,,1][[1]])
    firstMeasCell <- apply(Vel.E, 2, FUN = function(x) {
      rNames <- seq_along(x)
      x <- na.omit(data.frame(x, rNames))
      if (nrow(x) == 0) {
        return(NA)
      } else {
        return(min(x$rNames))
      }
    }
    )
    firstMeasCellDep <- (firstMeasCell*cellHeight)+startDepth
    bottomCellDepth <- cells * cellHeight + firstMeasCellDep
    
    gps <- data.frame(transect[[8]][1:2,,1],transect[[8]][7,,1],transect[[8]][,,1]$UTM) #transect[[8]][,,1]$UTM is used, because UTM was pushed back a position in riversurveyor 3.9.50, this allows UTM to be called no matter what the position
    suppressWarnings(gps[which(transect[[8]][[4]] == 0),] <- NA) #if the GPS has acquired zero satellites, convert coordinates to NA
    names(gps)[4:5] <- c("UTM_X", "UTM_Y")
    
    #expand objects with one measurement per ensemble to match individual cells within each ensemble
    Mean.Vel.E <- rep(meanVel$Mean.Vel.1, each=nrow(Vel.E))
    Mean.Vel.N <- rep(meanVel$Mean.Vel.2, each=nrow(Vel.E))
    depth <-rep(depth, each=nrow(Vel.E))
    UTM_X <- rep(gps$UTM_X, each=nrow(Vel.E))
    UTM_Y <- rep(gps$UTM_Y, each=nrow(Vel.E))
    Longitude <- rep(gps$Longitude, each=nrow(Vel.E))
    Latitude <- rep(gps$Latitude, each=nrow(Vel.E))
    Altitude <- rep(gps$Altitude, each=nrow(Vel.E))
    startDepth <- rep(startDepth, each=nrow(Vel.E))
    cellHeight <- rep(cellHeight, each=nrow(Vel.E))  
    bottomCellDepth <- rep(bottomCellDepth, each=nrow(Vel.E))  
    distance <- rep(distance, each=nrow(Vel.E))
    sampleNum <- rep(sampleNum, each=nrow(Vel.E))
    
    Vel.E <- as.vector(Vel.E)
    Vel.N <- as.vector(Vel.N)
    Vel.up <- as.vector(Vel.up)
    Vel.error <- as.vector(Vel.error)
    
    secondaryFlow.1[[i]] <- data.table(Mean.Vel.E, Mean.Vel.N, depth, Vel.E, Vel.N, Vel.up, Vel.error, startDepth, cellHeight, cellNumber, cellDepth, bottomCellDepth, UTM_X, UTM_Y, Longitude, Latitude, Altitude, distance, sampleNum)
    secondaryFlow.1[[i]]$transectName <- transectNames[[i]]
  }
  
  #row bind ADCP data compiled from each transect into a single data.table
  secondaryFlow <- do.call("rbind", secondaryFlow.1)
  
  secondaryFlow <- secondaryFlow[cellDepth<=bottomCellDepth | is.na(bottomCellDepth) ,,]
  secondaryFlow <- secondaryFlow[!is.na(bottomCellDepth) | !cellNumber>1 ,,]
  secondaryFlow$depth[secondaryFlow$depth<=0] <- NA 
  
  secondaryFlow$cellID <- seq_along(secondaryFlow[[1]])
  
  #if project is TRUE, spatial coordinates are passed to project.transect
  if(project == TRUE){
    utmProj <- project.transect(x = secondaryFlow$UTM_X, y = secondaryFlow$UTM_Y, transectNames = secondaryFlow$transectName)
    setnames(utmProj, c("xProj", "yProj"), c("UTM_X_Proj", "UTM_Y_Proj"))
    latLongProj <- project.transect(x = secondaryFlow$Longitude, y = secondaryFlow$Latitude, transectNames = secondaryFlow$transectName)
    setnames(latLongProj, c("xProj", "yProj"), c("Longitude_Proj", "Latitude_Proj"))
    xyProj <- data.table(transectName = latLongProj$transectNames, cellID = latLongProj$id, Longitude_Proj = latLongProj$Longitude_Proj, Latitude_Proj = latLongProj$Latitude_Proj, UTM_X_Proj = utmProj$UTM_X_Proj, UTM_Y_Proj = utmProj$UTM_Y_Proj)
    
    keycols = c("transectName", "cellID")
    setkeyv(secondaryFlow, keycols)
    setkeyv(xyProj, keycols)
    secondaryFlow <- merge(secondaryFlow, xyProj)
  } else {
  }
  
  class(secondaryFlow) <- c(class(secondaryFlow), "adcp.secondary")
  
  return(secondaryFlow)  
}