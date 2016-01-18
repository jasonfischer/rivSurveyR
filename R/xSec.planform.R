#' xSec.planform
#' 
#' Compiles planform ADCP data exported from SonTek RiverSurveyor as MATLAB file into an object of class "adcp.planform" which inherits from data.table
#' @param data	A list of MATLAB files exported from RiverSurveyor to be compiled
#' @param transectNames	A list of transect names corresponding to the transects represented by the MATLAB files. See 'details' for additional information.
#' @param depthReference	Defines the depth measurements to be used. The default ("unit") uses the method defined in RiverSurveyor, but \code{depthReference} can be set to use measurements from the vertical beam ("VB"), bottom-track beams ("BT"), or both ("composite").
#' @param layerReference	Defines layerAvg bounds x1 and x2 as a distance from surface ("surface", default), distance above bottom ("bottom"), or percent of total depth ("percent"). If \code{layerReference} is "bottom" x1 and x2 are reversed.
#' @param x1	Lower bound velocities are averaged over. Either a numeric representing a depth or percentage or a character where "minimum" equals the minimum distance from the reference point (i.e., 0 when \code{layerReference} = "surface" and bottomCellDepth when \code{layerReference} = "bottom") and "maximum" equals maximum distance from the reference point. Must be numeric if \code{layerReference} = "percent". Defaults to "minimum".
#' @param x2	Upper bound velocities are averaged over. Either a numeric representing a depth or percentage or a character where "minimum" equals the minimum distance from the reference point (i.e., 0 when \code{layerReference} = "surface" and bottomCellDepth when \code{layerReference} = "bottom") and "maximum" equals maximum distance from the reference point. Must be numeric if \code{layerReference} = "percent". Defaults to "maximum".
#' @param project	Logical.  If TRUE (default), transect replicates are projected to a mean transect line using an orthogonal projection of the x and y coordinates. If FALSE transects are not projected to a mean transect line.
#' @details	When assigning transectNames, order matters, the transect names must be in the same order of the MATLAB files assigned to the data argument. Ex., if the first two files in the list are Transect1a and Transect1b, which are replicates of Transect1, t first two files in the name list are "Transect1" and "Transect1".  The processing functions will spatially average and compile all data from the same transects, this tells the processing functions to average and compile these two files. Listing the names as "Transect1a" and "Transect1b" results in the transects being processed separately.  \code{transectNames} will only expect characters and names must be enclosed in "".
#' @return	A \link[data.table]{data.table} with the transect name (transectName), UTM x coordinates (UTM_X), UTM y coordinates (UTM_Y), depth, depth averaged velocity to the east (Mean.Vel.E), depth averaged velocity to the north (Mean.Vel.N), depth of the bottom most cell in an ensemble (bottomCellDepth), distance of the bottom cell to the bottom (distToBottom), bottom cell velocity to the east (BC.Vel.E), bottom cell velocity to the north (BC.Vel.N), bottom cell vertical velocity (BC.Vel.Up), bottom cell error velocity (BC.Vel.Error), layer averaged velocity to the east (layerVel.E), layer averaged velocity to the north (layerVel.N), layer averaged vertical velocity (layerVel.Up), longitude, latitude, altitude of water surface, distance from starting point, and temperature. If project = TRUE, orthogonally projected longitude, latitude, and UTM coordinates (Longitude_Proj, Latitude_Proj, UTM_X_Proj, UTM_Y_Proj) are also returned.
#' @export
#' @examples
#' data(mNine)
#' #mNine is a list of MATLAB files
#' names(mNine)
#' #Drop the last letter (l or r) from the MATLAB file names
#' tNames <- substr(names(mNine), 0, nchar(names(mNine))-1)
#' xSec.planform(mNine, tNames)
#' #Get velocities within a meter of the river bottom
#' xSec.planform(mNine, tNames, layerReference = "bottom", x1 = 1, x2 = "minimum")
#' @seealso	\link{layerAvg} and \link{project.transect}
 
xSec.planform <- function(data, transectNames, depthReference = "unit", layerReference = "surface", x1 = "minimum", x2 = "maximum", project = TRUE){
  if(missing(transectNames) == TRUE){
    transectNames <- as.list(as.character(seq_along(data)))
  } else {
    transectNames <- as.list(transectNames)
  }
  stopifnot(length(data) == length(transectNames), is.list(data) == TRUE)
  
  #create list to populate data from each transect with
  ADCP.1 <- list()
  
  #Check that coordinate system of each transect is ENU
  for (i in seq_along(data)){
    if(!data[[i]][[1]][4,,1][[1]]==2){
      stop("Coordinate System of transect ", i, " in data list is not ENU. Coordinate system can be changed to ENU in RiverSurveyor Live.")
    } else {      
    }
    
    #depth, velocity, cell info, and GPS data are iteratively extracted from original lists from MATLAB files and appended to ADCP.1
    transect <- data[[i]]    
    meanVel <- data.frame(transect[[6]][10,,1][[1]])
    names(meanVel) <- c("Mean.Vel.E", "Mean.Vel.N")
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
      stop("depthReference must either defined as unit, VB, BT, or composite")
    }
    gps <- data.frame(transect[[8]][1:2,,1][1:2], transect[[8]][7,,1][1],transect[[8]][11,,1][1])
    names(gps)[4:5] <- c("UTM_X", "UTM_Y")
    suppressWarnings(gps[which(transect[[8]][[4]] == 0),] <- NA) #if the GPS has acquired zero satellites, convert coordinates to NA, warnings are suppressed, because datasets with complete satellite coverage will trigger a warning about no non-missing arguments
    satellites <- transect[[8]][[4]]
    vel <- transect[[9]][1,,1][[1]]
    cells <- transect[[6]][8,,1][[1]]
    distance <- sqrt(transect[[6]][[9]][,1]^2 + transect[[6]][[9]][,2]^2)
    
    cellHeight <- transect[[5]][8,,1][[1]]
    startDepth <- transect[[5]][7,,1][[1]]
    temperature <- transect[[5]][[2]]
    bottomCellDepth <- cells * cellHeight + startDepth
    distToBottom <- depth - bottomCellDepth
    names(distToBottom) <- c("distToBottom")
    #set-up a matrix to hold east, north, up, and difference velocities of the deepest cell measured for each ensemble
    bottomCellVel <- matrix(ncol=4, nrow=nrow(cells))
    #iteratively populate the matrix with velocities from the deepest cell measured in each ensemble
    for (w in seq_along(cells)){
      if(cells[w,1] == 0){
        bottomCellVel[w,] <- NaN
      } else {
        bottomCellVel[w,] <- vel[cells[w,1],,w]
      }
      
    }
    bottomCellVel <- data.frame(bottomCellVel)
    names(bottomCellVel) <- c("BC.Vel.E" , "BC.Vel.N" , "BC.Vel.Up", "BC.Vel.Error")
    stopifnot(layerReference %in% c("percent", "bottom", "surface"))
    #determine the value of x1 and x2 to be passed to layerAvg, based on layerReference and distance from either the surface or bottom
    if(layerReference == "percent"){
      x1.1 <- depth * x1/100
      x2.1 <- depth * x2/100
    }else if(layerReference == "bottom"){
      if(x2 == "minimum"){
        x2.1 <- bottomCellDepth
      }else if(x2 == "maximum"){
        x2.1 <- cellHeight + startDepth
      }else{
        x2.1 <- depth - x2
      }
      if(x1 == "minimum"){
        x1.1 <- bottomCellDepth
      }else if(x1 == "maximum"){
        x1.1 <- matrix(rep(0, times = length(depth)))
      }else{
        x1.1 <- depth - x1
      }
    }else{ 
      if(x1 == "minimum"){
        x1.1 <- matrix(rep(0, times = length(depth)))
      }else if(x1 == "maximum"){
        x1.1 <- bottomCellDepth
      }else{
        x1.1 <- matrix(rep(x1, times = length(depth)))
      }
      if(x2 == "minimum"){
        x2.1 <- cellHeight + startDepth
      }else if(x2 == "maximum"){
        x2.1 <- bottomCellDepth
      }else{
        x2.1 <- matrix(rep(x2, times = length(depth)))
      }
    }
    #create a matrix containing the number of each cell, each column represents a ensemble and each row a cell
    cellNumber <- matrix(rep(seq_along(vel[,1,1]), times = length(vel[1,1,])), nrow = length(vel[,1,1]), ncol = length(vel[1,1,]))
    #create a matrix to populate with depths of each measured cell
    cellDepth <- matrix(0, nrow = length(vel[,1,1]), ncol = length(vel[1,1,]))
    for (w in seq_along(cellNumber[,1])) {
      cellDepth[w,] <- (cellNumber[w,] * cellHeight) + startDepth
    }
    
    #create list to populate with layer averaged velocities and calculate layer averaged velocity for each ensemble
    layAvg <- list()
    for(w in seq_along(cellDepth[1,])) {
      if(depth[w,] <= 0 | is.na(depth[w,]) == TRUE) {
        layAvg[w] <- NA
      } else {
        layerVel.E <- layerAvg(x = cellDepth[,w], y = vel[,1,w], x1 = x1.1[w,], x2 = x2.1[w,])
        layerVel.N <- layerAvg(x = cellDepth[,w], y = vel[,2,w], x1 = x1.1[w,], x2 = x2.1[w,])
        layerVel.Up <- layerAvg(x = cellDepth[,w], y = vel[,3,w], x1 = x1.1[w,], x2 = x2.1[w,])
        
        layAvg[[w]] <- data.frame(layerVel.E, layerVel.N, layerVel.Up)
      }
    }
    #row bind all layer average velocities together
    lav <- do.call("rbind", layAvg)
    ADCP.1[[i]] <- data.frame(depth, meanVel, bottomCellDepth, distToBottom, bottomCellVel, lav, gps, distance, temperature)
    ADCP.1[[i]]$transectName <- transectNames[[i]]    
  }
  
  #row bind ADCP data compiled from each transect into a single data.table
  ADCP <- do.call("rbind", ADCP.1)
  ADCP <- data.table(ADCP)
  
  ADCP$cellID <- seq_along(ADCP[[1]])
  
  #if project is TRUE, spatial coordinates are passed to project.transect
  if(project == TRUE){
    utmProj <- project.transect(x = ADCP$UTM_X, y = ADCP$UTM_Y, transectNames = ADCP$transectName)
    setnames(utmProj, c("xProj", "yProj"), c("UTM_X_Proj", "UTM_Y_Proj"))
    latLongProj <- project.transect(x = ADCP$Longitude, y = ADCP$Latitude, transectNames = ADCP$transectName)
    setnames(latLongProj, c("xProj", "yProj", "id"), c("Longitude_Proj", "Latitude_Proj", "cellID"))
    xyProj <- data.table(transectName = latLongProj$transectNames, cellID = latLongProj$cellID, Longitude_Proj = latLongProj$Longitude_Proj, Latitude_Proj = latLongProj$Latitude_Proj, UTM_X_Proj = utmProj$UTM_X_Proj, UTM_Y_Proj = utmProj$UTM_Y_Proj)
    
    keycols = c("transectName", "cellID")
    setkeyv(ADCP, keycols)
    setkeyv(xyProj, keycols) 
    ADCP <- merge(ADCP, xyProj)
  } else {    
    keycols = c("transectName", "cellID")
    setkeyv(ADCP, keycols)    
  }
  
  ADCP <- ADCP[,!"cellID", with = FALSE]
  
  class(ADCP) <- c(class(ADCP), "adcp.planform")
  return(ADCP)
}