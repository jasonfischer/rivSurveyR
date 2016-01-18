#' process.secondary
#' 
#' Summarizes ADCP cell data exported from SonTek RiverSurveyor, calculates primary and secondary velocities, and conducts spatial averaging.
#' @param data	A list of MATLAB files exported from RiverSurveyor to be compiled
#' @param transectNames	A list of transect names corresponding to the transects represented by the MATLAB files. See 'details' for additional information.
#' @param depthReference	Defines the depth measurements to be used. The default ("unit") uses the method defined in RiverSurveyor, but \code{depthReference} can be set to use measurements from the vertical beam ("VB"), bottom-track beams ("BT"), or both ("composite").
#' @param project	Logical.  If TRUE (default), transect replicates are projected to a mean transect line using an orthogonal projection of the x and y coordinates. If FALSE transects are not projected to a mean transect line.
#' @param FUN	function used to summarize ADCP data, default is \code{mean}
#' @param binWidth	Width used to bin data along a transect, for summarization. Default is the mean distance between samples.
#' @param binHeight	Depth interval used to bin data within an ensemble, for summarization. If missing, cell depth is used group cells for summarization.
#' @param rotation	Character string indicating which rotation method is used to calculate secondary velocities. Can be no rotation ("xSec", the default), Rozovskii rotation ("rozovskii"), and/or zero secondary discharge rotation ("zeroSecQ"). If multiple rotations are used, use \code{c} to combine multiple character strings.
#' @param xWindow	Number of neighbors to include in the x orientation of the moving window, including the observation the window centers on. If NULL (the default), spatial averaging is not conducted.
#' @param yWindow	Optional. Number of neighbors to include in the y orientation of the moving window, including the observation the window centers on.
#' @param spatialCoords	A vector of column name(s) containing spatial coordinates within the dataset. The x coordinates must be designated first, followed by the y coordinates. If \code{yWindow} is not provided, the y coordinates can be omitted, otherwise they will be ignored. If \code{spatialCoords} is omitted, a moving window will be used, but data will not be spatially referenced prior to calculations.
#' @param groups	Optional. Names of column(s) containing information used to subset the data into separate groups for processing. Useful if dataset contains multiple discrete study areas.
#' @param na.ignore	Logical. If FALSE (default) NA values are assumed to be the average of non-NA values within the moving window, if TRUE NA values remain NA.
#' @param spatialAvg	Logical. If TRUE (default) and xWindow is not NULL, spatial averaging is conducted, else if FALSE no spatial averaging is conducted.
#' @param spatialStdDev	Logical. If TRUE and xWindow is not NULL, spatial standard deviation is calculated, else if FALSE (default) no spatial standard deviation is calculated.
#' @details	When assigning transectNames, order matters, the transect names must be in the same order of the MATLAB files assigned to the data argument. Ex., if the first two files in the list are Transect1a and Transect1b, which are replicates of Transect1, t first two files in the name list are "Transect1" and "Transect1".  The processing functions will spatially average and compile all data from the same transects, this tells the processing functions to average and compile these two files. Listing the names as "Transect1a" and "Transect1b" results in the transects being processed separately.  \code{transectNames} will only expect characters and names must be enclosed in "".
#' @return	A \link[data.table]{data.table} with cell depth, distance from the right bank (tDist), transect name, cell id, depth averaged velocity to the east (Mean.Vel.E), depth averaged velocity to the north (Mean.Vel.N), depth, cell velocity to the east (Vel.E), cell velocity to the north (Vel.N), cell vertical velocity (Vel.up), cell error velocity (Vel.error), blanking distance or depth where measurements first begin (StartDepth), cell height, cell number within an ensemble (cell one is closest to the surface; cellNumber), UTM x coordinates (UTM_X), UTM y coordinates (UTM_Y), longitude, latitude, cell altitude, heading of transect (assumes transect starts at river right), cellWidth, speed within a cell (speed), flow heading within a cell (velHeading), and output from rotation method(s) used. If project = TRUE, orthogonally projected longitude, latitude, and UTM coordinates (Longitude_Proj, Latitude_Proj, UTM_X_Proj, UTM_Y_Proj) are also returned.
#' @export
#' @examples
#' data(mNine)
#' #mNine is a list of MATLAB files
#' names(mNine)
#' #Drop the last letter (l or r) from the MATLAB file names
#' tNames <- substr(names(mNine), 0, nchar(names(mNine))-1)
#' adcp.secondary <- process.secondary(mNine, tNames, rotation = c("xSec", "rozovskii", "zeroSecQ"), xWindow = 21, yWindow = 5)
#' @seealso	\link{average.secondary}, \link{xSec}, \link{rozovskii}, \link{zeroSecQ}, \link{spatialAverage} and \link{spatialSD}.

process.secondary <- function(data, transectNames, depthReference = "unit", project = TRUE, FUN = "mean", binWidth, binHeight, 
                              rotation = c("xSec"), 
                              xWindow = NULL, yWindow, spatialCoords = c("tDist","cellDepth"),
                              groups = "transectName", na.ignore = FALSE, spatialAvg = "TRUE", spatialStdDev = "FALSE") {
  message("Compiling and Averaging Data")
  flush.console()
  #pass data to average planform
  ADCP <- average.secondary(data = data, transectNames = transectNames, depthReference = depthReference, project = project, binWidth = binWidth, binHeight = binHeight, FUN = FUN)
  
  ADCP <- ADCP[order(ADCP$transectName, ADCP$tDist, ADCP$cellDepth),]
  
  if(as.character(substitute(FUN)) == "mean"){
    #if xSec is specified in rotation, iteratively calcualte vx and vy for each transect    
    if("xSec" %in% rotation){
      uniqueTransects <- unique(ADCP$transectName)
      xSec.Vel <- list()
      for(i in seq_along(uniqueTransects)){
        dt <- ADCP[transectName == uniqueTransects[i],]
        xSec.Vel[[i]] <- xSec(ve = dt$Vel.E, vn = dt$Vel.N, vu = dt$Vel.up, transectHeading = mean(dt$transectHeading, na.rm = TRUE))
        xSec.Vel[[i]] <- data.table(dt[, list(transectName, tDist, cellDepth),], xSec.Vel[[i]])
      }
      xSec.Vel <- data.table(do.call("rbind", xSec.Vel))
      xSec.Vel <- xSec.Vel[order(xSec.Vel$transectName, xSec.Vel$tDist, xSec.Vel$cellDepth),]
      stopifnot(all.equal(ADCP[, list(transectName, tDist, cellDepth),], xSec.Vel[, list(transectName, tDist, cellDepth),], check.attributes=FALSE))
      ADCP <- cbind(ADCP, xSec.Vel[, !c("transectName", "tDist", "cellDepth"), with = FALSE])
    } else{
      
    }
    
    #if zeroSecQ is specified in rotation, iteratively calculate primary and secondary velocities for each transect with the zero secondary discharge method
    if("zeroSecQ" %in% rotation){
      uniqueTransects <- unique(ADCP$transectName)
      zsq <- list()
      for (i in seq_along(uniqueTransects)){
        dt <- ADCP[transectName == uniqueTransects[i],]
        zsq[[i]] <- zeroSecQ(ve = dt$Vel.E, vn = dt$Vel.N, vu = dt$Vel.up, cellHeight = dt$cellHeight, cellWidth = dt$cellWidth, transectHeading = mean(dt$transectHeading, na.rm = TRUE))
        zsq[[i]] <- data.table(dt[, list(transectName, tDist, cellDepth)], zsq[[i]])
      }
      zsq <- data.table(do.call("rbind", zsq))
      zsq <- zsq[order(zsq$transectName, zsq$tDist, zsq$cellDepth),]
      stopifnot(all.equal(ADCP[, list(transectName, tDist, cellDepth)], zsq[, list(transectName, tDist, cellDepth)], check.attributes=FALSE))
      ADCP <- cbind(ADCP, zsq[, !c("transectName", "tDist", "cellDepth"), with = FALSE])
    } else {
      
    }
    
    #if Rozovskii is specified in rotation, iteratively calculate primary and secondary velocities for each transect with the Rozovskii method
    if("rozovskii" %in% rotation){
      uniqueTransects <- unique(ADCP$transectName)
      roz <- list()
      for(i in seq_along(uniqueTransects)){
        dt <- ADCP[transectName == uniqueTransects[i],]
        roz[[i]] <- rozovskii(ve = dt$Vel.E, vn = dt$Vel.N, vu = dt$Vel.up, meanVe = dt$Mean.Vel.E, meanVn = dt$Mean.Vel.N)
        roz[[i]] <- data.table(dt[, list(transectName, tDist, cellDepth),], roz[[i]])
      }
      roz <- data.table(do.call("rbind", roz))
      roz <- roz[order(roz$transectName, roz$tDist, roz$cellDepth),]
      stopifnot(all.equal(ADCP[, list(transectName, tDist, cellDepth)], roz[, list(transectName, tDist, cellDepth)], check.attributes=FALSE))
      ADCP <- cbind(ADCP, roz[, !c("transectName", "tDist", "cellDepth"), with = FALSE])    
    } else {
      
    }    
    
    if(FALSE %in% c(rotation %in% c("xSec", "rozovskii", "zeroSecQ"))){
      warning(print(rotation[which(!rotation %in% c("xSec", "rozovskii", "zeroSecQ"))]), " rotation method not found and was ignored")
    } else if(is.null(rotation) == TRUE){
      message("Rotation = NULL, velocities were not rotated")
    } else {
      
    }
    
    if(is.null(xWindow) == TRUE | spatialAvg == FALSE & spatialStdDev == FALSE){
      
    }else{
      if(spatialAvg == TRUE & spatialStdDev == FALSE){
	#pass data to spatialAverage
        ADCP <- spatialAverage(data = ADCP, xWindow = xWindow, yWindow = yWindow, spatialCoords = spatialCoords, groups = groups, na.ignore = na.ignore)
      } else if(spatialAvg == FALSE & spatialStdDev == TRUE){
        #pass data to spatialSD
        ADCP <- spatialSD(data = ADCP, xWindow = xWindow, yWindow = yWindow, spatialCoords = spatialCoords, groups = groups, na.ignore = na.ignore)
      } else {
        #create data.table of data passed to spatialSD and combine with data passed to spatialAverage
        ADCP_SD <- spatialSD(data = ADCP, xWindow = xWindow, yWindow = yWindow, spatialCoords = spatialCoords, groups = groups, na.ignore = na.ignore)
        setnames(ADCP_SD, names(ADCP_SD), paste(names(ADCP_SD), "_SD", sep = ""))
        ADCP <- data.table(ADCP, ADCP_SD)
      }
    }
  } else {
    
  }
  
  class(ADCP) <- c(class(ADCP), "adcp.secondary")
  
  return(ADCP)
}