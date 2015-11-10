#' average.secondary
#' 
#' Applies a function summarizing replicated ADCP data and calculates transect direction (from river right) and distance along transect.
#' @param data	Either a list of MATLAB files exported from RiverSurveyor to be compiled or an object of class "adcp.planform".
#' @param transectNames	A list of transect names corresponding to the transects represented by the MATLAB files. See 'details' for additional information.
#' @param FUN	Function used to summarize ADCP data, default is \link{mean}
#' @param binWidth	Width used to bin data along a transect, for summarization. Default is the mean distance between samples.
#' @param binHeight	Depth interval used to bin data within an ensemble, for summarization. If missing, cell depth is used group cells for summarization.
#' @param project	Logical.  If TRUE (default), transect replicates are projected to a mean transect line using an orthogonal projection of the x and y coordinates. If FALSE transects are not projected to a mean transect line.
#' @details	When assigning transectNames, order matters, the transect names must be in the same order of the MATLAB files assigned to the data argument. Ex., if the first two files in the list are Transect1a and Transect1b, which are replicates of Transect1, t first two files in the name list are "Transect1" and "Transect1".  The processing functions will spatially average and compile all data from the same transects, this tells the processing functions to average and compile these two files. Listing the names as "Transect1a" and "Transect1b" results in the transects being processed separately.  \code{transectNames} will only expect characters and names must be enclosed in "".
#' @return	A \link[data.table]{data.table} with the transect name (transectName), the distance from the right bank (tDist), cell depth, cell id, depth averaged velocity to the east (Mean.Vel.E), depth averaged velocity to the north (Mean.Vel.N), depth, east velocity within a cell (Vel.E), north velocity within a cell (Vel.N), vertical velocity within a cell (Vel.up), error velocity within a cell (Vel.error), blanking distance or depth where measurements first begin (StartDepth), cell height, cell number within an ensemble (cell one is closest to the surface; cellNumber), cell depth, UTM x coordinates (UTM_X), UTM y coordinates (UTM_Y), longitude, latitude, heading of transect (assumes transect starts at river right), cellWidth, speed within a cell (speed), flow heading within a cell (velHeading), and cell altitude. If project = TRUE, orthogonally projected longitude, latitude, and UTM coordinates (Longitude_Proj, Latitude_Proj, UTM_X_Proj, UTM_Y_Proj) are also returned.
#' @export
#' @examples
#' data(mNine)
#' #mNine is a list of MATLAB files
#' names(mNine)
#' #Drop the last letter (l or r) from the MATLAB file names
#' tNames <- substr(names(mNine), 0, nchar(names(mNine))-1)
#' average.secondary(mNine, tNames)
#' @seealso	\link{xSec.secondary} and \link{transectHeading}.
 
average.secondary <- function(data, transectNames, project = TRUE, binWidth, binHeight, FUN = mean){
  #if data passed to average.secondary is class adcp.secondary (indicating it has been processed by xSec.secondary), proceed to summarization, otherwise, data is passed to xSec.secondary
  if("adcp.secondary" %in% class(data)){
    x <- data
  } else if(is.list(data) == TRUE){
    x <- xSec.secondary(data, transectNames = transectNames, project = project)
  } else {
    stop("data must be list of objects from RiverSurveyor .mat output or of class adcp.secondary")
  }
  
  #create a list to populate with summarized values for each transect
  ADCP <- list()
  #iteratively summarize values for each transect
  for(i in seq_along(unique(x$transectName))){ 
    ADCP[[i]] <- x[which(transectName == unique(x$transectName)[i]),]
    #create a object for projected UTM coordinates or raw UTM coordinates (if transects were not projected) to be evaluated in following functions
    if("UTM_X_Proj" %in% names(ADCP[[i]]) & "UTM_Y_Proj" %in% names(ADCP[[i]])){    
      xUTM <- quote(ADCP[[i]]$UTM_X_Proj)
      yUTM <- quote(ADCP[[i]]$UTM_Y_Proj)
    } else {
      xUTM <- quote(ADCP[[i]]$UTM_X)
      yUTM <- quote(ADCP[[i]]$UTM_Y)
    }
    
    tHead <- transectHeading(x = eval(xUTM), y = eval(yUTM), velE = ADCP[[i]]$Vel.E, velN = ADCP[[i]]$Vel.N, depth = ADCP[[i]]$cellHeight)
    
    ADCP[[i]]$transectHeading <- tHead
    
    #determine the starting coordinates of the transect, assuming transect begins at river right
    if(tHead > 90 & tHead < 270){
      yStart <- max(eval(yUTM), na.rm=TRUE)
    } else {
      yStart <- min(eval(yUTM), na.rm=TRUE)
    }
    if(tHead > 180){
      xStart <- max(eval(xUTM), na.rm=TRUE)
    } else {
      xStart <- min(eval(xUTM), na.rm=TRUE)
    }
    
    if(missing(binWidth)==TRUE){
      sampleDistance <- ADCP[[i]][,mean(distance, na.rm = TRUE), by = sampleNum]
      binWidth <- round(mean(diff(sampleDistance$V1), na.rm = TRUE), digits = 1)
    } else {
    }
    ADCP[[i]] <- ADCP[[i]][,!c("distance"), with = FALSE]
    
    ADCP[[i]]$tDist <- sqrt((eval(xUTM) - xStart)^2 + (eval(yUTM) - yStart)^2)
    ADCP[[i]]$tDist <- round(ADCP[[i]]$tDist/binWidth) * binWidth
    
    if(missing(binHeight) == TRUE){
      binHeight <- round(mean(ADCP[[i]]$cellHeight, na.rm = TRUE), digits = 1)
    } else {      
    }
    ADCP[[i]]$cellDepth <- round(ADCP[[i]]$cellDepth/binHeight) * binHeight
    
    ADCP[[i]]$cellWidth <- mean(diff(sort(unique(ADCP[[i]]$tDist))))

  }
  #compile all summarized transect data into data.table  
  ADCP <- data.table(do.call("rbind", ADCP))

  #if the summarizing function allows additional arguments to be passed to it (...), include na.rm=TRUE in the argument line, otherwise, na.rm is excluded and NA values contribute to the summary. If a function has ... in the accepted arguments, but does use na.rm, this argument is ignored by the function  
  if(grepl("...", deparse(formals(FUN)), fixed=TRUE)==TRUE){
    ADCP <- ADCP[,lapply(.SD, match.fun(FUN), na.rm = TRUE), by = list(transectName, tDist, cellDepth)]
  } else {
    ADCP <- ADCP[,lapply(.SD, match.fun(FUN)), by = list(transectName, tDist, cellDepth)]
  }
  
  #Calculate water speed and heading (0 = North, 90 = East, 180 = South, 270 = West)
  ADCP$speed <- sqrt(ADCP$Vel.E^2 + ADCP$Vel.N^2 + ADCP$Vel.up^2)
  ADCP$velHeading <- heading(ADCP$Vel.E, ADCP$Vel.N)
  
  ADCP$cellAlt <- ADCP$Altitude - ADCP$cellDepth

  ADCP <- ADCP[order(ADCP$transectName, ADCP$tDist, ADCP$cellDepth),]
  
  class(ADCP) <- c(class(ADCP), "adcp.secondary")
  
  return(ADCP)
}