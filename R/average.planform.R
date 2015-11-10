#' average.planform
#' 
#' Applys a function summarizing replicated ADCP data and calculates water speeds, flow headings, transect direction (from river right), and distance along transect.
#' @param data	Either a list of MATLAB files exported from RiverSurveyor to be compiled or an object of class "adcp.planform".
#' @param transectNames	A list of transect names corresponding to the transects represented by the MATLAB files. See 'details' for additional information.
#' @param FUN	Function used to summarize ADCP data, default is \link{mean}
#' @param binWidth	Width used to bin data along a transect, for summarization. Default is the mean distance between samples.
#' @param layerReference	Defines layerAvg bounds x1 and x2 as a distance from surface ("surface", default), distance above bottom ("bottom"), or percent of total depth ("percent").
#' @param x1	Lower bound velocities are averaged over. Either a numeric representing a depth or percentage or a character where "minimum" equals the minimum distance from the reference point (i.e., 0 when \code{layerReference} = "surface" and bottomCellDepth when \code{layerReference} = "bottom") and "maximum" equals maximum distance from the reference point. Must be numeric if \code{layerReference} = "percent". Defaults to "minimum".
#' @param x2	Upper bound velocities are averaged over. Either a numeric representing a depth or percentage or a character where "minimum" equals the minimum distance from the reference point (i.e., 0 when \code{layerReference} = "surface" and bottomCellDepth when \code{layerReference} = "bottom") and "maximum" equals maximum distance from the reference point. Must be numeric if \code{layerReference} = "percent". Defaults to "maximum".
#' @param project	Logical.  If TRUE (default), transect replicates are projected to a mean transect line using an orthogonal projection of the x and y coordinates. If FALSE transects are not projected to a mean transect line.
#' @details	When assigning transectNames, order matters, the transect names must be in the same order of the MATLAB files assigned to the data argument. Ex., if the first two files in the list are Transect1a and Transect1b, which are replicates of Transect1, t first two files in the name list are "Transect1" and "Transect1".  The processing functions will spatially average and compile all data from the same transects, this tells the processing functions to average and compile these two files. Listing the names as "Transect1a" and "Transect1b" results in the transects being processed separately.  \code{transectNames} will only expect characters and names must be enclosed in "".
#' @return	A \link[data.table]{data.table} with the distance from the right bank (tDist), the transect name (transectName), UTM x coordinates (UTM_X), UTM y coordinates (UTM_Y), depth, depth averaged velocity to the east (Mean.Vel.E), depth averaged velocity to the north (Mean.Vel.N), depth of the bottom most cell in an ensemble (bottomCellDepth), distance of the bottom cell to the bottom (distToBottom), bottom cell velocity to the east (BC.Vel.E), bottom cell velocity to the north (BC.Vel.N), bottom cell vertical velocity (BC.Vel.Up), bottom cell error velocity (BC.Vel.Error), layer averaged velocity to the east (layerVel.E), layer averaged velocity to the north (layerVel.N), layer averaged vertical velocity (layerVel.Up), longitude, latitude, altitude of water surface, temperature, depth averaged speed, depth averaged flow heading (Heading), layer averaged speed, layer averaged flow heading (LayerAvgSpeedHeading), and discharge of ensembles (SpecificQ). If project = TRUE, orthogonally projected longitude, latitude, and UTM coordinates (Longitude_Proj, Latitude_Proj, UTM_X_Proj, UTM_Y_Proj) are also returned.
#' @export
#' @examples
#' data(mNine)
#' #mNine is a list of MATLAB files
#' names(mNine)
#' #Drop the last letter (l or r) from the MATLAB file names
#' tNames <- substr(names(mNine), 0, nchar(names(mNine))-1)
#' average.planform(mNine, tNames)
#' @seealso	\link{xSec.planform} and \link{transectHeading}.
 
average.planform <- function(data, transectNames, FUN = mean, binWidth, layerReference = "surface", x1 = "minimum", x2 = "maximum", project = TRUE){
  #if data passed to average.planform is class adcp.planform (indicating it has been processed by xSec.planform), proceed to summarization, otherwise, data is passed to xSec.planform
  if("adcp.planform" %in% class(data)){
    dt <- data
  } else if(is.list(data) == TRUE){
    dt <- xSec.planform(data, transectNames = transectNames, layerReference = layerReference, x1 = x1, x2 = x2, project = project)
  } else {
    stop("data must be list of objects from RiverSurveyor .mat output or of class adcp.planform")
  }
  
  #create a list to populate with summarized values for each transect
  ADCP <- list()
  #iteratively summarize values for each transect
  for(i in 1:length(unique(dt$transectName))){ 
    ADCP[[i]] <- dt[which(transectName == unique(dt$transectName)[i]),]
    #create an object for projected UTM coordinates or raw UTM coordinates (if transects were not projected) to be evaluated in following functions
    if("UTM_X_Proj" %in% names(ADCP[[i]]) & "UTM_Y_Proj" %in% names(ADCP[[i]])){
      xUTM <- quote(ADCP[[i]]$UTM_X_Proj)
      yUTM <- quote(ADCP[[i]]$UTM_Y_Proj)
    } else {
      xUTM <- quote(ADCP[[i]]$UTM_X)
      yUTM <- quote(ADCP[[i]]$UTM_Y)
    }
    tHead <- transectHeading(eval(xUTM), eval(yUTM), ADCP[[i]]$Mean.Vel.E, ADCP[[i]]$Mean.Vel.N, ADCP[[i]]$depth)
    
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
      binWidth <- round(mean(diff(ADCP[[i]]$distance), na.rm = TRUE), digits = 1)
    } else {
    }
    ADCP[[i]] <- ADCP[[i]][,!c("distance"), with = FALSE]
    
    #calculate the distance each sample is from the transect origin
    ADCP[[i]]$tDist <- sqrt((eval(xUTM) - xStart)^2 + (eval(yUTM) - yStart)^2)
    ADCP[[i]]$tDist <- round(ADCP[[i]]$tDist/binWidth) * binWidth
    
  }
  #compile all summarized transect data into data.table
  ADCP <- data.table(do.call("rbind", ADCP))
  keycols = c("transectName", "tDist")
  setkeyv(ADCP, keycols)
  
  #if the summarizing function allows additional arguments to be passed to it (...), include na.rm=TRUE in the argument line, otherwise, na.rm is excluded and NA values contribute to the summary. If a function has ... in the accepted arguments, but does use na.rm, this argument is ignored by the function
  if(grepl("...", deparse(formals(FUN)), fixed=TRUE)){
    ADCP <- ADCP[,lapply(.SD, match.fun(FUN), na.rm = TRUE), by = list(transectName, tDist)]
  } else {
    ADCP <- ADCP[,lapply(.SD, match.fun(FUN)), by = list(transectName, tDist)]
  }
 
  #Calculate depth averaged water speed and heading (0 = North, 90 = East, 180 = South, 270 = West)
  ADCP$DepthAveragedSpeed <- sqrt(ADCP$Mean.Vel.E^2 + ADCP$Mean.Vel.N^2)
  ADCP$Heading <- heading(ADCP$Mean.Vel.E, ADCP$Mean.Vel.N)

  ADCP$LayerAveragedSpeed <- sqrt(ADCP$layerVel.E^2 + ADCP$layerVel.N^2)
  ADCP$LayerAvgSpeedHeading <- heading(ADCP$layerVel.E, ADCP$layerVel.N)
   
  ADCP$SpecificQ <- ADCP$DepthAveragedSpeed * ADCP$depth

  ADCP <- ADCP[order(ADCP$transectName, ADCP$tDist),]

  class(ADCP) <- c(class(ADCP), "adcp.planform") 

  return(ADCP)
}