#' process.planform
#' 
#' Summarizes replicated ADCP data, calculates water speeds, flow headings, transect direction (from river right), distance along transect, shear velocities, and conducts spatial averaging.
#' @param data	A list of MATLAB files exported from RiverSurveyor to be compiled
#' @param transectNames	A list of transect names corresponding to the transects represented by the MATLAB files. See 'details' for additional information.
#' @param FUN	Function used to summarize ADCP data, default is \link{mean}
#' @param binWidth	Width used to bin data along a transect, for summarization. Default is the mean distance between samples.
#' @param depthReference	Defines the depth measurements to be used. The default ("unit") uses the method defined in RiverSurveyor, but \code{depthReference} can be set to use measurements from the vertical beam ("VB"), bottom-track beams ("BT"), or both ("composite").
#' @param layerReference	Defines layerAvg bounds x1 and x2 as a distance from surface ("surface", default), distance above bottom ("bottom"), or percent of total depth ("percent").
#' @param x1	Lower bound velocities are averaged over. Either a numeric representing a depth or percentage or a character where "minimum" equals the minimum distance from the reference point (i.e., 0 when \code{layerReference} = "surface" and bottomCellDepth when \code{layerReference} = "bottom") and "maximum" equals maximum distance from the reference point. Must be numeric if \code{layerReference} = "percent". Defaults to "minimum".
#' @param x2	Upper bound velocities are averaged over. Either a numeric representing a depth or percentage or a character where "minimum" equals the minimum distance from the reference point (i.e., 0 when \code{layerReference} = "surface" and bottomCellDepth when \code{layerReference} = "bottom") and "maximum" equals maximum distance from the reference point. Must be numeric if \code{layerReference} = "percent". Defaults to "maximum".
#' @param project	Logical.  If TRUE (default), transect replicates are projected to a mean transect line using an orthogonal projection of the x and y coordinates. If FALSE transects are not projected to a mean transect line.
#' @param Dc	Sediment diameter of 90th percentile. If NULL (default) shear velocities are not calculated.
#' @param n	If NULL (default) shear velocities are not calculated.
#' @param p	Water density (kg/m^3)
#' @param reference	Method used to estimate shear velocity. If "mean" the Keulegan method is used where shear velocity = U/(1/k * ln(11* (depth/kc))), where U is the depth averaged speed, k is von Karman's constant (0.41), and, kc=n/Dc. For details see Garcia (2008). If "bottomCell" (the default) shear velocity = u/(9.5 * (h/kc)^(1/6)), where u is speed at h m above bottom. For details see Cheng-Lung (1991) and Simpson and Oltmann (1990).
#' @param units	Unit system of velocities, depths, and sediment diameters. Can either be "metric" (default) or "standard".
#' @param spatialAvg	Logical. If TRUE (default) and xWindow is not NULL, spatial averaging is conducted, else if FALSE no spatial averaging is conducted.
#' @param spatialStdDev	Logical. If TRUE and xWindow is not NULL, spatial standard deviation is calculated, else if FALSE (default) no spatial standard deviation is calculated.
#' @param xWindow	Number of neighbors to include in the x orientation of the moving window, including the observation the window centers on. If NULL (the default), spatial averaging is not conducted.
#' @param groups	Optional. Names of column(s) containing information used to subset the data into separate groups for processing. Useful if dataset contains multiple discrete study areas.
#' @param na.ignore	Logical. If FALSE (default) NA values are assumed to be the average of non-NA values within the moving window, if TRUE NA values remain NA.
#' @details	When assigning transectNames, order matters, the transect names must be in the same order of the MATLAB files assigned to the data argument. Ex., if the first two files in the list are Transect1a and Transect1b, which are replicates of Transect1, t first two files in the name list are "Transect1" and "Transect1".  The processing functions will spatially average and compile all data from the same transects, this tells the processing functions to average and compile these two files. Listing the names as "Transect1a" and "Transect1b" results in the transects being processed separately.  \code{transectNames} will only expect characters and names must be enclosed in "".
#' \code{FUN} must be \link{mean} to calculate shear velocities and conduct spatial averaging or spatial standard deviation
#' @return	A \link[data.table]{data.table} with the distance from the right bank (tDist), the transect name (transectName), UTM x coordinates (UTM_X), UTM y coordinates (UTM_Y), depth, depth averaged velocity to the east (Mean.Vel.E), depth averaged velocity to the north (Mean.Vel.N), depth of the bottom most cell in an ensemble (bottomCellDepth), distance of the bottom cell to the bottom (distToBottom), bottom cell velocity to the east (BC.Vel.E), bottom cell velocity to the north (BC.Vel.N), bottom cell vertical velocity (BC.Vel.Up), bottom cell error velocity (BC.Vel.Error), layer averaged velocity to the east (layerVel.E), layer averaged velocity to the north (layerVel.N), layer averaged vertical velocity (layerVel.Up), longitude, latitude, altitude of water surface, temperature, depth averaged speed, depth averaged flow heading (Heading), layer averaged speed, layer averaged flow heading (LayerAvgSpeedHeading), discharge of ensembles (SpecificQ), shear velocity to the east (ustarE), shear velocity to the north (ustarN), shear speed (ustar), and shear speed heading (ustarHeading). If project = TRUE, orthogonally projected longitude, latitude, and UTM coordinates (Longitude_Proj, Latitude_Proj, UTM_X_Proj, UTM_Y_Proj) are also returned.
#' @export
#' @examples
#' data(mNine)
#' #mNine is a list of MATLAB files
#' names(mNine)
#' #Drop the last letter (l or r) from the MATLAB file names
#' tNames <- substr(names(mNine), 0, nchar(names(mNine))-1)
#' adcp.planform <- process.planform(mNine, tNames, Dc = 0.0375, n = 2, xWindow = 21)
#' @seealso	\link{average.planform}, \link{shearVel}, \link{spatialAverage}, and \link{spatialSD}.
 
process.planform <- function(data, transectNames, FUN = mean, binWidth, depthReference = "unit", layerReference = "surface", x1 = "minimum", x2 = "maximum", project = TRUE,
                             Dc = NULL, n  = NULL, p = 1000, reference = "bottomCell", units = "metric",
                             spatialAvg = TRUE, spatialStdDev = FALSE,
                             xWindow = NULL, spatialCoords = c("tDist"), groups = c("transectName"), na.ignore = FALSE){  
  message("Compiling and Averaging Data")
  flush.console()
  #Pass data to average.planform
  ADCP <- average.planform(data = data, transectNames = transectNames, FUN = FUN, binWidth = binWidth, depthReference = depthReference, layerReference = layerReference, x1 = x1, x2 = x2, project = project)

  ADCP <- ADCP[order(ADCP$transectName, ADCP$tDist),]

  if(as.character(substitute(FUN)) == "mean") {
    if(is.null(Dc) == TRUE | is.null(n) == TRUE) {
      ADCP <- ADCP
    } else {
      #if Dc and n are provided, pass data to shearVel and calculate shear velocity
      message("Calculating Bottom Velocities")
      flush.console()
      if(reference == "bottomCell") {
        botVel <- shearVel(ve = ADCP$BC.Vel.E, vn = ADCP$BC.Vel.N, vu = ADCP$BC.Vel.Up, depth = ADCP$depth, bcDepth = ADCP$bottomCellDepth, Dc = Dc, n = n, p = p, reference = reference)
      } else {
        botVel <- shearVel(ve = ADCP$Mean.Vel.E, vn = ADCP$Mean.Vel.N, depth = ADCP$depth, Dc = Dc, n = n, p = p, reference = reference)
      }
      ADCP <- data.table(ADCP, botVel)
    }
    if(is.null(xWindow) == TRUE | spatialAvg == FALSE & spatialStdDev == FALSE){
      ADCP <- ADCP
    } else {
      if(spatialStdDev == FALSE & spatialAvg == TRUE){
      #pass data to spatialAverage
      ADCP <- spatialAverage(ADCP, xWindow = xWindow, spatialCoords = spatialCoords, groups = groups, na.ignore = na.ignore)
      } else if(spatialStdDev == TRUE & spatialAvg == FALSE){
        #pass data to spatialSD
        ADCP <- spatialSD(ADCP, xWindow = xWindow, spatialCoords = spatialCoords, groups = groups, na.ignore = na.ignore)
      }else{
        #create data.table of data passed to spatialSD and combine with data passed to spatialAverage
        ADCP_SD <- spatialSD(ADCP, xWindow = xWindow, spatialCoords = spatialCoords, groups = groups, na.ignore = na.ignore)
        setnames(ADCP_SD, names(ADCP_SD), paste(names(ADCP_SD), "_SD", sep = ""))
        ADCP <- spatialAverage(ADCP, xWindow = xWindow, spatialCoords = spatialCoords, groups = groups, na.ignore = na.ignore)
        ADCP <- data.table(ADCP, ADCP_SD)
      }
    }
  } else {
    ADCP <- ADCP
  }
  
  class(ADCP) <- c(class(ADCP), "adcp.planform")
  
  return(ADCP)
} 