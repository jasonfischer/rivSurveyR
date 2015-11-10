#' spatialAverage
#' 
#' Averages spatially referenced data using a moving window
#' @param data	Dataset containing data to be averaged and x, y coordinates
#' @param xWindow	Number of neighbors to include in the x orientation of the moving window, including the observation the window centers on.
#' @param yWindow	Optional. Number of neighbors to include in the y orientation of the moving window, including the observation the window centers on.
#' @param spatialCoords	A vector of column name(s) containing spatial coordinates within the dataset. The x coordinates must be designated first, followed by the y coordinates. If \code{yWindow} is not provided, the y coordinates can be omitted, otherwise they will be ignored. If \code{spatialCoords} is omitted, a moving window will be used, but data will not be spatially referenced prior to calculations.
#' @param groups	Optional. Names of column(s) containing information used to subset the data into separate groups for processing. Useful if dataset contains multiple discrete study areas.
#' @param na.ignore	Logical. If FALSE (default) NA values are assumed to be the average of non-NA values within the moving window, if TRUE NA values remain NA.
#' @return	A \link[data.table]{data.table} of spatially averaged values
#' @export
#' @seealso	\link[caTools]{runmean}
#' @examples
#' data(vels)
#' #Average velocities along ADCP transects
#' velSpatialMean <- spatialAverage(vels, 21, spatialCoords = "tDist", groups = "transectName")
#' #Average velocities along ADCP transects and ensembles
#' data(cellVels)
#' velSpatialMean <- spatialAverage(cellVels, 21, 5, spatialCoords = c("tDist","cellDepth"), groups = "transectName")
 
spatialAverage <- function(data, xWindow, yWindow, spatialCoords, groups, na.ignore = FALSE) {
  message("Spatial Averaging in Process, Please Be Patient")
  flush.console()
  datT <- data.table(data)
  #create a column containing a unique value for each row, this will be used to ensure NA values remain NA if na.ignore=TRUE
  datT$index <- seq_along(datT[[1]])
  columnNames <- names(datT)
  
  #order the dataset by the first list item supplied to spatialCoords, ensures neighbors in space are neighbors in the dataset
  if(missing(spatialCoords) == FALSE){
    datT <- datT[order(datT[[spatialCoords[1]]]),]    
    if(missing(yWindow) == TRUE){
      spatialCoords <- spatialCoords[1]
    }else{      
    }
  } else {
    datT <- datT
  }

  #create an object containing the names of objects to be spatially averaged 
  if(missing(spatialCoords) == FALSE){
    if(missing(groups) == FALSE){
      spatialCoords <- c(spatialCoords, groups)
      avgNames <- which(!names(datT)%in% spatialCoords & !names(datT)=="index")
      } else {        
        avgNames <- which(!names(datT)%in% spatialCoords & !names(datT)=="index")
      }
  } else {
    avgNames <- names(datT)[which(!names(datT)=="index")]
  }

  #Determine if spatial averaging should be conducted in one or two directions
  if(missing(xWindow)==FALSE & missing(yWindow)==FALSE & length(spatialCoords)>2){
    #spatial averaging in two directions
    #create a data.table with columns containing the spatialCoords values and data spatially averaged in the x direction
    datT <- data.table(
      setnames(datT[, list(index, get(spatialCoords[1])), by = eval(bquote(.(spatialCoords[-1])))], c( spatialCoords[-1], "index", spatialCoords[1])),
      datT[,lapply(.SD, runmean, k= xWindow, alg="C", endrule="mean", align = "center"), .SDcols = avgNames, by = eval(bquote(.(spatialCoords[-1])))]
    )

    datT <- datT[,!(1:(length(spatialCoords)-1)), with=FALSE] #remove the duplicate group variables

    #order by the second list item of spatialCoords
    datT <- datT[order(datT[[spatialCoords[2]]]),]

    avgNames <- which(!names(datT)%in% spatialCoords & !names(datT)=="index")
    #create a data.table with columns containing the spatialCoords values and data spatially averaged in the y direction
    datT <- data.table(
      setnames(datT[, list(index, get(spatialCoords[2])), by = eval(bquote(.(spatialCoords[-2])))], c(spatialCoords[-2], "index", spatialCoords[2])),
      datT[,lapply(.SD, runmean, k=yWindow, alg="C", endrule="mean", align = "center"), .SDcols = avgNames, by = eval(bquote(.(spatialCoords[-2])))]
    )

    datT <- datT[,!(1:(length(spatialCoords)-1)), with=FALSE] #remove the duplicate group variables
    
    message("Two Dimensional Spatial Averaging Complete")
    message(c("Window Size: ", xWindow, "x", yWindow, " cells"))
    flush.console()
  }else{
    #spatial averaging in one direction
    if(missing(spatialCoords) == FALSE & missing(groups) == FALSE){
      #spatially average each group separately
      datT <- data.table(
        setnames(datT[, list(index, get(spatialCoords[1])), by = eval(bquote(.(spatialCoords[-1])))], c(spatialCoords[-1], "index", spatialCoords[1])),
        datT[,lapply(.SD, runmean, k=xWindow, alg="C", endrule="mean", align = "center"), .SDcols = avgNames, by = eval(bquote(.(spatialCoords[-1])))]
      )

      datT <- datT[,!(1:(length(spatialCoords)-1)), with=FALSE] #remove the duplicate group variables
    } else if(missing(groups) == TRUE & missing(spatialCoords) == FALSE){
      #spatially average all data together
      datT <- data.table(
        setnames(datT[, list(index, get(spatialCoords[1])), by = eval(bquote(.(spatialCoords[-1])))], c("index", spatialCoords[1])),
        datT[,lapply(.SD, runmean, k=xWindow, alg="C", endrule="mean", align = "center"), .SDcols = avgNames]
      )
    } else {
    #create a data.table with columns containing the spatialCoords values and data spatially averaged by the index value, assumes data supplied has been ordered correctly
      datT <- data.table(
        datT[, index,],
        datT[,lapply(.SD, runmean, k=xWindow, alg="C", endrule="mean", align = "center"), .SDcols = avgNames]
      )
    }

    message("One Dimensional Spatial Averaging Complete")
    message(c("Window Size: ", xWindow, " cells"))
    flush.console()
  }
  #If na.ignore is TRUE NA values are backfilled into the spatially averaged data, otherwise, NA values are assumed to be the average of the values within the moving window
  if(na.ignore == TRUE){ 
    datT <- datT[order(datT$index),]
    datT <- datT[,!c("index"), with=FALSE]
    for (i in 1:ncol(datT)){
      datT[which(is.na(data[[i]])==TRUE), eval(bquote(.(columnNames[i])))] <- NA
    }
  }else{
    datT <- datT[order(datT$index),]
    datT <- datT[,!c("index"), with=FALSE]
  }
  return(datT)
}