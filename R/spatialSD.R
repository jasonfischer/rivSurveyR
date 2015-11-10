#' spatialSD
#' 
#' Computes the standard deviation of spatially referenced data using a moving window
#' @param data	Dataset containing data used to calculate standard deviation and x, y coordinates.
#' @param xWindow	Number of neighbors to include in the x orientation of the moving window, including the observation the window centers on.
#' @param yWindow	Optional. Number of neighbors to include in the y orientation of the moving window, including the observation the window centers on.
#' @param spatialCoords	A vector of column name(s) containing spatial coordinates within the dataset. The x coordinates must be designated first, followed by the y coordinates. If \code{yWindow} is not provided, the y coordinates can be omitted, otherwise they will be ignored. If \code{spatialCoords} is omitted, a moving window will be used, but data will not be spatially referenced prior to calculations.
#' @param groups	Optional. Names of column(s) containing information used to subset the data into seperate groups for processing. Useful if dataset contains multiple discrete study areas.
#' @param na.ignore	Logical. If FALSE (default) NA values are assumed to be the standard deviation of non-NA values within the moving window, if TRUE NA values remain NA.
#' @return	A \link[data.table]{data.table} of spatially calculated standard deviation values
#' @export
#' @seealso	\link[caTools]{runsd}
#' @examples
#' data(vels)
#' #Standard deviation of velocities along ADCP transects
#' velSpatialSD <- spatialSD(vels, 21, spatialCoords = "tDist", groups = "transectName")
#' #Standard devation of velocities along ADCP transects and ensambles
#' data(cellVels)
#' velSpatialSD <- spatialSD(cellVels, 21, 5, spatialCoords = c("tDist","cellDepth"), groups = "transectName")

spatialSD <- function(data, xWindow, yWindow, spatialCoords, groups, na.ignore = FALSE) {
  message("Spatial Standard Deviation Being Calculated, Please Be Patient")
  flush.console()
  datT <- data.table(data) #NA values are assumed to be the standard deviation of the values within the moving window...minimizes data loss
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
  #create an object containing the names of objects to calculate spatial standard deviation
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
    #create a data.table with columns containing the spatialCoords values and data where spatial standard deviation is calculated in the x direction
    datT <- data.table(
      setnames(datT[, list(index, get(spatialCoords[1])), by = eval(bquote(.(spatialCoords[-1])))], c( spatialCoords[-1], "index", spatialCoords[1])),
      datT[,lapply(.SD, runsd, k= xWindow, endrule="sd", align = "center"), .SDcols = avgNames, by = eval(bquote(.(spatialCoords[-1])))]
    )
    datT <- datT[,!(1:(length(spatialCoords)-1)), with=FALSE] #remove the duplicate group variables


    #order by the second list item of spatialCoords    
    datT <- datT[order(datT[[spatialCoords[2]]]),]
    
    avgNames <- which(!names(datT)%in% spatialCoords & !names(datT)=="index")

    #create a data.table with columns containing the spatialCoords values and data where spatial standard deviation is calculated in the y direction    
    datT <- data.table(
      setnames(datT[, list(index, get(spatialCoords[2])), by = eval(bquote(.(spatialCoords[-2])))], c(spatialCoords[-2], "index", spatialCoords[2])),
      datT[,lapply(.SD, runsd, k=yWindow, endrule="sd", align = "center"), .SDcols = avgNames, by = eval(bquote(.(spatialCoords[-2])))]
    )
    datT <- datT[,!(1:(length(spatialCoords)-1)), with=FALSE] #remove the duplicate group variables
    
    message("Two Dimensional Standard Deviation Calculated")
    message(c("Window Size: ", xWindow, "x", yWindow, " cells"))
    flush.console()
  } else {
    #spatial averaging in one direction
    #create a data.table with columns containing the spatialCoords values and data where spatial standard deviation is calculated in the x direction
    if(missing(spatialCoords) == FALSE & missing(groups) == FALSE){
      #calculate spatial standard deviation for each group separately
      datT <- data.table(
        setnames(datT[, list(index, get(spatialCoords[1])), by = eval(bquote(.(spatialCoords[-1])))], c(spatialCoords[-1], "index", spatialCoords[1])),
        datT[,lapply(.SD, runsd, k=xWindow, endrule="sd", align = "center"), .SDcols = avgNames, by = eval(bquote(.(spatialCoords[-1])))]
      )
      datT <- datT[,!(1:(length(spatialCoords)-1)), with=FALSE] #remove the duplicate group variables
    } else if(missing(groups) == TRUE & missing(spatialCoords) == FALSE){
      #calculate spatial standard deviation for all data together
      #create a data.table with columns containing the spatialCoords values and data where spatial standard deviation is calculated in the x direction
      datT <- data.table(
        setnames(datT[, list(index, get(spatialCoords[1])), by = eval(bquote(.(spatialCoords[-1])))], c("index", spatialCoords[1])),
        datT[,lapply(.SD, runsd, k=xWindow, endrule="sd", align = "center"), .SDcols = avgNames]
      )
    } else {
    #create a data.table with columns containing the index values and spatial standard deviation data by the index value, assumes data supplied has been ordered correctly
      datT <- data.table(
        datT[, index,],
        datT[,lapply(.SD, runsd, k=xWindow, endrule="sd", align = "center"), .SDcols = avgNames]
      )
    }
    
    message("One Dimensional Standard Deviation Calculated")
    message(c("Window Size: ", xWindow, " cells"))
    flush.console()
  } 
  if(na.ignore == TRUE){
    datT <- datT[order(datT$index),]
    datT <- datT[,!c("index"), with=FALSE]
    for (i in 1:ncol(datT)){
      datT[which(is.na(data[[i]])==TRUE), eval(bquote(.(columnNames[i])))] <- NA
    }
  }else{
    datT <- datT[,!c("index"), with=FALSE]
  }
  return(datT)
}