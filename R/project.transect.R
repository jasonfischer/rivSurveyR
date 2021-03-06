#' project.transect
#' 
#' Orthogonally projects x and y coordinates of replicated ADCP transects to a mean transect line.
#' @param x	A vector of X geographic coordinates
#' @param y	A vector of Y geographic coordinates
#' @param transectNames	A vector identifying the transect corresponding to each x and y coordinate. Individual replicates should not be identified, see 'details' for additional information.
#' @details	\code{transectNames} is used to group replicate transects together during processing, if each replicate has a unique name all replicates will be handled as unique transects and projected coordinates will not represent the mean transect line. If \code{transectNames} is missing, all coordinates are assumed to correspond to a single transect.
#' @return	A \link[data.table]{data.table} of sample id, transect names, x coordinates fitted to an orthogonal plane (xProj), and y coordinates fitted to an orthogonal plane (yProj).
#' @export
#' @examples
#' data(vels)
#' plot(UTM_Y~UTM_X, vels)
#' mNineProject <- project.transect(vels$UTM_X, vels$UTM_Y, vels$transectName)
#' plot(yProj~xProj, mNineProject)
#' @seealso	\link{ortho.proj}
 
project.transect <- function(x, y, transectNames){
  if(missing(transectNames) == TRUE){
    transectNames <- "transect"
  } else {
    
  }
  projection <- list()
  dataset <- data.table(id = seq_along(x), x, y, transectNames)
  #create a column containing both x and y spatial coordinates, this will be used to remove overlapping cells, which is particular useful when processing secondary velocities
  dataset$xy <- paste(dataset$x, ",", dataset$y)
  uniqueNames <- unique(unlist(transectNames))
  #iteratively determine mean transect line for each replicated transect
  for(i in seq_along(uniqueNames)){
    transectData <- dataset[which(transectNames == uniqueNames[i]),]
    #condense dataset to remove redundancy from overlapping cells
    transectData <- transectData[,lapply(.SD, mean, na.rm = TRUE), by = list(xy, transectNames)]
    proj <- lm(transectData$y ~ transectData$x)
    #add orthogonally projected cell coordinates to a data.table with raw spatial coordinates
    if (is.na(proj$coefficients[2]) == TRUE | summary(proj)$sigma > diff(range(transectData$x, na.rm=TRUE))){
      #x and y coords are supplied backwards to ensure regression runs along the length of the transect
      proj <- lm(transectData$x ~ transectData$y)
      orthoProj <- ortho.proj(x = dataset[which(transectNames == uniqueNames[i]), y], y = dataset[which(transectNames == uniqueNames[i]), x], m = proj$coefficients[2], b = proj$coefficients[1])
      setnames(orthoProj, c("xProj", "yProj"), c("yProj", "xProj")) 
      projection[[i]] <- data.table(dataset[which(transectNames == uniqueNames[i]), c("id", "transectNames"), with = FALSE], xProj=orthoProj$xProj, yProj=orthoProj$yProj)
      warning(paste("Transect:", transectData$transectNames[1],
                    "Assumed to be Oriented North-South"))
    } else if(summary(proj)$coefficients[2,4] > 0.05){
      projection[[i]] <- data.table(dataset[which(transectNames == uniqueNames[i]), c("id", "transectNames"), with = FALSE], ortho.proj(x = dataset[which(transectNames == uniqueNames[i]), x], y = dataset[which(transectNames == uniqueNames[i]), y], m = proj$coefficients[2], b = proj$coefficients[1]))
      warning(paste("p > 0.05 for Linear Model of Transect:", transectData$transectNames[1],
                    "Projection May be Poor Fit or Transect is Oriented East-West"))
    } else {     
      projection[[i]] <- data.table(dataset[which(transectNames == uniqueNames[i]), c("id", "transectNames"), with = FALSE], ortho.proj(x = dataset[which(transectNames == uniqueNames[i]), x], y = dataset[which(transectNames == uniqueNames[i]), y], m = proj$coefficients[2], b = proj$coefficients[1]))
    }
  }  
    xyProj <- do.call("rbind", projection)
    xyProj <- xyProj[order(id), , ]
  return(xyProj)
}