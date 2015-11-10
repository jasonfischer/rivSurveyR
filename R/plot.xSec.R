#' plot.xSec
#' 
#' Plots cross-section velocities given an object of class "adcp.secondary"
#' @param data	An object of class "adcp.secondary"
#' @param x	Parameter to be plotted
#' @param transects	A list of transects to be plotted, if missing all transects are plotted
#' @param ...	Arguments to be passed to \link{plot.xSec.default}.
#' @details	Use transects argument to subset data by transect, may produce error if data is subset prior to being supplied to function.
#' @return	A list of levelplot plots. Transect names contained within the dataset are used to name each element (i.e., plot) within the list.
#' @export
#' @examples
#' data(mNine)
#' #mNine is a list of MATLAB files
#' names(mNine)
#' #Drop the last letter (l or r) from the MATLAB file names
#' tNames <- substr(names(mNine), 0, nchar(names(mNine))-1)
#' adcp.secondary <- process.secondary(mNine, tNames, rotation = c("xSec", "rozovskii", "zeroSecQ"), xWindow = 21, yWindow = 5)
#' secondaryQ.zsq <- plot.xSec(adcp.secondary, vp.zsq)
#' #Display first cross-section
#' secondaryQ.zsq[[1]]
#' #Alternatively
#' secondaryQ.zsq$t1
#' @seealso	\link{plot.xSec.default}, \link[lattice]{levelplot}, \link{process.secondary}.

plot.xSec <- function(data, x, transects, ...){
  #if transects is omitted, plot all unique transects in the dataset, otherwise, only plot the transects specified in the transects argument
  if(missing(transects) == TRUE){
    transects <- sort(unique(data$transectName))
  } else {
    transects <- transects
  }
  #create a list to populate with plots of class trellis
  tplots <- list()
  #iteratively pass each transect to plot.xSec.default and add to tplots
  for(i in seq_along(transects)){
    tplot <- eval(substitute(plot.xSec.default(data = data, x=x, transect = transects[[i]], ...)))  
    tplots[[transects[[i]]]] <- tplot
  }
  names(tplots) <- transects
  return(tplots)
}