#' plot.planform
#' 
#' Plots spatial data from an object of class "adcp.planform"
#' @param data	An object of class "adcp.planform"
#' @param x	Parameter to be plotted as points using \code{xyplot}
#' @param interpolate	Optional. Parameter to plot using \code{levelplot}, if provided in addition to \code{x} the interpolated parameter is overlaid by \code{x}.
#' @param coord	Spatial coordinates used to define x and y-axes, can be orthogonally "projected.utm" (default), "utm", orthogonally "projected.latlong", or "latlong".
#' @param point.spacing	Distance between points, default is 1.
#' @param point.color	Color of points displayed, maybe numeric, character string, or function. Default is \code{jet.colors}.
#' @param arrows	If TRUE (default) and x is a velocity vector, arrows indicating flow direction are plotted
#' @param arrow_scale	Multiplication factor for sizing secondary flow arrows, default is 1
#' @param arrowwd	Line width of arrows
#' @param arrow.key	Numeric vector of length 2 defining the X and y coordinates (c(x,y)) of where scaling key for points or arrows (if \code{arrows} = TRUE) should be written. If omitted, a key is placed in the lower left.
#' @param keyCol	Color of points or arrows (if \code{arrows} = TRUE) in arrow.key, default is "black".
#' @param textCol	Colors of text in arrow.key, default is "black".
#' @param xUnits	Character string. If provided, the units of \code{x} (e.g., "m" or "m/s") to be displayed in arrow.key.
#' @param cellWidth	Width of cell size used to grid interpolate. If missing, defaults to 1/100 of the range of \code{xlim}.
#' @param cellHeight	Height of cell size used to grid interpolate. If missing, defaults to 1/100 of the range of \code{ylim}.
#' @param interpolate.color	Color function used to display \code{interpolate}, see \link[lattice]{levelplot} for details. Default is \code{jet.colots}.
#' @param xlab	X-axis label, see \link[lattice]{levelplot} for details. Default is "X Coordinates" in size 18 font.
#' @param ylab	Y-axis label, see \code{levelplot} for details. Default is "Y Coordinates" in size 18 font.
#' @param ylab.right	\code{colorkey} label, see \link[lattice]{xyplot} for details. Default is as.character(interpolate) in size 18 font.
#' @param colorkey	A list of arguments for the color key drawn alongside the plot, see \link[lattice]{levelplot} for details. Default writes the tick labels with cex 1.5.
#' @param scales	A list of arguments determining how the axes are drawn, see \link[lattice]{xyplot} for details. Default is to draw tick labels with cex 1.5.
#' @param key.cex	Font size of text in arrow.key, default is 1.25
#' @param xlim	Numeric vector of length 2 defining the range of the x-axis.
#' @param ylim	Numeric vector of length 2 defining the range of the y-axis.
#' @param aspect	Controls the aspect ratio of the panels, see \link[lattice]{xyplot} for details. Default is "iso".
#' @param par.settings	A list of arguments for fine-tuned control of display, see \link[lattice]{xyplot} and \link[lattice]{trellis.par.set} for details. Default is not to pad the color key, offset the color key label by 2.
#' @param ...	Additional arguments to be passed to \link[lattice]{xyplot} and \link[lattice]{levelplot}.
#' @param contour	Logical, if TRUE (default) and \code{interpolate} is defined, contour lines are drawn
#' @param cuts	Number of regions contour lines denote
#' @param labels	Logical, if TRUE (default) contour lines are labeled
#' @param contourwd	Line width of contours
#' @return A spatially referenced plot of \code{x} and \code{interpolate} (if provided).
#' @export
#' @examples
#' data(mNine)
#' #mNine is a list of MATLAB files
#' names(mNine)
#' #Drop the last letter (l or r) from the MATLAB file names
#' tNames <- substr(names(mNine), 0, nchar(names(mNine))-1)
#' adcp.planform <- process.planform(mNine, tNames, Dc = 0.0375, n = 2, xWindow = 21)
#' #Plot depth averaged speed with arrows displaying flow heading
#' plot.planform(adcp.planform, DepthAveragedSpeed, arrow.key = c(325540,4679130), xUnits = "m/s", arrow_scale = 10)
#' #Plot depth averaged speed as points
#' plot.planform(adcp.planform, DepthAveragedSpeed, arrow.key = c(325540,4679130), xUnits = "m/s", arrows = FALSE)
#' #Plot depth averaged speed overlaying depth
#' plot.planform(adcp.planform, DepthAveragedSpeed, interpolate = depth, arrow.key = c(325540,4679130), xUnits = "m/s", point.spacing = 10, point.color = "black", ylab.right = list(label = "Depth (m)", fontsize = 18), arrow_scale=10)
#' @seealso \link{process.planform},\link[lattice]{levelplot}, and \link[lattice]{xyplot}.

plot.planform <- function(data, x, interpolate, coord = "projected.utm",                       
                          point.spacing = 1, point.color = jet.colors,
                          arrows = TRUE, arrow_scale = 1, arrowwd = 1, 
                          arrow.key, keyCol = "black", textCol = "black", xUnits = NULL,
                          cellWidth, cellHeight, interpolate.color = jet.colors,                                                   
                          xlab = list(label = "X Coordinates", fontsize = 18), 
                          ylab = list(label = "Y Coordinates", fontsize = 18),
                          ylab.right = list(label = as.character(interpolate), fontsize = 18),
                          colorkey=list(labels=list(cex=1.25)), scales = list(cex = c(1.25, 1.25)),
                          key.cex = 1.25, 
                          xlim, ylim,                          
                          aspect = "iso", 
                          par.settings = list(layout.widths = list(axis.key.padding = 0, ylab.right = 2)),
                          ..., contour = TRUE, cuts = 7, labels = TRUE, contourwd = 1){
  stopifnot("adcp.planform" %in% class(data))
  
  #create an object for the x and y coordinates to be evaluated in following functions  
  if(coord == "projected.utm"){
    xCoordinate <- as.character(substitute(UTM_X_Proj))
    yCoordinate <- as.character(substitute(UTM_Y_Proj))
  } else if(coord == "utm"){
    xCoordinate <- as.character(substitute(UTM_X))
    yCoordinate <- as.character(substitute(UTM_Y))
  } else if(coord == "projected.latlong"){
    xCoordinate <- as.character(substitute(Longitude_Proj))
    yCoordinate <- as.character(substitute(Latitude_Proj))
  } else if(coord == "latlong"){
    xCoordinate <- as.character(substitute(Longitude))
    yCoordinate <- as.character(substitute(Latitude))
  } else {
    stop("coord must be projected.utm, utm, projected.latlong, or latlong")
  }  
  xCoord <- eval(quote(xCoordinate))
  yCoord <- eval(quote(yCoordinate))
  
  #Set scaling factor for legend arrows
  if(!"ylim" %in% (names(match.call(expand.dots=T)))){
    ylim <- range(data[[yCoord]])
  } else {
    ylim <- ylim
  }
  if(!"xlim" %in% (names(match.call(expand.dots=T)))){
    xlim <- range(data[[xCoord]])
  } else {
    xlim <- xlim
  }
  
  if(missing(x)==FALSE){
    #create an object for x to be evaluated in following functions
    x <- as.character(substitute(x))
    x <- eval(quote(x))
    data <- data.table(data)
    
    #determine the values to be displayed in the arrow key
    keyValues <- summary(data[,abs(get(x))])
    keyValues <- keyValues[which(!names(keyValues) == "Mean")]
    keyValues <- c(keyValues["Min."][[1]], keyValues["Median"][[1]], keyValues["Max."][[1]])
    keyValues <- round(as.numeric(keyValues), digits = 2)
    keyValues <- as.factor(keyValues)
    keyValues <- c(min(as.numeric(as.character(keyValues))), median(as.numeric(as.character(keyValues))), max(as.numeric(as.character(keyValues))))
  
    #create a data.table of the point.data to be displayed, with the defined spacing between points. Because points near (within 1 unit (e.g. meter) of the designated spacing distance are selected, multiple points may be included (e.g., points 9.6m, 9.8m, 10m, 10.2m, and 10.4m will all be considered near the 10m spacing interval), each spacing interval may have multiple data points. Thus the mean value of all data at each spacing interval is taken, so that only one point is plotted at each interval.
    point.data <- data[round(tDist)%%point.spacing==0, lapply(.SD, mean, na.rm=TRUE), by=list(transectName, round(tDist))]
    point.data <- point.data[!is.na(get(x)),,]

    #determine the color of points
    if(is.function(point.color)==TRUE){
      point.data <- point.data[base:::order(get(x)),,] #currently using the order function in base R, because data.table order requires the column name without quotations
      point.col <- point.color(nrow(unique(point.data[,x, with=FALSE]))) 
      point.col <- data.table(point.col, unique(point.data[,x, with=FALSE]))
      setkeyv(point.col, as.character(substitute(x)))
      setkeyv(point.data, as.character(substitute(x)))
      point.data<-merge(point.data, point.col)
      
      keyCol <- point.color(length(keyValues))
    }else{
      point.col <- point.color
      keyCol <- point.color
    }

    #determine location for the key  
    if("arrow.key" %in% (names(match.call(expand.dots=T))) == FALSE){
      arrow.key <- c(min(xlim), min(ylim))
      key.x <- min(xlim)+(diff(xlim)*0.03)
      key.y <- c((min(ylim)+(diff(ylim)*0.10*key.cex)), (min(ylim)+(diff(ylim)*0.06*key.cex)), min(ylim)+(diff(ylim)*0.02*key.cex))
    } else if(is.numeric(arrow.key) == TRUE){
      key.x <- arrow.key[1]+diff(xlim)*0.05
      key.y <- c((arrow.key[2] + (diff(ylim)*0.10*key.cex)), (arrow.key[2]+(diff(ylim)*0.06*key.cex)), arrow.key[2]+(diff(ylim)*0.02*key.cex))
    } else {
    }
    
    #if arrows is true and x is a velocity parameter, create a data vector of the distance moved to the east and north in one second, otherwise, x and y displacement is assumed to be zero    
    if(arrows == TRUE){
      if(x == "DepthAveragedSpeed"){
        xDis <- point.data$Mean.Vel.E
        yDis <- point.data$Mean.Vel.N
      } else if(x == "ustar"){
        xDis <- point.data$ustarE
        yDis <- point.data$ustarN
      } else if(x == "LayerAveragedSpeed"){
        xDis <- point.data$layerVel.E
        yDis <- point.data$layerVel.N
      } else {
        xDis <- 0
        yDis <- 0
        warning("x and y velocity vectors not found for ", x, ", unable to plot flow arrows")
      }
      
      #plot spatially referenced x using arrows
      planform.plot <- xyplot(get(yCoord) ~ get(xCoord), data = point.data, col = point.data$point.col, pch = 19,
                              ylim = ylim, xlim = xlim, aspect = aspect, ylab = ylab, xlab = xlab,
                              scales = scales, 
                              panel = function(...){
                                panel.arrows(x0 = point.data[[xCoord]], y0 = point.data[[yCoord]], x1 = (point.data[[xCoord]] + (xDis * arrow_scale)), y1 = (point.data[[yCoord]] + arrow_scale * yDis), length = 0.05,
                                             col = point.data$point.col, lwd = arrowwd,
                                             at = seq(min(point.data[[x]]), max(point.data[[x]]), length = length(unique(point.data[[x]])))
                                )
                                #plot scale arrows
                                if(is.numeric(arrow.key) == TRUE){
                                  panel.arrows(x0 = key.x, y0 = key.y, x1 = key.x + keyValues * arrow_scale, y1 = key.y, length = 0.05,
                                               col = keyCol, lwd = arrowwd)
                                  panel.text(x = key.x + diff(xlim)*0.02 ,y = key.y + diff(ylim)*0.02*key.cex, label = paste(keyValues, xUnits), cex = key.cex, col = textCol, adj = c(0,0.5))                              
                                }else{
                                  
                                }
                              },
                              ...)
    }else{
      #plot spatially referenced x using points and add a color scale
      if(is.function(point.color)==TRUE){
        planform.plot <- xyplot(get(yCoord) ~ get(xCoord), data = point.data, col = point.data$point.col, pch = 19, cex = arrowwd,
                                ylim = ylim, xlim = xlim, aspect = aspect, ylab = ylab, xlab = xlab,
                                at = seq(min(point.data[[x]]), max(point.data[[x]]), length = length(unique(point.data[[x]]))),
                                scales = scales,
                                panel = function(...){
                                  panel.xyplot(...)
                                  panel.points(x = key.x, y = key.y, col = keyCol, pch = 19, cex = arrowwd)
                                  panel.text(x = key.x + diff(xlim)*0.02 ,y = key.y + diff(ylim)*0.02*key.cex, label = paste(keyValues, xUnits), cex = key.cex, col = textCol, adj = c(0,0.5))
                                },
                                ...)
        
      }else{
        #plot spatially referenced x with single color and no key
        planform.plot <- xyplot(get(yCoord) ~ get(xCoord), data = point.data, col = point.col, pch = 19,
                                ylim = ylim, xlim = xlim, aspect = aspect, ylab = ylab, xlab = xlab, cex = arrowwd,
                                scales = scales,
                                ...)
      }
    }
  }
  
  #if interpolate is provided develop a plot of the interpolated data
  if(missing(interpolate) == FALSE){ 
    #if cellWidth was not provided, set it to 1/100 of the x range
    if(missing(cellWidth) == TRUE){
      cellWidth <- abs(diff(xlim))/100
    } else {
    }
    #if cellHeight was not provided set it to 1/100 of the y range
    if(missing(cellHeight) == TRUE){
      cellHeight <- abs(diff(ylim))/100
    }
    
    #create an object with x and y coordinates of sampled area, which will be used to interpolate interpolate to a grid
    coorADCP <- cbind(data[[xCoord]], data[[yCoord]])
    spADCP <- SpatialPointsDataFrame(coorADCP, data)
    
    
    #set grid for interpolation that is within the measured area
    bb <- bbox(spADCP)
    cs <- c(cellWidth, cellHeight)
    cc <- bb[,1] + (cs/2)
    cd <- ceiling(diff(t(bb))/cs)
    gridADCP <- GridTopology(cellcentre.offset = cc, cellsize = cs, cells.dim = cd)
    p4s <- CRS(proj4string(spADCP))
    sgADCP <- SpatialGrid(gridADCP, proj4string =p4s)
    
    #create an object for x to be evaluated in following functions
    interpolate <- as.character(substitute(interpolate))  
    interpolate <- eval(quote(interpolate))  
    #use inverse distance weighting to interpolate interpolate
    idwADCP <- krige(get(interpolate) ~ 1, spADCP[!is.na(spADCP@data[[interpolate]]),], sgADCP)
    
    #create polygon outlining sampled area, this will be used to remove interpolated data that falls outside the measured area
    suppressWarnings(data[,startDist := tDist - min(tDist, na.rm = TRUE), by = transectName])
    rightBank <- data[startDist == min(data$tDist, na.rm = TRUE), .SD, ]
    rightBank <- cbind(rightBank[[xCoord]], rightBank[[yCoord]])
    suppressWarnings(data[,inverseDist := tDist - max(tDist, na.rm = TRUE), by = transectName])
    leftBank <- data[inverseDist == max(data$inverseDist, na.rm = TRUE), .SD, ]
    leftBank <- cbind(leftBank[[xCoord]], leftBank[[yCoord]])
    if(abs(diff(range(rightBank[[1]]))) > abs(diff(range(rightBank[[2]])))){
      bound <- rbind(rightBank <- rightBank[order(rightBank[,1]), ], leftBank <- leftBank[order(leftBank[,1], decreasing = TRUE), ])
    } else {
      bound <- rbind(rightBank <- rightBank[order(rightBank[,2]), ], leftBank <- leftBank[order(leftBank[,2], decreasing = TRUE), ])
    }
    bound <- rbind(bound, bound[1,])
    bound <- Polygon(bound)
    bound <- Polygons(list(bound), 1)
    bound <- SpatialPolygons(list(bound))

    #clip the interpolated data set to the sampled area
    idwADCP <- idwADCP[!is.na(over(idwADCP, bound)),]    
    
    #create data.table of interpolated data
    idw.data <- data.frame(idwADCP)
    idw.data <- data.table(idw.data[,-(2)])
    setnames(idw.data, names(idw.data), c(interpolate, xCoord, yCoord))
    
    #create plot of interpolated data
    idwPlot <- levelplot(get(interpolate) ~ get(xCoord) * get(yCoord), data = idw.data,
                         xlim = xlim, ylim = ylim, aspect = aspect, ylab = ylab, xlab = xlab,
                         scales = scales,
                         col.regions = interpolate.color(length(unique(idw.data[[interpolate]]))),
                         at = seq(min(idw.data[[interpolate]], na.rm = TRUE), max(idw.data[[interpolate]], na.rm = TRUE), length = length(unique(na.omit(idw.data[[interpolate]])))),
                         colorkey = colorkey,
                         ylab.right = ylab.right,
                         par.settings = par.settings,
                         ...)
    #is x was provided lay the plot of x over the interpolated data
    if(missing(x)==FALSE){
      planform.plot <- idwPlot + as.layer(planform.plot, under = FALSE) 
    }else{      
      planform.plot <- idwPlot
    }
    #if contour is TRUE, add isopleths of interpolate to the plot
    if(contour==TRUE){
      planform.plot <- planform.plot +
        as.layer(contourplot(get(interpolate) ~ get(xCoord) * get(yCoord), data = idw.data, cuts = cuts, labels = labels, lwd=contourwd))
    } else {    
    }
  }else{
    
  }
  return(planform.plot)
}