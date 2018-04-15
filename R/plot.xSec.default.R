#' plot.xSec.default
#' 
#' Plots cross-section velocities given an object of class "adcp.secondary"
#' @param data	An object of class "adcp.secondary"
#' @param x	Parameter to be plotted
#' @param transect	Transect to be plotted
#' @param cellWidth	Width of cell size used to grid x. If missing, defaults to mean cell width from dataset
#' @param cellHeight	Height of cell size used to grid x. If missing, defaults to mean cell height from dataset
#' @param arrows	If TRUE (default) and x is a primary velocity vector, arrows indicating secondary flow velocity are plotted
#' @param arrow.x.spacing	Horizontal spacing between secondary flow arrows. Units equal to units of x-axis (i.e., distance from start of transect).
#' @param arrow.y.spacing	Vertical spacing secondary flow arrows. Units equal to units of y-axis (i.e., depth).
#' @param arrow_scale	Multiplication factor for sizing secondary flow arrows, default is 1
#' @param arrowCol	Color of secondary flow arrows, default is "black"
#' @param arrowwd	Line width of arrows
#' @param lwd	Line width of surface and bottom outlines
#' @param arrow.key	Numeric vector of length 2 defining the X and y coordinates (c(x,y)) of where scaling key for secondary flow arrows should be written. If omitted, a key is placed in the lower left.
#' @param keyCol	Color of text and arrows in arrow.key, default is "white"
#' @param lineCol	Line color of surface and bottom outlines, default is white
#' @param col.regions	Color function used to display x, see \link[lattice]{levelplot} for details. Default is \link{jet.colors}.
#' @param scale.by	Character string. If "individual" (default) variables in plot are scaled to their range in the cross-section plotted. If "all" variables in plot are scaled to their range in the dataset provided. See 'details' for additional information
#' @param xlab	X-axis label, see \link[lattice]{levelplot} for details. Default is "Distance (m)" in size 18 font.
#' @param ylab	Y-axis label, see \link[lattice]{levelplot} for details. Default is "Depth (m)" in size 18 font.
#' @param ylab.right	\code{colorkey} label, see \link[lattice]{xyplot} for details. Default is "Velocity m/s" in size 18 font.
#' @param header	Title of plot, default is transect
#' @param headersize	Font size of header, default is 18
#' @param colorkey	A list of arguments for the color key drawn alongside the plot, see \link[lattice]{levelplot} for details. Default writes the tick labels with cex 1.5.
#' @param scales	A list of arguments determining how the axes are drawn, see \link[lattice]{xyplot} for details. Default is to draw tick labels with cex 1.5.
#' @param key.cex	Font size of text in arrow.key, default is 1.25
#' @param xlim	Numeric vector of length 2 defining the range of the x-axis.
#' @param ylim	Numeric vector of length 2 defining the range of the y-axis. See 'details' for additional information.
#' @param aspect	Controls the aspect ratio of the panels, see \link[lattice]{xyplot} for details. Default is "iso".
#' @param par.settings	A list of arguments for fine-tuned control of display, see \link[lattice]{xyplot} and \link[lattice]{trellis.par.set} for details. Default is not to pad the color key, offset the color key label by 2, and draw a black background.
#' @param ...	Additional arguments to be passed to \link[lattice]{levelplot}.
#' @param contour	Logical, if TRUE (default), contour lines are drawn
#' @param cuts	Number of regions contour lines denote
#' @param labels	Logical, if TRUE (default) contour lines are labeled
#' @details	When producing multiple plots of cross-sections, setting scale.by to "all" eases comparison among plots by using the same scale for all plots.  Setting scale.by to "individual" exaggerates differences within individual cross-sections. 
#' If scale.by equals "individual" default for \code{ylim} is the max depth of the cross-section to -0.5. If scale.by equals "all" default for \code{ylim} is the max depth of the dataset to -0.5.
#' Use transect argument to subset data by transect, may produce error if data is subset prior to being supplied to function.
#' @return	A levelplot plot of the cross-section
#' @export
#' @examples
#' data(mNine)
#' #mNine is a list of MATLAB files
#' names(mNine)
#' #Drop the last letter (l or r) from the MATLAB file names
#' tNames <- substr(names(mNine), 0, nchar(names(mNine))-1)
#' adcp.secondary <- process.secondary(mNine, tNames, rotation = c("xSec", "rozovskii", "zeroSecQ"), xWindow = 21, yWindow = 5)
#' plot.xSec(adcp.secondary, vp.zsq, transect = "t1")
#' #Use the same scale to display primary and secondary velocities of two cross-sections
#' plot.xSec(adcp.secondary, vp.zsq, transect = "t1", scale.by = "all")
#' @seealso \link{plot.xSec}, \link[lattice]{levelplot}, \link{process.secondary}, and \link[lattice]{xyplot}.
 
plot.xSec.default <- function(data, x, transect, cellWidth, cellHeight, arrows = TRUE, arrow.x.spacing, arrow.y.spacing,
                              arrow_scale = 1, arrowCol = "black", arrowwd = 1, lwd = 1, 
                              arrow.key, keyCol = "white", 
                              lineCol = "white", col.regions = jet.colors,                              
                              scale.by = "individual",
                              xlab = list(label = "Distance (m)", fontsize = 18), 
                              ylab = list(label = "Depth (m)", fontsize = 18),
                              ylab.right = list(label = "Velocity (m/s)", fontsize = 18),
                              header = transect, headersize = 18,
                              colorkey=list(labels=list(cex=1.5)), scales = list(cex = c(1.5, 1.5)),                              
                              key.cex = 1.25,
                              xlim, ylim, 
                              aspect = "iso", 
                              par.settings = list(layout.widths = list(axis.key.padding = 0, ylab.right = 2), 
                                                  panel.background = list(col = "black")),
                              ..., contour = TRUE, cuts = 7, labels = TRUE){
  stopifnot("adcp.secondary" %in% class(data))
  
  #subset data to only include transect of interest
  dt <- data[transectName == transect,,]
  main = list(label = as.character(header), fontsize = headersize)
  #create an object for x to be evaluated in following functions
  x <- as.character(substitute(x))
  x <- eval(quote(x)) 
 
  #if cellWidth is missing, set cellWidth to be equal to the mean width of cells in the transect 
  if(missing(cellWidth) == TRUE){
    cellWidth <- mean(dt$cellWidth, na.rm = TRUE)
    message("Cell Width: ", round(cellWidth, digits = 2))
  } else {
    
  }
  #if cellHeight is missing, set cellHeight to be equal to the mean height of cells in the transect   
  if(missing(cellHeight) == TRUE){
    cellHeight <- mean(dt$cellHeight, na.rm = TRUE)
    message("Cell Height: ", round(cellHeight, digits = 2))
  } else {
    
  }
  #if arrow.x.spacing is missing, set arrow.x.spacing to be equal to 1/20 the length of the transect 
  if(missing(arrow.x.spacing) == TRUE){
    arrow.x.spacing <- max(dt$tDist, na.rm = TRUE)/20
    message("Arrow X Spacing : ", round(arrow.x.spacing, digits = 2))
  } else {
    
  }  
  #if arrow.y.spacing is missing, set arrow.y.spacing to be equal to 1/20 the maximum depth measured along the transect 
  if(missing(arrow.y.spacing) == TRUE){
    arrow.y.spacing <- max(dt$depth, na.rm = TRUE)/20
    message("Arrow Y Spacing : ", round(arrow.y.spacing, digits = 2))
  } else {
    
  }
  
  #create an object with depths and distances along transect, which will be used to interpolate x to a grid
  dtInterp <- dt[!is.na(tDist) & !is.na(cellDepth),,]
  coorADCP <- cbind(tDist = dtInterp$tDist, cellDepth = dtInterp$cellDepth)
  
  #create a spatial points data.frame
  spADCP <- SpatialPointsDataFrame(coorADCP, dtInterp) 
  spoints <- SpatialPoints(coorADCP + (c(cellWidth, cellHeight)/2))
  
  #set grid for interpolation that is within the measured area (includes cells, below max depths)
  bb <- bbox(spADCP)
  cs <- c(cellWidth, cellHeight)
  cc <- bb[,1] + (cs/2)
  cd <- ceiling(diff(t(bb))/cs)
  gridADCP <- GridTopology(cellcentre.offset = cc, cellsize = cs, cells.dim = cd)
  p4s <- CRS(proj4string(spADCP))
  sgADCP <- SpatialGrid(gridADCP, proj4string = p4s)
  
  #use inverse distance weighting to interpolate x
  idwADCP <- krige(get(x) ~ 1, spADCP[!is.na(spADCP@data[,x]),], sgADCP)
  
  #create polygon outlining sampled area, this will be used to remove interpolated data that falls outside the measured area
  shallow <- dt[,min(cellDepth, na.rm = TRUE), by = tDist]
  dep <- rbind(dt[,mean(depth, na.rm = TRUE), by = tDist], dt[,mean(bottomCellDepth + mean(cellHeight, na.rm = TRUE)/2, na.rm = TRUE), by = tDist])
  dep <- dep[,min(V1), by=tDist] #keep whichever is less, vb depth or bottom cell depth
  bound <- rbind(dep[order(tDist)], shallow[order(tDist, decreasing = TRUE)])
  bound <- na.omit(bound)
  bound <- rbind(bound, bound[1,])
  bound <- Polygon(bound)
  bound <- Polygons(list(bound), 1)
  bound <- SpatialPolygons(list(bound))
  
  #clip the interpolated data set to the sampled area
  idwADCP <- idwADCP[!is.na(over(idwADCP, bound)),]
  
  #create a data.table of data to be plotted
  plot.data <- data.frame(idwADCP)
  plot.data <- data.table(plot.data[,-(2)])
  
  #change the column in plot.data containing x data to x
  setnames(plot.data, names(plot.data)[1], x)
  #create a column of distances along transect rounded down to nearest integer and make it the data.table key
  plot.data$tDist2 <- floor(plot.data$tDist)
  setkey(plot.data, tDist2)
  
  #create a data.table of depths at each ensemble, this will be used to define the bottom in the plot
  vb.data <- dt[,mean(depth, na.rm = TRUE), by = tDist]
  vb.data <- vb.data[order(vb.data[,tDist,])]
  setnames(vb.data, c("V1", "tDist"), c("depth", "tDist2"))
  setkey(vb.data, tDist2)
  plot.data <- merge(plot.data, vb.data) 
  setkey(plot.data, tDist)
  
  #set the x and y limits of the plot
  if("ylim" %in% (names(match.call(expand.dots=F)))){
    ylim = ylim
  } else {
    stopifnot(scale.by %in% c("individual", "all"))
    if(scale.by == "individual"){
      ylim = c(max(dt$depth, na.rm = TRUE)+diff(c(-0.5,max(dt$depth, na.rm = TRUE)))*key.cex*0.2, -0.5) #scales y axis to max depth in plot        
    } else {
      ylim = c(max(data$depth, na.rm = TRUE)+diff(c(-0.5,max(data$depth, na.rm = TRUE)))*key.cex*0.2, -0.5) #scales y axis to max depth in data set   
    }
  }
  if("xlim" %in% (names(match.call(expand.dots=F)))){
    xlim = xlim 
  } else {    
    xlim = c(min(plot.data$tDist, na.rm = TRUE), max(plot.data$tDist, na.rm = TRUE))      
  }
  
  #if arrows is true and a primary velocity parameter is defined for x, create an object used to call the appropriate secondary velocity parameter
  if(arrows == TRUE){
    if(x %in% c("vp.roz", "vpx.roz", "vsx.roz", "vp.zsq", "vx")){
      if(x == "vp.roz"){
        vs <- "vs.roz"
      } else if(x == "vpx.roz"){
        vs <- "vpy.roz"
      } else if(x == "vsx.roz"){
        vs <- "vsy.roz"
      } else if(x == "vp.zsq"){
        vs <- "vs.zsq"
      } else {
        vs <- "vy"
      }
      vs <- eval(quote(vs))
      
      #determine location for the arrow key
      if("arrow.key" %in% (names(match.call(expand.dots=F))) == FALSE){
        arrow.key <- c(min(xlim), max(ylim))
        key.x <- min(xlim)+(diff(xlim)*0.025)
        key.y <- c((max(ylim)+(diff(ylim)*0.14*key.cex)), (max(ylim)+(diff(ylim)*0.08*key.cex)), max(ylim)+(diff(ylim)*0.02*key.cex))
      } else if(is.numeric(arrow.key) == TRUE){
        key.x <- arrow.key[1]+diff(xlim)*0.05
        key.y <- c((arrow.key[2] + (diff(ylim)*0.14*key.cex)), (arrow.key[2]+(diff(ylim)*0.08*key.cex)), arrow.key[2]+(diff(ylim)*0.02*key.cex))
      } else {
      }
      
      #determine the values to be displayed in the arrow key
      if(scale.by == "individual"){
        keyValues <- summary(dt[,abs(na.omit(get(vs)))]) #This scales arrow size to individual transects
      } else {
        keyValues <- summary(data[,abs(na.omit(get(vs)))]) #This scales arrow size to all transects of interest similarly
      } 
      keyValues <- keyValues[which(!names(keyValues) == "Mean")]
      keyValues <- c(keyValues["Min."][[1]], keyValues["Median"][[1]], keyValues["Max."][[1]])
      keyValues <- round(as.numeric(keyValues), digits = 2)
      keyValues <- as.factor(keyValues)
      keyValues <- c(min(as.numeric(as.character(keyValues))), median(as.numeric(as.character(keyValues))), max(as.numeric(as.character(keyValues))))
      
      #set a grid and interpolate the horizontal and vertical secondary velocities
      spADCP <- SpatialPointsDataFrame(coorADCP, dtInterp) 
      bb <- bbox(spADCP)
      cs <- c(arrow.x.spacing, arrow.y.spacing)
      cc <- bb[,1] + (cs/2)
      cd <- ceiling(diff(t(bb))/cs)
      gridADCP <- GridTopology(cellcentre.offset = cc, cellsize = cs, cells.dim = cd)
      p4s <- CRS(proj4string(spADCP))
      sgADCP <- SpatialGrid(gridADCP, proj4string = p4s)
      
      x.arrow <- krige(get(vs) ~ 1, spADCP[!is.na(spADCP@data[,vs]),], sgADCP)    
      x.arrow <- x.arrow[!is.na(over(x.arrow, bound)),]    
      x.arrow <- data.frame(x.arrow) 
      x.arrow <- data.table(x.arrow[,(-2)])
      keycols <- c("tDist", "cellDepth")
      setkeyv(x.arrow, keycols)
      
      y.arrow <- krige(Vel.up ~ 1, spADCP[!is.na(spADCP$Vel.up),], sgADCP)
      y.arrow <- y.arrow[!is.na(over(y.arrow, bound)),]    
      y.arrow <- data.frame(y.arrow)
      y.arrow <- data.table(y.arrow[,(-2)])
      setkeyv(y.arrow, keycols)
      
      #merge the horizontal and vertical secondary velocities into one data.table
      flow.arrows <- (merge(x.arrow, y.arrow))
      setnames(flow.arrows, names(flow.arrows)[3:4], c("vs", "vu"))
      
      #scale arrow x and y vectors to the x and y range of the figure
      flow.arrows$scaleFactor <- diff(range(dt$tDist, na.rm = TRUE))/diff(range(dt$cellDepth, na.rm = TRUE))
      flow.arrows$arrow_scale <- arrow_scale
      if(is.numeric(arrow.key) == TRUE){
        key.scale <- data.frame(key.x, key.y, keyValues, scaleFactor = flow.arrows$scaleFactor[1], arrow_scale, keyCol, xlim = diff(xlim), ylim = diff(ylim))
      }else{
        
      }
    } else {
      warning("x must be vp.roz, vpx.roz, vsx.roz, vp.zsq, or vx to plot directional arrows")
    }   
  } else {
  }
  #define the color scale used to display x
  if(scale.by == "individual"){
    color.regions = col.regions(length(unique(dt[[x]]))) #scales everything to min and max of transect
    at.interval = seq(min(dt[[x]], na.rm = TRUE), max(dt[[x]], na.rm = TRUE), length = length(unique(na.omit(dt[[x]])))) #scales everything to min and max of transect
  } else {
    color.regions = col.regions(length(unique(data[[x]]))) #scales everything to min and max of dataset instead of transects
    at.interval = seq(min(data[[x]], na.rm = TRUE), max(data[[x]], na.rm = TRUE), length = length(unique(na.omit(data[[x]])))) #scales everything to min and max of dataset instead of transects
  }  
  
  #create plot of the data
  xSec.plot <- levelplot(get(x) ~ tDist * cellDepth, data = plot.data,
                         xlab = xlab, ylab = ylab, main = main, ylim = ylim, xlim = xlim,
                         col.regions = color.regions, 
                         at = at.interval, 
                         colorkey = colorkey,
                         scales = scales, 
                         panel = function(...){
                           panel.levelplot(...)
                           panel.abline(h = 0, col = lineCol, lwd = lwd, ...) #add line displaying water surface
                           if(exists("flow.arrows") == TRUE){ #if flow.arrows has been defined, add them and their key to the plot                          
                             panel.arrows(x0 = flow.arrows$tDist, y0 = flow.arrows$cellDepth,
                                          x1 = (flow.arrows$tDist + (flow.arrows$scaleFactor * flow.arrows$vs * flow.arrows$arrow_scale)), y1 = (flow.arrows$cellDepth - (flow.arrows$arrow_scale * flow.arrows$vu)),
                                          length = arrowwd*0.05, col = arrowCol, lwd = arrowwd)
                             if(is.numeric(arrow.key) == TRUE){
                               panel.arrows(x0 = key.scale$key.x, y0 = key.scale$key.y, x1 = (key.scale$key.x + (key.scale$scaleFactor * key.scale$keyValues * key.scale$arrow_scale)), y1 = key.scale$key.y, length = arrowwd*0.05,
                                            col = keyCol, lwd = arrowwd)
                               panel.text(x = key.scale$key.x ,y = key.scale$key.y + key.scale$ylim*0.02*key.cex, label = paste(key.scale$keyValues, "m/s"), cex= key.cex, col = keyCol, adj = c(0,0.5))
                             } else {                                  
                             }
                           } else {                                
                           }
                         }, 
                         ylab.right = ylab.right, 
                         par.settings = par.settings,
                         ...
  )
   #add bottom contour to the plot
   xSec.plot <- xSec.plot + 
    as.layer(xyplot(depth~tDist2, data = vb.data, type = 'l', col = lineCol,  lwd = lwd, ...), under = FALSE)
  #if contour is TRUE, add isopleths of x to the plot
  if(contour==TRUE){
    xSec.plot <- xSec.plot +
      as.layer(contourplot(get(x) ~ tDist * cellDepth, data = plot.data, cuts = cuts, labels = labels, lwd = lwd))
  } else {    
  }
  
  return(xSec.plot)
}  