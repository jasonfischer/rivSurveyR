#' transectHeading
#'
#' Determines the compass heading of a cross-section transect, assuming transect starts at river right
#' @param x	A vector of x coordinates.
#' @param y	A vector of y coordinates. Must be the same length as \code{x}.
#' @param velE	A vector of velocities to the east. Must be the same length as \code{x}.
#' @param velN	A vector of velocities to the north. Must be the same length as \code{x}.
#' @param depth	A vector of water depths. Must be the same length as \code{x}.
#' @param flowHeading	Optional. Mean direction of flow over the cross-section. If provided velE, velN, and depth may be omitted.
#' @return	Returns the compass heading of a transect
#' @export
#' @examples
#' data(vels)
#' velSub <- vels[transectName == "t1",,]
#' transectHeading(velSub$UTM_X_Proj, velSub$UTM_Y_Proj, velSub$Mean.Vel.E, velSub$Mean.Vel.N, velSub$depth)
#' #If flow heading is known
#' transectHeading(velSub$UTM_X_Proj, velSub$UTM_Y_Proj, flowHeading = 149)

transectHeading <- function(x, y, velE, velN, depth, flowHeading){
  xy <- na.omit(data.table(x,y))  

  #determine the direction (x or y) of the greatest dispersion and then determine the transect heading based on the distance traveled in the x and y directions
  if (diff(range(xy$x)) >= diff(range(xy$y))){
    xy <- xy[order(x),,]
  }else{
    xy <- xy[order(y),,]
  }
  tHeading <- heading(xy$x[nrow(xy)]-xy$x[1], xy$y[nrow(xy)]-xy$y[1])

  if (missing(flowHeading) == TRUE){
    flowHeading <- heading(mean(velE * depth, na.rm = TRUE), mean(velN * depth, na.rm = TRUE))
  }else{
    flowHeading <- flowHeading
  }
  
  #use flow heading to determine river right and rotate transect heading 180 degrees if transect is heading from river left to right
  switch(tHeading <= 180 & (tHeading-flowHeading >=0 | tHeading-flowHeading < -180), tHeading <- tHeading + 180)
  switch(tHeading > 180 & tHeading-flowHeading > 0 & tHeading-flowHeading <= 180, tHeading <- tHeading + 180)
  switch(tHeading >= 360, tHeading <- tHeading - 360)
  
  return(tHeading)
}