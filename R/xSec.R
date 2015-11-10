#' xSec
#' 
#' Rotates velocities to be parallel and perpendicular to a cross-section
#' @param ve	Vector of east velocities within cells
#' @param vn	Vector of north velocities within cells. Must be the same length as \code{ve}.
#' @param vu	Vector of vertical velocities within cells. Must be the same length as \code{ve}.
#' @param transectHeading	Compass heading of transect, assumes transect begins at river right and heads towards river left.
#' @return	Returns a \link[base]{data.frame} with compass heading of flow perpendicular to cross-section (vxHeading), compass heading of flow parallel to cross-section (vyCompassHeading), arithmetic heading of flow parallel to cross-section (vyArithmeticHeading), difference between flow perpendicular to cross-section and heading of flow within a cell (diffHeading), velocity perpendicular to cross-section (vx), and velocity parallel to cross-section where positive values indicate flow from river right to river left (vy).
#' @export
#' @examples
#' data(cellVels)
#' velSub <- cellVels[transectName == "t1",,]
#' tHeading <- transectHeading(velSub$UTM_X_Proj, velSub$UTM_Y_Proj, velSub$Vel.E, velSub$Vel.N, velSub$cellDepth)
#' secondary <- xSec(velSub$Vel.E, velSub$Vel.N, velSub$Vel.up, tHeading)
#' @seealso	\link{transectHeading}
 
xSec <- function(ve, vn, vu, transectHeading) {
  vxHeading <- transectHeading + 90 #assumes facing upstream
  switch(vxHeading < 0, vxHeading <- vxHeading + 360)
  switch(vxHeading >= 360, vxHeading <- vxHeading - 360)
  
  velHeading <- heading(ve, vn)
  
  #calculate the difference between flow heading perpendicular to transect and depth averaged flow heading
  diffHeading <- (vxHeading - velHeading)
  
  speed.xSec <- sqrt(ve^2 + vn^2 + vu^2)
  
  #if transect is oriented north and south, vx equals east velocity, if transect is oriented east and west, vx equals north velocity, otherwise vx equal cos(diffHeading) * water speed in cell
  if(transectHeading == 0 | transectHeading == 180){
    vx <- ve
    vy <- vn
  } else if(transectHeading == 90 | transectHeading == 270){
    vx <- vn
    vy <- ve
  } else {
    vx <- cos(diffHeading * pi/180) * speed.xSec
    vy <- sin(diffHeading * pi/180) * speed.xSec
  }
  
  #determine the compass heading of velocities parallel to transect
  vyCompassHeading <- vxHeading - 90
  switch(vyCompassHeading < 0, vyCompassHeading <- vyCompassHeading + 360)
  switch(vyCompassHeading > 360, vyCompassHeading <- vyCompassHeading - 360)
  
  vyArithmeticHeading <- heading(vu, vy)
  
  xSecVel <- data.frame(vxHeading, vyCompassHeading, vyArithmeticHeading, diffHeading, vx, vy)
    
  return(xSecVel)
}