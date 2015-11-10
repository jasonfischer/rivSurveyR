#' rozovskii
#' 
#' Uses the Rozovskii method of rotation to calculate secondary velocities within a cross-section
#' @param ve	Vector of east velocities within cells
#' @param vn	Vector of north velocities within cells. Must be the same length as \code{ve}.
#' @param vu	Vector of vertical velocities within cells. Must be the same length as \code{ve}.
#' @param meanVe	Vector of depth averaged east velocities of the ensemble a cell belongs to. Must be the same length as \code{ve}.
#' @param meanVn	Vector of depth averaged north velocities of the ensemble a cell belongs to. Must be the same length as \code{ve}.
#' @return	Returns a \link[base]{data.frame} with heading of primary velocity (VpHeading.roz), compass heading of secondary velocity (vsCompassHeading.roz), arithmetic heading of secondary velocity (vsArithmeticHeading.roz), difference between cell heading and primary velocity heading (Theta.roz), primary velocity (vp.roz), secondary velocity (vs.roz), cross-stream component of primary velocity (vpy.roz), down-stream component of primary velocity (vpx.roz), cross-stream component secondary velocity (vsy.roz), and down-stream component of secondary velocity (vsx.roz)
#' @export
#' @examples
#' data(cellVels)
#' velSub <- cellVels[transectName == "t1",,]
#' roz <- rozovskii(velSub$Vel.E, velSub$Vel.N, velSub$Vel.up, velSub$Mean.Vel.E, velSub$Mean.Vel.N)
#' @references	Rhoads, B.L. and S.T. Kenworthy. 1998. Time-averaged flow structure in the central region of a stream confluence. Earth Surface Processes and Landforms 23:171-191.

rozovskii <- function(ve, vn, vu, meanVe, meanVn) {
  #calculate water speed in each cell
  Speed.roz=sqrt(ve^2 + vn^2 + vu^2)
  #determine heading of flow within each cell and depth averaged flow heading
  VelHeading.roz <- heading(ve, vn)
  vpHeading.roz <- heading(meanVe, meanVn)
  #calculate the difference between cell heading and depth averaged flow heading
  Theta.roz <- vpHeading.roz - VelHeading.roz
  #vp = primary vel, vs = secondary vel, vpu = primary up vel, vsu = secondary up vel, vpy = primary vel component parallel to cross section, vpx = primary vel component orthogonal to cross section, vsy = secondary vel component parallel to cross section, vsx = secondary vel component orthogonal to cross section
  vp.roz <- cos(Theta.roz * pi/180) * Speed.roz
  vs.roz <- sin(Theta.roz * pi/180) * Speed.roz #When referenced to the transect, positive values assumed to be to the right, but the inverse is true when referenced by geographic heading
  vpy.roz <- vp.roz * sin(vpHeading.roz * pi/180)
  vpx.roz <- vp.roz * cos(vpHeading.roz * pi/180)
  vsy.roz <- vs.roz * cos(vpHeading.roz * pi/180)
  vsx.roz <- vs.roz * sin(vpHeading.roz * pi/180)
  
  #determine the compass heading of secondary velocities in compass degrees
  vsCompassHeading.roz <- vpHeading.roz - 90
  #if secondary velocities are negative, rotate the heading 180 degrees
  vsCompassHeading.roz[which(vs.roz < 0)] <- vsCompassHeading.roz[which(vs.roz < 0)] + 180
  vsCompassHeading.roz[which(vsCompassHeading.roz < 0)] <- vsCompassHeading.roz[which(vsCompassHeading.roz < 0)] + 360
  vsCompassHeading.roz[which(vsCompassHeading.roz > 360)] <- vsCompassHeading.roz[which(vsCompassHeading.roz > 360)] - 360
  
  vsArithmeticHeading.roz <- heading(vu, vs.roz)
  
  rozVel <- data.frame(vpHeading.roz, vsCompassHeading.roz, vsArithmeticHeading.roz, Theta.roz, vp.roz, vs.roz, vpy.roz, vpx.roz, vsy.roz, vsx.roz)
  
  return(rozVel)
}