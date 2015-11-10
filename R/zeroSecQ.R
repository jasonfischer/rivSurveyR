#' zeroSecQ
#' 
#' Uses the zero secondary discharge method of rotation to calculate secondary velocities within a cross-section
#' @param ve	Vector of east velocities within cells
#' @param vn	Vector of north velocities within cells. Must be the same length as \code{ve}.
#' @param vu	Vector of vertical velocities within cells. Must be the same length as \code{ve}.
#' @param cellHeight	Vector of cell heights. Must be the same length as \code{ve}.
#' @param cellWidth	Vector of cell Widths. Must be the same length as \code{ve}.
#' @param transectHeading	Compass heading of transect, assumes transect begins at river right and heads towards river left.
#' @return	Returns a \link[base]{data.frame} with within cell discharge perpendicular to cross-section (qx), within cell discharge parallel to cross-section (qy), within cell discharge in the primary flow direction (qp), within cell discharge in the secondary flow direction (qs), velocity in the primary direction (vp.zsq), velocity in the secondary direction (vs.zsq), compass heading of primary velocity (vpHeading.zsq), within cell compass heading of secondary velocity (vsHeading.zsq), and arithmetic heading of within cell secondary velocity (vsArithmeticHeading.zsq).
#' @export
#' @examples
#' data(cellVels)
#' velSub <- cellVels[transectName == "t1",,]
#' tHeading <- transectHeading(velSub$UTM_X_Proj, velSub$UTM_Y_Proj, velSub$Vel.E, velSub$Vel.N, velSub$cellDepth)
#' zsq <- zeroSecQ(velSub$Vel.E, velSub$Vel.N, velSub$Vel.up, velSub$cellHeight, velSub$cellWidth, tHeading)
#' @seealso	\link{xSec}.
#' @references	Lane, S.N., K.F. Bradbrook, K.S. Richards, P.M. Biron, and A.G. Roy. 2000. Secondary circulation cells in river channel confluences: Measurement artefacts or coherent flow structures? Hydrological Processes 14:2047-2071.
 
zeroSecQ <- function(ve, vn, vu, cellHeight, cellWidth, transectHeading){
  #pass data to xSec
  xSecVels <- xSec(ve, vn, vu, transectHeading)
  vx <- xSecVels$vx
  vy <- xSecVels$vy

  qy <- (vy * cellHeight * cellWidth) #cross-stream discharge in a cell
  qx <- (vx * cellHeight * cellWidth) #downstream discharge in a cell
  
  Qy <- sum(qy, na.rm = TRUE) #cross-stream discharge
  Qx <- sum(qx, na.rm = TRUE) #downstream discharge
  
  #deviation from the flow perpendicular to the cross-section, in geographic degrees
  primaryHeading.ZeroNetQ.Dev <- heading(Qx,Qy)

  #calculate primary and secondary discharge in each cell
  qp <- (qx * cos((90-primaryHeading.ZeroNetQ.Dev) * pi/180)) + (qy * sin((90-primaryHeading.ZeroNetQ.Dev) * pi/180))
  qs <- (-1 * qx * sin((90-primaryHeading.ZeroNetQ.Dev) * pi/180)) + (qy * cos((90-primaryHeading.ZeroNetQ.Dev) * pi/180))

  #calculate total primary and secondary discharge 
  Qp <- sum(qp, na.rm = TRUE)
  Qs <- sum(qs, na.rm = TRUE)
    
  #convert primary and secondary cell discharge into primary and secondary velocity
  vs.zsq <- qs/(cellHeight * cellWidth)
  vp.zsq <- qp/(cellHeight * cellWidth)
  
  #determine the compass heading of primary velocity
  vpHeading.zsq <- transectHeading + primaryHeading.ZeroNetQ.Dev
  switch(vpHeading.zsq >= 360, vpHeading.zsq <- vpHeading.zsq - 360)

  #determine the compass heading of secondary velocities 
  vsHeading.zsq <- rep(NA, times = length(vs.zsq))
  vsHeading.zsq[which(vs.zsq < 0)] <- vpHeading.zsq + 90
  vsHeading.zsq[which(vs.zsq >= 0)] <- vpHeading.zsq - 90
  
  vsArithmeticHeading.zsq <- heading(vu, vs.zsq)
  
  noSecQ <- data.frame(qx, qy, qp, qs, vp.zsq, vs.zsq, vpHeading.zsq, vsHeading.zsq, vsArithmeticHeading.zsq)
  
  return(noSecQ)
}