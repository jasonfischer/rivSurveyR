#' shearVel
#' 
#' Calculates shear velocity from ADCP measurements
#' @param ve	A vector of velocities to the east, should be bottom cell velocity is reference is "bottomCell" and depth averaged velocity if reference is "mean"
#' @param vn	A vector of velocities to the north, should be bottom cell velocity is reference is "bottomCell" and depth averaged velocity if reference is "mean". Must be the same length as \code{ve}.
#' @param vu	A vector of vertical velocities, default is 0, should be bottom cell velocity is reference is "bottomCell" and depth averaged velocity if reference is "mean". Must be the same length as \code{ve}.
#' @param depth	A vector of depths, the same length as \code{ve}.
#' @param h	A vector containing the depth of the bottom cell, only needed if reference is bottomCell. Must be the same length as \code{ve}.
#' @param Dc	Sediment diameter of 90th percentile
#' @param n	Used to calculate local roughness height (Kc) Kc=n*Dc
#' @param p	Water density (kg/m^3)
#' @param reference	Method used to estimate shear velocity. If "mean" the Keulegan method is used where shear velocity = U/(1/k * ln(11* (depth/kc))), where U is the depth averaged speed, k is von Karman's constant (0.41), and, kc=n/Dc. For details see Garcia (2008). If "bottomCell" (the default) shear velocity = u/(9.5 * (h/kc)^(1/6)), where u is speed at h m above bottom. For details see Cheng-Lung (1991) and Simpson and Oltmann (1990).
#' @param units	Unit system of velocities, depths, and sediment diameters. Can either be "metric" (default) or "standard".
#' @return	A \link[data.table]{data.table} containing shear velocity to the east (ustarE), shear velocity to the north (ustarN), shear speed (ustar), and shear speed heading (ustarHeading). The number of rows equals the number of observations in \code{ve}.
#' @export
#' @examples
#' data(vels)
#' #Use bottom cell velocities
#' shearVel(vels$BC.Vel.E, vels$BC.Vel.N, vels$BC.Vel.Up, vels$depth, vels$bottomCellDepth, Dc = 0.2032, n = 2)
#' #Use depth averaged velocities
#' shearVel(vels$Mean.Vel.E, vels$Mean.Vel.N, depth = vels$depth, Dc = 0.2032, n = 2, reference = "mean")
#' @references	Garcia, M.H. 2008. Sediment transport and morphodynamics, chap. 2. In: Garcia, M.H. (Ed.), Sedimentation Engineering: Processes, Measurments, Modeling, and Practice No. 110. American Society of Civil Engineers, Reston, Virginia, pp. 21-163.
#' Cheng-Lung 1991. “Unified Theory on Power Laws for Flow Resistance.” Journal of Hydraulic Engineering, Vol. 117, No. 3, March 1991, 371-389. 
#' Simpson, M.R. and Oltmann, R.N. 1990. “An Acoustic Doppler Discharge Measurement System.” Proceedings of the 1990 National Conference on Hydraulic Engineering, Vol. 2, 903-908.

shearVel <- function(ve, vn, vu = 0, depth, h, Dc, n, p = 1000, reference = "bottomCell", units = "metric") {
  stopifnot(reference %in% c("bottomCell", "mean"), units %in% c("metric", "standard"))
  #convert standard units to metric
  if(units == "standard"){
    ve <- ve * 0.3048
    vn <- vn * 0.3048
    vu <- vu * 0.3048
    depth <- depth * 0.3048
    h <- h * 0.3048
    Dc <- Dc * 0.3048
  } else {
    
  }
  #calculate intermediary values for shear stress calculation	
  kc <- n * Dc
  Cf <- ((1/0.41) * log(11 * (depth/kc)))^-2
  
  if (reference == "mean") {
    TbE <- p * Cf * ve^2
    ustarE <- sqrt(TbE / p) * ve/abs(ve)
    TbN <- p * Cf * vn^2
    ustarN <- sqrt(TbN / p) * vn/abs(vn)
    TbU <- p * Cf * vu^2
    ustarU <- sqrt(TbU / p) * vn/abs(vu)
  } else if(reference == "bottomCell") {
    h <- depth - h
    ustarE <- ve/(9.5 * (h/kc)^(1/6))
    ustarN <- vn/(9.5 * (h/kc)^(1/6))
    ustarU <- vu/(9.5 * (h/kc)^(1/6))
  } else {
    stop ("'reference' not found")
  }
  
  ustar <- sqrt(ustarE^2 + ustarN^2 + ustarU^2)
  ustarHeading <- heading(ustarE, ustarN)

  shearVel <- data.table(ustarE, ustarN, ustar, ustarHeading)
  if(units == "standard"){
    shearVel[,-4] <- shearVel[,-4] * 3.28084
  } else {
    
  }
  return(shearVel)
}