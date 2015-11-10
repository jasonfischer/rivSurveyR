#' ortho.proj
#' 
#' Calculates the x and y coordinates along an orthogonal plane, given the linear regression coefficients of y~x
#' @param x	A vector of X geographic coordinates
#' @param y	A vector of Y geographic coordinates
#' @param m	Slope parameter from linear regression of y~x
#' @param b	Intercept parameter from linear regression of y~x
#' @return	A \link[base]{data.frame} containing projected x and y coordinates (xProj and yProj). The number of rows equals the number of observations in \code{x}.
#' @export
#' @examples
#' data(vels)
#' vels.t1 <- vels[transectName=="t1",,]
#' plot(UTM_Y~UTM_X, vels.t1)
#' linmod <- lm(vels.t1$UTM_Y~vels.t1$UTM_X)
#' projected <- ortho.proj(vels.t1$UTM_X, vels.t1$UTM_Y, linmod$coefficients[2], linmod$coefficients[1])
#' plot(projected)

ortho.proj <- function(x, y, m, b){
  xProj <- (x - (m * b) + (m * y))/((m^2) + 1)
  yProj <- (b + (m * x) + ((m^2) * y))/((m^2) + 1)
  orthoProj <- data.frame(xProj, yProj)
  return(orthoProj)
}	