#' heading
#'
#' Calculates the degree heading of two velocity vectors
#' @param ve	A vector of velocities to the east.
#' @param vn	A vector of velocities to the north. Must be the same length as ve
#' @return	A vector the same length as \code{ve}.
#' @export
#' @examples
#' heading(rnorm(5), rnorm(5))
#' easting <- c(0,1,1,1,0,-1,-1,-1)
#' northing <- c(1,1,0,-1,-1,-1,0,1)
#' heading(easting,northing)

heading <- function(ve, vn) {
  direction <- atan2(ve, vn) * 180/pi
  direction[which(direction < 0)] <- direction[which(direction < 0)] + 360
  return(direction)
}