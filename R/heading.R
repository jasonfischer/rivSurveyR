#' heading
#'
#' Calculates the heading in degrees of two velocities
#' @param ve	A vector of velocities to the east.
#' @param vn	A vector of velocities to the north. Must be the same length as ve
#' @return	A vector the same length as \code{ve}.
#' @export
#' @examples
#' heading(rnorm(5), rnorm(5))

heading <- function(ve, vn) {
  direction <- atan(abs(ve)/abs(vn)) * 180/pi
  for (i in seq_along(direction)) {
    switch(ve[i] > 0 & vn[i] < 0, direction[i] <- 180 - direction[i])
    switch(ve[i] <= 0 & vn[i] < 0, direction[i] <- direction[i] + 180)
    switch(ve[i] < 0 & vn[i] > 0,   direction[i] <- 90 - direction[i] + 270)
    switch(ve[i] < 0 & vn[i] == 0, direction[i] <- 270)  	
  }
  return(direction)
} 