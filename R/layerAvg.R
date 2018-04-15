#' layerAvg
#' 
#' Calculates the integrated mean of y within the defined bounds of x
#' @param x	Sorted vector of the independent variable
#' @param y	A vector of the dependent variable
#' @param x1	Lower bound of x
#' @param x2  Upper bound of x
#' @param na.rm If TRUE NA values are removed before computation. Default is FALSE
#' @return  Returns the integrated mean of y
#' @export
#' @examples
#' x <- c(1:10)
#' y <- x^2
#' layerAvg(x, y, 2, 8)
 
layerAvg <- function(x, y, x1, x2, na.rm = FALSE){
  if(x2 < min(x)){
    layAvg <- NA
    warning("x2 < min(x), NA produced")    
  }else if(x1 > max(x)){
    layAvg <- NA
    warning("x1 > max(x), NA produced")
  }else{
    if(x1 < min(x, na.rm=TRUE)){
      x1 <- min(x, na.rm=TRUE)
    } else {
    }
    if(x2 > max(x, na.rm=TRUE)){
      x2 <- max(x, na.rm=TRUE)
    } else {
    }
    dt <- data.frame(x,y)
    if(na.rm==TRUE){
      dt <- na.omit(dt)
    } else {
    }
    r1 <- which.min(abs(dt$x - x1))
    r2 <- which.min(abs(dt$x - x2))
    if(length(r1)==0 |length(r2)==0){
      layAvg <- NA
    }else if(r1 == r2){
      layAvg <- dt$y[r1]
    } else {
      intg <- trapz(dt$x[r1:r2], dt$y[r1:r2])
      layAvg <- intg/abs(x2 - x1)      
    }
  }
  return(layAvg)
}