#' white.black
#' 
#' Create a vector of n continuous colors from white to black
#' @param n	Number of colors in palette.
#' @return A	Character vector of colors
#' @export
#' @examples
#' plot(1:1000,1:1000, col=white.black(1000), pch=19)
#' @seealso	\link[grDevices]{colorRampPalette}

white.black <- colorRampPalette(c("ivory", "gray90","gray80", "gray70", "gray60", "gray50", "gray40", "gray30", "gray20", "gray10", "black"))
