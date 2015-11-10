#' black.white
#' 
#' Create a vector of n continuous colors from black to white
#' @param n	Number of colors in palette.
#' @return A	Character vector of colors
#' @export
#' @examples
#' plot(1:1000,1:1000, col=black.white(1000), pch=19)
#' @seealso	\link[grDevices]{colorRampPalette}

black.white <- colorRampPalette(c("black", "gray10","gray20", "gray30", "gray40", "gray50", "gray60", "gray70", "gray80", "gray90", "ivory"))
