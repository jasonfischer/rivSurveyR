#' jet.colors
#' 
#' Create a vector of n continuous colors from blue to red
#' @param n	Number of colors in palette.
#' @return A	Character vector of colors
#' @export
#' @examples
#' plot(1:1000,1:1000, col=jet.colors(1000), pch=19)
#' @seealso	\link[grDevices]{colorRampPalette}

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
