#' cellVels
#'
#' Dataset of processed cross-section data collected with a SonTek m9 ADCP
#' @details  \itemize{
#' \item \code{transectName}	Transect velocities were measured along
#' \item \code{Vel.E}	Cell velocity to the east in meters/second
#' \item \code{Vel.N}	Cell velocity to the north in meters/second
#' \item \code{Vel.up}	Cell vertical velocity in meters/second
#' \item \code{Mean.vel.E}	Depth averaged velocity to the east in meters/second
#' \item \code{Mean.Vel.N}	Depth averaged velocity to the north in meters/second
#' \item \code{cellDepth}	Depth of cell in water column in meters
#' \item \code{cellHeight}	Height of cell in meters
#' \item \code{cellWidth}	Width of cell in meters
#' \item \code{tDist}	Distance along transect in meters
#' \item \code{UTM_X_Proj}	Projected Universal Transverse Mercator (UTM) east coordinates 
#' \item \code{UTM_Y_Proj}	Projected UTM north coordinates
#' }
#'
#' @docType data
#' @keywords dataset
#' @name cellVels
#' @usage	\code{data(cellVels)}
#' @format	\code{cellVels} is a \link[data.table]{data.table} of 26383 rows and 12 columns, each row corresponds to a cell average. Although data were collected with an m9, they have been edited to represent fictitious values and locations.
NULL