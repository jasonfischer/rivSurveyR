#' vels
#'
#' Dataset of processed planform data collected with a SonTek m9 ADCP
#' @details  \itemize{
#' \item \code{transectName}	Transect velocities were measured along
#' \item \code{UTM_X}	Universal Transverse Mercator (UTM) east coordinates
#' \item \code{UTM_Y}	UTM north coordinates
#' \item \code{Mean.vel.E}	Depth averaged velocity to the east in meters/second
#' \item \code{Mean.Vel.N}	Depth averaged velocity to the north in meters/second
#' \item \code{BC.Vel.E}	Velocity to the east in the deepest cell measured in a vertical in meters/second
#' \item \code{BC.Vel.N}	Velocity to the north in the deepest cell measured in a vertical in meters/second
#' \item \code{BC.Vel.up}	Vertical velocity in the deepest cell measured in a vertical in meters/second
#' \item \code{depth}	Depth of sample in meters
#' \item \code{bottomCellDepth}	Depth of deepest cell measured in meters
#' \item \code{tDist}	Distance along transect in meters
#' \item \code{UTM_X_Proj}	Projected UTM east coordinates
#' \item \code{UTM_Y_Proj}	Projected UTM north coordinates
#' }
#'
#' @docType data
#' @keywords dataset
#' @name vels
#' @usage	\code{data(vels)}
#' @format	\code{vels} is a \link[data.table]{data.table} of 1957 rows and 13 columns, each row corresponds to a sample average. Although data were collected with an m9, they have been edited to represent fictitious values and locations.
NULL