#' Normalize raster data
#' @description This function normalizes Raster* data to have a mean 0 and standard deviation 1.

#' @param in.raster a Raster* object
#' @return a standardized Raster* object
#' @examples
#' env.var.norm <- Normalize(env.vars)
#' @export

Normalize <- function(in.raster){
  out.raster <- (in.raster - cellStats(in.raster,stat="mean"))/cellStats(in.raster, stat="sd")
}
