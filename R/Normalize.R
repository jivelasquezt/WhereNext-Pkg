#' Normalize raster data
#' @description This function normalizes Raster* data to have a mean 0 and standard deviation 1.

#' @param in.raster a Raster* object
#' @return a standardized Raster* object
#' @examples
#' env.var.norm <- Normalize(env.vars)
#' @export

Normalize <- function(in.raster){
  # Calculate mean and standard deviation using terra functions
  mean_val <- global(in.raster, fun = "mean", na.rm = TRUE)
  sd_val <- global(in.raster, fun = "sd", na.rm = TRUE)
  
  # Convert mean and sd values to the same format as the raster (single value)
  mean_val <- mean_val$mean
  sd_val <- sd_val$sd
  
  # Apply normalization to the raster
  out.raster <- (in.raster - mean_val) / sd_val
  return(out.raster)
}

