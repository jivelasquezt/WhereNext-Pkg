#' Remove correlated layers
#'
#' @description For a pair of correlated raster layers, this function will remove one. The output stack will have layers with a correlation between them below the user specified threshold.
#'
#' @param in.stack a RasterStack object
#' @param t.cor the maximum allowed correlation between layers in the RasterStack. Must be a value between 0 and 1.
#'
#' @return A RasterStack of layers with correlation below the user specified correlation threshold.
#'
#' @examples
#' env.vars.uncorr <- RemCorrLayers(env.vars, 0.8)
#' @export


RemCorrLayers <- function(in.stack, t.cor) {
  # Convert RasterBrick or RasterStack to SpatRaster if needed
  if (inherits(in.stack, "RasterBrick") || inherits(in.stack, "RasterStack")) {
    in.stack <- terra::rast(in.stack)
  }
  
  # Define sample size based on number of cells
  ssize <- if (ncell(in.stack[[1]]) < 10000) {
    ncell(in.stack[[1]])
  } else {
    10000
  }
  
  # Sample values with spatSample
  sample.vals <- terra::spatSample(in.stack, size = ssize, method = "random", na.rm = TRUE)
  
  # Calculate correlation matrix
  cor.matrix <- cor(sample.vals, use = "pairwise.complete.obs")
  
  i <- 1
  while (i <= ncol(cor.matrix)) {
    rem.ind <- which(abs(cor.matrix[, i]) > t.cor & abs(cor.matrix[, i]) < 1)
    if (length(rem.ind) == 0) {
      i <- i + 1
      next
    } else {
      print(paste("removing ", colnames(cor.matrix)[rem.ind]))
      cor.matrix <- cor.matrix[-rem.ind, -rem.ind]
      i <- i + 1
    }
  }
  
  # Select layers based on remaining columns in correlation matrix
  out.stack <- in.stack[[which(names(in.stack) %in% colnames(cor.matrix))]]
  return(out.stack)
}
