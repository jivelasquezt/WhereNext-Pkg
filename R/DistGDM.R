#' Compute the environmental or biological distance between pairs of sites
#'
#' @param site.grid A two column data.frame with the cell id of pairs of sites to compare. Cell ids must match the env.vars raster* object.
#' @param gdm.rast A gdm model object.
#' @param env.vars A rasterStack object of gdm transformed environmental predictors. Obtained by using the gdm.transform option in the gdm package.
#' @param dist Either "bio" to compute the biological distance between pairs of sites or "mhtn" to compute the manhattan environmental distance between pairs of sites.
#' @return A vector of biological distances between pairs of sites.
#'
#' @export
#' @import raster
#'
DistGDM <- function(site.grid, gdm.rast, env.vars, dist="mhtn"){
  colnames(site.grid) <- c("Var1", "Var2")
  if(dist=="bio"){
    res <- rep(gdm.rast$intercept, nrow(site.grid))
    for (i in 1:nlayers(env.vars)){
      res <- res+abs(env.vars[[i]][site.grid$Var1] - env.vars[[i]][site.grid$Var2])
    }
    res <- 1 - 1/(exp(res))
  }
  if(dist=="mhtn"){
    res <- rep(0, nrow(site.grid))
    for (i in 1:nlayers(env.vars)){
      res <- res+abs(env.vars[[i]][site.grid$Var1] - env.vars[[i]][site.grid$Var2])
    }
  }
  return(res)
}
