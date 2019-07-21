#Modified from function viewRGB in mapview (Appelhans et al. 2019)
rgbPalette<-function(xout, r=1, g=2, b=3){
  quantiles = c(0.02, 0.98)
  mat <- cbind(xout[[r]][],
               xout[[g]][],
               xout[[b]][])
  for(i in seq(ncol(mat))){
    z <- mat[, i]
    lwr <- stats::quantile(z, quantiles[1], na.rm = TRUE)
    upr <- stats::quantile(z, quantiles[2], na.rm = TRUE)
    z <- (z - lwr) / (upr - lwr)
    z[z < 0] <- 0
    z[z > 1] <- 1
    mat[, i] <- z
  }
  na_indx <- apply(mat, 1, anyNA)
  cols <- mat[, 1]
  cols[na_indx] <- rgb(0, 0, 0, alpha=0)
  cols[!na_indx] <- rgb(mat[!na_indx, ], alpha = 1)
  p <- function(x) cols
  return(p)
}
