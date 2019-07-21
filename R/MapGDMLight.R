#' Create a PCA transformed set of variables to map GDM
#' @description An internal function to create a PCA transformed set of variables to map spatially the results of GDM analysis using an RGB palette.
#'
#' @param rast.trans A RasterStack of gdm transformed environmental predictors returned by the gdm.transform function in gdm-package
#' @return A list with the following objects:
#' \describe{
#'   \item{pcaRast}{A RasterStack with the first three principal components}
#'   \item{pcaRast.all}{A RasterStack with the all the principal components}
#' }
#' @examples
#' gdm.map <- MapGDMLight(rast.trans)
#' plotRGB(gdm.map)

MapGDMLight<-function(rast.trans){
  #Maps general dissimilarity

  #Args:
  #   rastTrans: a rasterStack with GDM transformed environmental variables
  #Returns:
  #   pcaRast.all: a rasterStack with the pca transformed rastTrans
  #   pcaRast: a rasterStack of the first three pca transformed rastTrans, scaled from 0 to 255
  #            to plot as RGB.

  require(gdm)
  require(raster)
  #Plot dissimilarity
  cellsWData <- raster::Which(!is.na(rast.trans[[1]]),cells=T)
  rastDat <- rast.trans[cellsWData]
  pcaSamp <- prcomp(rastDat)
  pcaRast.all <- raster::predict(rast.trans, pcaSamp, index=1:nlayers(rast.trans))
  pcaRast <- pcaRast.all[[1:3]]

  # scale rasters
  pcaRast[[1]] <- (pcaRast[[1]]-cellStats(pcaRast[[1]],min)) /
    (cellStats(pcaRast[[1]],max)-cellStats(pcaRast[[1]],min))*255
  pcaRast[[2]] <- (pcaRast[[2]]-cellStats(pcaRast[[2]],min)) /
    (cellStats(pcaRast[[2]],max)-cellStats(pcaRast[[2]],min))*255
  pcaRast[[3]] <- (pcaRast[[3]]-cellStats(pcaRast[[3]],min)) /
    (cellStats(pcaRast[[3]],max)-cellStats(pcaRast[[3]],min))*255

  return(list(pcaRast=pcaRast,pcaRast.all=pcaRast.all))
}

