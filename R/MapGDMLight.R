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
  require(terra)
  #Plot dissimilarity
  cellsWData <- which(!is.na(terra::values(rast.trans[[1]])))
  rastDat <- rast.trans[cellsWData]
  pcaSamp <- prcomp(rastDat)
  pcaRast.all <- terra::predict(rast.trans, pcaSamp, index = 1:nlyr(rast.trans))
  pcaRast <- pcaRast.all[[1:3]]
  
  # scale rasters
  pcaRast[[1]] <- (pcaRast[[1]] - min(terra::values(pcaRast[[1]], na.rm = TRUE))) /
    (max(terra::values(pcaRast[[1]], na.rm = TRUE)) - min(terra::values(pcaRast[[1]], na.rm = TRUE))) * 255
  pcaRast[[2]] <- (pcaRast[[2]] - min(terra::values(pcaRast[[2]], na.rm = TRUE))) /
    (max(terra::values(pcaRast[[2]], na.rm = TRUE)) - min(terra::values(pcaRast[[2]], na.rm = TRUE))) * 255
  pcaRast[[3]] <- (pcaRast[[3]] - min(terra::values(pcaRast[[3]], na.rm = TRUE))) /
    (max(terra::values(pcaRast[[3]], na.rm = TRUE)) - min(terra::values(pcaRast[[3]], na.rm = TRUE))) * 255
  
  return(list(pcaRast=pcaRast,pcaRast.all=pcaRast.all))
}
