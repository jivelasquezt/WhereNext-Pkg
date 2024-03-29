#' Download GBIF data
#' @description This function retrieves GBIF data in SIMPLE_CSV for a user provided taxon key and area, the latter either as a ISO2
#' country code or a SpatialPolygons* object in WGS84 (EPSG: 4326) projection. When area is defined by a 
#' SpatialPolygons* object, records for the entire extent of the shapefile are retrieved.
#'
#' @param key a gbif taxon key obtained through the rgbif::name_backbone function. Required
#' @param user your GBIF account user name. Required.
#' @param user.email your GBIF account email. Required.
#' @param pwd your GBIF password. Required.
#' @param country ISO2 code of the coutry to download data from.
#' @param custom.shp a SpatialPolygons* object to download data from. Ignored if country is not NULL.
#' @return a list with:
#' \describe{
#'   \item{occ.table}{Occurrence table for the defined taxon key and area}
#'   \item{citation}{The citation of your GBIF download}
#' }
#' @note By default this function will sample at random 2 million records for files
#' with more than 2 milion rows. Change the default n.max and subset.size to change
#' this behavior.
#' @examples
#'
#' \dontrun{
#' library(raster)
#' library(rgbif)
#' library(WhereNext)
#'
#' #Example 1: Get occurrence data for colombian bats by country code
#' gbif.key <- name_backbone(name = "Chiroptera")
#' gbif.res <- DownloadGBIF(gbif.key$orderKey, "your username", "your email", "your password", "CO") #Enter your GBIF credentials here
#'
#' #Example2: Get occurrence data for bats from the Antioquia province in Colombia by using a shapefile
#' col.shp <- getData(name="GADM",country="CO", level=1) #Get country shapefile with level 1 administrative boundaries
#' antioquia <- col.shp[2, ]
#' gbif.res <- DownloadGBIF(gbif.key$orderKey, "your username", "your email", "your password", NULL, antioquia)
#' 
#' #Refine data to input shapefile
#' ind.nna<- !is.na(extract(antioquia, gbif.res$occ.table[,c("decimalLongitude", "decimalLatitude")])$poly.ID)
#' occ.table <- gbif.res$occ.table[ind.nna, ]
#' plot(antioquia)
#' points(occ.table$decimalLongitude, occ.table$decimalLatitude)
#' }
#' @export

DownloadGBIF <- function(key, user, user.email, pwd, country, custom.shp){
  if(!is.null(country)){
    if(country %in% ISOcodes::ISO_3166_1$Alpha_2){
      data.request <- rgbif::occ_download(rgbif::pred_in("taxonKey", key),
                                          rgbif::pred_in("country",country),
                                          rgbif::pred_in("hasCoordinate",TRUE),
                                          user=user,
                                          pwd=pwd,
                                          email=user.email,
                                          format="SIMPLE_CSV")
    } else {
      print(paste("Input country ISO2 code", country, "not valid"))
      return()
    }
  } else {
    if(is.null(custom.shp)){
      print(paste("Must provide either a valid country ISO2 code or SpatialPolygons* object"))
      return()
    } else {
      aoi.wkt <- rgbif::gbif_bbox2wkt(bbox = sp::bbox(custom.shp))
      data.request <- rgbif::occ_download(rgbif::pred_in("taxonKey", key),
                                          rgbif::pred_in("hasCoordinate",TRUE),
                                          rgbif::pred_within(aoi.wkt),
                                          user=user,
                                          pwd=pwd,
                                          email=user.email,
                                          format="SIMPLE_CSV")
    }
  }

  #Sleep until query returns total number of records or it's status changes to KILLED
  current.status <- rgbif::occ_download_meta(data.request)$status

  #Check every 30 seconds if data is ready
  while(current.status!="SUCCEEDED" & current.status!="KILLED"){
    s <- system.time(Sys.sleep(10))
    current.status <- rgbif::occ_download_meta(data.request)$status
    print(paste("Download status is:", current.status))
  }
  data.grab <- rgbif::occ_download_get(data.request, overwrite = T)
  citation <- renderPrint({rgbif::gbif_citation(data.grab)$download})
  occ.table <- GetGBIFData(data.request)
  return(list(occ.table=occ.table, citation=citation))
}

