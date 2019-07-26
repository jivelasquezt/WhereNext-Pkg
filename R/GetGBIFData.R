#' Get GBIF data
#' @description This is an internal function to retrieve compressed (.zip) downloaded GBIF data in SIMPLE_CSV based on a occ_download object from the rgbif package. It's similar to rgbif:::occ_download_import but allows random sampling of rows and is more stable (so far).

#' @param data.request an occ_download object returned by the rgbif occ_download function
#' @return a data.frame of downloaded gbif data

GetGBIFData <- function(data.request){
  key <- rgbif::occ_download_meta(data.request)$key
  filepath <- paste0(key, ".zip")
  if(file.exists(filepath)){
    unzip(filepath)
    n.recs <- bigreadr::nlines(paste0(key,".csv"))
    print(paste("occurrence file has", n.recs, "records"))
    occs<-bigreadr::big_fread1(paste0(key,".csv"), 
                      every_nlines=1000000, 
                      .transform = function(x) {
                        res<-dplyr::select(x, c("gbifID","taxonRank", "species", "decimalLongitude" ,"decimalLatitude", "eventDate","countryCode","locality","recordedBy","coordinateUncertaintyInMeters"))
                        res<-dplyr::distinct(res, species, decimalLongitude, decimalLatitude, eventDate, countryCode, .keep_all=T)
                        return(res)
                      })
    file.remove(filepath)
    return(occs)
  } else {
    print(paste0("File ", filepath, " does not exist in working directory"))
    return(NULL)
  }
}
