#' Get GBIF data
#' @description This is an internal function to retrieve compressed (.zip) downloaded GBIF data in SIMPLE_CSV based on a occ_download object from the rgbif package. It's similar to rgbif:::occ_download_import but allows random sampling of rows and is more stable (so far).

#' @param data.request an occ_download object returned by the rgbif occ_download function
#' @param n.max integer. The maximum number of records before triggering the subset option.
#' @param subset.size integer. The number of records to subset randomly from file
#' @return a data.frame of downloaded gbif data

GetGBIFData <- function(data.request, n.max = 2e6, subset.size = 2e6){
  key <- rgbif::occ_download_meta(data.request)$key
  filepath <- paste0(key, ".zip")
  if(file.exists(filepath)){
    unzip(filepath)
    n.recs <- bigreadr::nlines(paste0(key,".csv"))
    print(paste("occurrence file has", n.recs, "records"))
    if(n.recs > n.max){
      warning("subsetting to a maximum of 2 million records")
      occs <- SampleCSV(paste0(key, ".csv"), subset.size)
    } else {
      occs <- read.delim(paste0(key, ".csv"), sep="\t", quote="", fill=FALSE)
    }
    file.remove(filepath)
    return(occs)
  } else {
    print(paste0("File ", filepath, " does not exist in working directory"))
    return(NULL)
  }
}
