#' Estimate richness and sampling completeness
#'
#' @description This code will estimate richness using non-parametric estimators for sampled cells, defined by using the raster provided in the arguments. Completeness will be estimated based on the proportion of estimated species that are actually recorded.
#'
#' @param occ.table data frame with species occurrence and sampling information. Must contain coordinates x, coordinates y, species, a sampling event id (e.g. date) and count field.
#' @param in.raster a raster in the same projection of coordinates provided in the occurrence table and with the desired cell size for analyses.
#' @param field.names a vector of field names in the occurrence table in the following order: coordinates x, coordinates y, event id, species name, count and cell field.
#'
#' @return A data frame of observed species richness(Species), number of surveys (n), estimated species
#' richness (chao, jack1, jack2, boot) and completeness of sampling (C_chao, C_jack1, C_jack2, C_boot) for
#' each sampled cell. XY coordinates of each cell and standard errors (suffix .se) are also provided.
#'
#' @examples
#' \dontrun{
#' #To compute cell stats for colombian bats
#' library(lubridate)
#' library(raster)
#' library(rgbif)
#' library(WhereNext)
#'
#' #Get occurrence data
#' gbif.key <- name_backbone(name = "Chiroptera")
#' gbif.res <- DownloadGBIF(gbif.key$orderKey, "your username", "your email", "your password", "CO") #Enter your GBIF credentials here
#'
#' #Get environmental data
#' col <- getData("GADM", country="COL",level=0)
#' env.vars <- getData("worldclim", var="bio", res=5)
#' env.vars <- crop(env.vars, col)
#' env.vars <- mask(env.vars, col)
#' env.vars <- Normalize(env.vars) #Normalize environmental variables
#' env.vars <- RemCorrLayers(env.vars, 0.8) #Remove variables correlated more than r=0.8.
#'
#' #Do minimal occ.table cleaning
#' occ.table <- gbif.res$occ.table
#' occ.table$eventDate <- as_date(occ.table$eventDate)
#' occ.table$individualCount <- 1 #Data is presence only
#' occ.table.clean <- subset(occ.table, !is.na(eventDate) & taxonRank=="SPECIES")
#' row.names(occ.table.clean) <- 1:nrow(occ.table.clean)
#' occ.table.clean <- CoordinateCleaner::clean_coordinates(occ.table.clean,
#'                                                         lon="decimalLongitude",
#'                                                         lat="decimalLatitude",
#'                                                         species="species",
#'                                                         countries = "countryCode",
#'                                                         value="clean",
#'                                                         tests=c("capitals","centroids", "equal", "gbif",
#'                                                                 "institutions", "outliers", "seas","zeros"))
#'
#' #Estimate cell sampling stats & filter occurrence data
#' occ.table.clean$cell <- cellFromXY(env.vars, occ.table.clean[, c("decimalLongitude","decimalLatitude")])
#' cell.stats <- RichSamp(occ.table.clean, env.vars, c("decimalLongitude","decimalLatitude","eventDate","species","individualCount","cell"))
#' selected.cells <- cell.stats$cell[which(cell.stats$n>3&cell.stats$Species>5)] #Consider places with at least 3 sampling events and 5 species recorded as well sampled
#' occ.table.sel <- occ.table.clean[which(occ.table.clean$cell %in% selected.cells), ] #Use only ocurrences of well sampled cells
#' }
#' @export


RichSamp <- function(occ.table, in.raster, field.names){
  #Reformat matrix
  M <- reshape2::dcast(occ.table, as.formula(paste(field.names[3], "+ cell", "~", field.names[4])),
             fun.aggregate=function(x) as.numeric(length(x)>0), value.var=field.names[5])
  M <- na.omit(M)

  #Find cells sufficiently sampled
  singNdoub<-function(siteBySp, cell){
    sub.m<-siteBySp[siteBySp$cell==cell, ]
    nrecs<-apply(sub.m[, 3:ncol(sub.m)], 2, sum) #number of records
    nobs<-sum(nrecs>0) #Species observed
    smps<-nrow(sub.m) #Samples
    sing<-sum(nrecs==1) #Singletons
    doub<-sum(nrecs==2) #Doubletons
    suf<-(sing/nobs<0.5)*1 #Criteria to define a place as sufficiently sampled according to Anne Chao in EstimateS User Guide
    return(list(nobs=nobs, smps=smps, sing=sing, doub=doub, suf=suf))
  }

  results <- sapply(unique(M$cell), singNdoub, siteBySp=M)
  results <- t(results)
  results <- as.data.frame(results)
  results$cell <- unique(M$cell)
  results[,"nobs"] <- unlist(results[,"nobs"])
  results[,"smps"] <- unlist(results[,"smps"])
  results[,"sing"] <- unlist(results[,"sing"])
  results[,"doub"] <- unlist(results[,"doub"])
  results[,"suf"] <- unlist(results[,"suf"])
  results[,"cell"] <- unlist(results[,"cell"])

  selected.cells<-results$cell[results$suf==1]

  #Estimate richness
  Msub <- M[M$cell%in%selected.cells, ] # Exclude from analysis cells not sufficiently sampled according to Chao criterion
  result <- vegan::specpool(Msub[,3:ncol(M)], Msub[,2])
  result$cell <- as.numeric(row.names(result))
  xyR <- raster::xyFromCell(in.raster, result$cell)
  result <- cbind(result, xyR)

  #Estimate sampling completeness
  result$C_chao <- result$Species/result$chao
  result$C_jack1 <- result$Species/result$jack1
  result$C_jack2 <- result$Species/result$jack2
  result$C_boot <- result$Species/result$boot

  return(result)
}
