#'Finds the next most complementary or disimilar site on average to sampled sites within a set of candidate sites.
#'
#' @param pre.params The object resulting of running PreFindNext
#' @param action One of the following
#' \describe{
#'   \item{add}{Add previously suggested site to the list of sampled sites and find next suggestion}
#'   \item{reject}{Reject previously suggested site and find next best site}
#'   \item{modify}{Modify previous suggestion with a new site with coordinates specified in custom.coords}
#'   }
#' @param custom.coords A vector with the XY coordinates of a site to add to the list of selected sites
#'
#' @details This function implements stage 2 of the algorithm outlined by Manion & Ridges (2009).
#' Specifically, it will take the results from PreFindNext or a previous iteration of
#' FindNext and for option "add" will update the distance (nnDistance in params) from demand points to candidate
#' sites considering that the sites suggested in the previous iteration is accepted (m2 in params). Then it
#' finds the cell that reduces the most the environmental distance to demand points (XY coordinates in selCoords)
#' compared to the previous iteration (nnDistanceUpdated in params object). The "reject" option
#' will exclude the suggestion found in the previous iteration and find the next best site. The
#' "modify" option will select the cell corresponding to custom.coords and find the next best site.
#'
#'
#' @return A list with the following objects:
#' \describe{
#'   \item{out.raster}{A raster object displaying the complementarity of candidate sites}
#'   \item{initED}{The total complementarity of existing sites}
#'   \item{outED}{The total complementarity after selecting the suggested site}
#'   \item{selCoords}{Coordinates of the suggested survey site (i.e. the most complementary)}
#'   \item{params}{A list object to input in function FindNext}
#' }
#' #' @references
#' Manion, G., & Ridges, M. (2009). An optimisation of the survey gap analysis technique to
#' minimise computational complexity and memory resources in order to accommodate fine grain environmental and
#' site data. In 18th World IMACS/MODSIM Congress, Cairns, Australia (pp. 13-17).
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
#'
#' #Run and map GDM
#' m1 <- RunGDM(occ.table.sel, env.vars, "bray", TRUE, TRUE, c("species", "decimalLongitude", "decimalLatitude"))
#' plotRGB(m1$gdm.map$pcaRast)
#'
#' #Find initial survey suggestions
#' init <- PreFindNext(m1$gdm.res, m1$occ.table, m1$gdm.rasters, subset.size=1000, search.type="1")
#' plot(init$out.raster)
#' points(init$selCoords)
#'
#' #Find additional survey suggestions
#' update <- FindNext(init, "add")
#' plot(update$out.raster)
#' points(update$selCoords)
#' }
#' @export


FindNext <- function(pre.params, action, custom.coords=NULL){
  start.time <- Sys.time()
  m2 <- pre.params$params$m2
  if(ncol(m2) < 4){
    return("No sites left to prioritize")
  }
  if(action == "add"){
    ind.cell<-cellFromXY(pre.params$out.raster, pre.params$selCoords)
    m2[, which(colnames(m2) == ind.cell)] <- NULL #Remove added cell
    m3 <- m2[, 2:ncol(m2)]
    imdv.mat <- matrix(rep(pre.params$params$nnDistanceUpdated, ncol(m3)), ncol = ncol(m3))
    ind.mat <- m3 > pre.params$params$nnDistanceUpdated
    m3 <- (ind.mat)*imdv.mat + m3*(!ind.mat)
    ed.comp <- apply(m3, 2, sum)
    iniTotalED <- pre.params$outED
    ed.comp.dif <- iniTotalED - ed.comp
    nnDistance <- pre.params$params$nnDistanceUpdated
    out.raster <- pre.params$out.raster
  }

  if(action=="reject"){
    ind.cell <- cellFromXY(pre.params$out.raster, pre.params$selCoords)
    out.raster <- pre.params$out.raster
    out.raster[ind.cell] <- NA
    m2[, as.character(ind.cell)] <- NULL #Remove rejected cell
    m3 <- m2[, 2:ncol(m2)]
    imdv.mat <- matrix(rep(pre.params$params$nnDistance, ncol(m3)), ncol = ncol(m3))
    ind.mat <- m3 > pre.params$params$nnDistance
    m3 <- (ind.mat)*imdv.mat + m3*(!ind.mat)
    ed.comp <- apply(m3, 2, sum)
    iniTotalED <- pre.params$initED
    ed.comp.dif <- iniTotalED - ed.comp
    nnDistance <- pre.params$params$nnDistance
  }

  if(action=="modify"){
    ind.cell <- cellFromXY(pre.params$out.raster, custom.coords)
    pre.params$params$nnDistanceUpdated <- m2[, as.character(ind.cell)]
    m3 <- m2[, 2:ncol(m2)]
    imdv.mat <- matrix(rep(pre.params$params$nnDistanceUpdated, ncol(m3)), ncol = ncol(m3))
    ind.mat <- m3 > pre.params$params$nnDistanceUpdated
    m3 <- (ind.mat)*imdv.mat + m3*(!ind.mat)
    ed.comp <- apply(m3, 2, sum)
    iniTotalED <- sum(pre.params$params$nnDistanceUpdated)
    ed.comp.dif <- iniTotalED - ed.comp
    nnDistance <- pre.params$params$nnDistanceUpdated
    out.raster <- pre.params$out.raster
  }
  #Nearest neighbor distance of demand points to survey site + site that minimizes total ED complementarity
  nnDistanceUpdated <- m3[, which.min(ed.comp)]
  outTotalED <- sum(nnDistanceUpdated)

  #Create ED complementarity raster
  out.raster[as.numeric(names(ed.comp.dif))] <- ed.comp.dif
  selCoords <- raster::xyFromCell(out.raster, as.numeric(names(ed.comp)[which.min(ed.comp)])) #Coordinates of site that is most complementary

  print(difftime(Sys.time(), start.time, units = "mins"))
  return(list(out.raster = out.raster, initED = iniTotalED, outED = outTotalED, selCoords=selCoords, params = list(nnDistance = nnDistance, nnDistanceUpdated = nnDistanceUpdated, m2 = m2)))
}
