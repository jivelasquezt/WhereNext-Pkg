#' Initial computation of complementarity and matrices for iterative site selection
#'
#' @param gdm.rast A gdm model object
#' @param occ.table A data frame with species occurrence information. Must contain a decimalLongitude and decimalLatitude field (x and y coordinates, respectively)ate field
#' @param env.vars A rasterStack object of gdm transformed environmental predictors. Obtained by using the gdm.transform option in the gdm package.
#' @param subset.size An integer specifying the number of demand points for continuous p-median evaluation. Default = 10000. See Faith & Walker (1996) and Manion (2009).
#' @param search.type Either a "1" or "2" string specifying whether to identify potential survey sites across the entire study area or from a set of candidate sites, respectively.
#' @param candidate.sites A data frame of candidate survey sites. Only necessary when using search.type="2". Must only have columns decimalLongitude and decimalLatitude, in that order.
#' @details This function implements the stage 1 algorithm proposed by Manion & Ridges (2009). Specifically, it computes the minimum
#' environmental distance from a set of demand points to survey sites (nnDistance in params object) and the
#' environmental distance from demand points to candidate sites, either all non-sampled cells in the input grids or
#' a user-defined set of candidate sites (m2 in params object). Then it finds the cell
#' (XY coordinates provided in selCoords) that reduces the most the environmental distance
#'  to demand points compared to the survey sites (nnDistanceUpdated in params object).
#'  The raster output measures each cells' relative contribution to improve survey
#'  coverage or complementarity if the cell is selected (Funk et al. 2005).
#'
#' @return A list with the following objects:
#' \describe{
#'   \item{out.raster}{A raster object displaying the complementarity of candidate sites}
#'   \item{initED}{A relative measure of total current survey coverage}
#'   \item{outED}{A relative measure of total survey coverage if the suggested site is surveyed}
#'   \item{selCoords}{Coordinates of the suggested survey site (i.e. the most complementary)}
#'   \item{params}{A list object to input in function FindNext}
#' }
#'
#' @note This is a very memory intensive function as it will try to fit a matrix of
#' roughly subset.size rows by ncell(env.vars) columns. PreFindNext will request that
#' you have at least 3 times the necessary memory available in your system to run.
#'
#' @references
#' Funk, V. A., Richardson, K. S., & Ferrier, S. (2005). Survey-gap analysis in expeditionary research: where
#' do we go from here?. Biological Journal of the Linnean Society, 85(4), 549-567.
#'
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
#' }
#' @export
#' @import raster


PreFindNext<-function(gdm.rast, occ.table, env.vars, subset.size=1000, search.type="1", candidate.sites){
  start.time <- Sys.time()
  #Create template and output raster
  t.raster <- !is.na(prod(env.vars)) #Cells that ARE NOT NA in any single layer become 1. NA cells become 0
  out.raster <- t.raster
  out.raster[out.raster==0] <- NA

  #Prepare table for gdm prediction in unsampled sites
  ind.na <- raster::extract(t.raster, occ.table[, c("decimalLongitude", "decimalLatitude")]) #identify whether any of the sampled cells have NA data
  occ.table <- occ.table[ind.na==1, ] #Extract only rows that are not NA in any environmental layer
  sampled <- unique(cellFromXY(t.raster, occ.table[, c("decimalLongitude", "decimalLatitude")]))
  usmpled <- Which(t.raster==1, cells=T)
  usmpled.ind <- which(!usmpled%in%sampled)
  usmpled <- usmpled[usmpled.ind] # Remove sampled cells from the unsampled cells set

  if(search.type=="2"){
    if(!is.null(candidate.sites)){
      cell2add <- na.omit(unique(cellFromXY(env.vars, candidate.sites)))
      usmpled.ind <- which(!usmpled%in%cell2add)
      usmpled <- usmpled[usmpled.ind] # Remove candidate sites from the unsampled cells set
      sampled.ind <- which(!sampled%in%cell2add)
      sampled <- sampled[sampled.ind] #Remove candidate sites from the sampled cells set
    } else {
      print("Must provide a candidate sites data frame to use search type 2")
      return()
    }
  }
  subset.size <- subset.size * ncell(out.raster)/length(usmpled) #multiply subset size for a factor to account for NA cells in raster
  sample.grid <- sampleRegular(out.raster, size=min(ncell(out.raster), subset.size), cells = T)
  sample.grid.nna <- apply(sample.grid, 1, function(x) !anyNA(x))
  demand.pts <- sample.grid[sample.grid.nna, "cell"]
  print(paste("Created",length(demand.pts),"regularly spaced demand points"))

  #Test if objects will fit in memory
#  lst.obj.sz <- object.size(matrix(rnorm(length(demand.pts)) * length(usmpled), nrow=length(demand.pts), ncol=length(usmpled))) * 3
#  test.mem <- lst.obj.sz > benchmarkme::get_ram()
#  if(test.mem){
#    print(paste("Cannot allocate extra", format(lst.obj.sz, standard="SI", units="GB"), "in memory"))
#    return(lst.obj.sz)
#  }

  #Calculate distance from demand points (rows) to sampled sites
  print("Computing distances from demand points to survey sites")
  site.grid <- expand.grid(demand.pts, sampled)
  gdm.pred <- DistGDM(site.grid, gdm.rast, env.vars)

  m1 <- data.frame(c.unsampl = site.grid[, 1], c.sampled = site.grid[, 2],  dis=gdm.pred)
  m1 <- tidyr::spread(m1, "c.sampled", "dis")

  rm(gdm.pred, site.grid)
  if(ncol(m1)<=2){
    nnDistance <- m1[, 2]
  } else {
    nnDistance <- apply(m1[, 2:ncol(m1)], 1, min)  #For all demand points, find the minimum distance to survey sites
  }

  print("Computing distances from demand points to candidate sites")
  if(search.type=="1"){
    #Calculate distance from demand points (rows) to all unsampled sites
      site.grid2 <- expand.grid(demand.pts, usmpled)
      gdm.pred2 <- DistGDM(site.grid2, gdm.rast, env.vars)
      m2 <- data.frame(c.unsampl = site.grid2[, 1], c.sampled = site.grid2[, 2],  dis=gdm.pred2)
      m2 <- tidyr::spread(m2, "c.sampled", "dis")


  } else if (search.type=="2"){
    site.grid2 <- expand.grid(demand.pts, cell2add)
    gdm.pred2 <- DistGDM(site.grid2, gdm.rast, env.vars)
    m2 <- data.frame(c.unsampl = site.grid2[, 1], c.sampled = site.grid2[, 2],  dis=gdm.pred2)
    m2 <- tidyr::spread(m2, "c.sampled", "dis")
  }
  rm(m1, gdm.pred2, site.grid2)
  gc()
  iniTotalED <- sum(nnDistance)   #Compute initial ED value

  #Compute ED complementarity of unsampled cells.
  m3 <- m2[, 2:ncol(m2)]
  imdv.mat <- matrix(rep(nnDistance, ncol(m3)), ncol=ncol(m3))
  ind.mat <- m3 > nnDistance
  m3 <- (ind.mat)*imdv.mat + m3*(!ind.mat)
  ed.comp <- apply(m3, 2, sum)
  ed.comp.dif <- iniTotalED - ed.comp
  #Nearest neighbor distance of demand points to survey site + site that minimizes total ED complementarity
  nnDistanceUpdated <- m3[, which.min(ed.comp)]
  outTotalED <- sum(nnDistanceUpdated)

  #Create ED complementarity raster
  out.raster[as.numeric(names(ed.comp.dif))] <- ed.comp.dif

  if(search.type==2){ #For search type 2, remove values for cells not sampled and not candidate
    out.raster.tmp <- out.raster
    out.raster.tmp[] <- NA
    out.raster.tmp[c(cell2add, sampled)] <- out.raster[c(cell2add, sampled)]
    out.raster <- out.raster.tmp
    plot(out.raster.tmp)
    rm(out.raster.tmp)
  }
  selCoords <- xyFromCell(out.raster, raster::which.max(out.raster)) #Coordinates of site that is most complementary

  print(difftime(Sys.time(), start.time, units = "mins"))
  return(list(out.raster=out.raster, initED=iniTotalED, outED=outTotalED, selCoords=selCoords, params=list(nnDistance=nnDistance, nnDistanceUpdated=nnDistanceUpdated, m2=m2)))
}

