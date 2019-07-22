#' Wrapper function to run a Generalized Dissimilarity Model
#'
#' @param occ.table data frame with species occurrence and sampling information. Each occurrence must contain XY coordinates species fields.
#' @param env.vars a rasterStack object of meaningful environmental variables to predict biological dissimilarity.
#' @param index the dissimilarity index to use. Valid options are "jaccard", "bray" or "betasim".
#' @param use.geo logical. If TRUE geographic distance is included in model as a predictor.
#' @param do.VarSel logical. If TRUE variable selection based on explained deviance is performed.
#' @param field.names a vector of field names in the occurrence table in the following order: species name, coordinates x, coordinates y.
#'
#' @details RunGDM adds functionality to the gdm function available in the gdm package. Specifically, 1) it adds
#' betasim, a richness independent measure of turnover (Lenon et al. 2001), 2) perfoms variable selection
#' by excluding from models variables with coefficients of 0 or that upon removal decrease the explained
#' deviance less than 5% and 3) computes principal components from the gdm transformed environmental variables
#' to aid in the visualization of turnover.
#'
#' @return A list with the following objects:
#' \describe{
#'   \item{occ.table}{The input occurrence table}
#'   \item{gdm.res}{A gdm object returned by the gdm function in gdm-package}
#'   \item{gdm.rasters}{A RasterStack of gdm transformed environmental predictors returned by the gdm.transform function in gdm-package}
#'   \item{gdm.map}{A list returned by the MapGDMLight function}
#' }
#' @references
#' Lennon, J. J., Koleff, P., Greenwood, J. J. D., & Gaston, K. J. (2001). The geographical structure of
#' British bird distributions: diversity, spatial turnover and scale. Journal of Animal Ecology, 70(6),
#' 966-979.
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
#' }
#' @export

RunGDM <- function(occ.table, env.vars, index, use.geo, do.VarSel, field.names){
  occ.table <- subset(occ.table, select=field.names[1:3])
  occ.table$cell <- raster::cellFromXY(env.vars, occ.table[,2:3])
  occ.table <- unique(na.omit(occ.table))

  if(index=="betasim"){
    #Create dissimilarity matrix using Beta-Sim distance measure (Kreft & Jetz 2010)
    sp.matrix <- reshape2::dcast(occ.table, as.formula(paste("cell ~",field.names[1])),
                       fun.aggregate = function(x){as.numeric(length(x)>0)})
    betasim.dist <- vegan::designdist(sp.matrix[, 2:ncol(sp.matrix)], method = "1-(a/(pmin(b,c)+a))", abcd = TRUE)
    betasim.dist <- data.frame(site = sp.matrix[, 1], as.matrix(betasim.dist))
    predData <- data.frame(site=betasim.dist$site, xFromCell(env.vars, betasim.dist$site),
                           yFromCell(env.vars, betasim.dist$site), env.vars[betasim.dist$site])
    colnames(predData)[2:3] <- field.names[2:3]
    gdm.table <- gdm::formatsitepair(betasim.dist, bioFormat=3, XColumn = field.names[2], YColumn = field.names[3],
                                sppColumn = field.names[1], siteColumn = "site",
                                predData = predData)
  } else {
    gdm.table <- gdm::formatsitepair(occ.table, bioFormat=2, XColumn = field.names[2], YColumn = field.names[3],
                                sppColumn = field.names[1], siteColumn = "cell",
                                predData=env.vars, dist=index)

  }
  gdm.table <- na.omit(gdm.table)
  gdm.rast <- gdm::gdm(gdm.table, geo=use.geo)

  if(do.VarSel){
    coef.df <- data.frame(name=rep(gdm.rast$predictors, each=3), coef=gdm.rast$coefficients)
    coef.df <- reshape2::dcast(coef.df, "name ~ .",value.var="coef", fun=sum)

    #Force geographic to remain even if not significant
    if(use.geo){
      coef.df <-coef.df[coef.df$name!="Geographic", ]
    }
    ind0 <- which(coef.df$. > 0)
    print(paste("removing variable", coef.df$name[-ind0], "with 0 coefficient"))
    vars2keep <- coef.df$name[ind0]
    test.names <- c(colnames(gdm.table)[1:6], paste0("s1.",vars2keep), paste0("s2.",vars2keep))
    ind.cols <- colnames(gdm.table)%in%test.names
    gdm.table <- gdm.table[, ind.cols]
    m1 <- gdm::gdm(gdm.table, geo=use.geo)

    #Remove iteratively variables. CSIRO method
    dev.dif <- 0
    while(dev.dif < 0.05){
      coef.df <- data.frame(name=rep(m1$predictors, each=3), coef=m1$coefficients)
      coef.df <- reshape2::dcast(coef.df, "name ~ .", value.var="coef", fun=sum)
      if(use.geo){
        coef.df <-coef.df[coef.df$name!="Geographic", ]
      }

      test.var <- which.min(coef.df$.)

      rem.vars <- paste0(c("s1.","s2."), coef.df$name[test.var])
      gdm.table.test <- gdm.table[, !(colnames(gdm.table)%in%rem.vars)]

      m2 <- gdm::gdm(gdm.table.test, geo=use.geo)

      if(is.null(m2$explained)|m2$explained==0){
        print("Model update has zero or null deviance. Undoing model update and exiting variable selection module")
        break()
      }

      dev.dif <- m1$explained - m2$explained
      if(dev.dif < 0.05){
        print(paste("Variable", coef.df$name[test.var],"not significant. Removing."))
        gdm.table <- gdm.table.test
        m1 <- m2
      } else {
        print("No more variables to remove from model")
      }
    }
    gdm.rast <- m1
  }

  env.vars.red <- env.vars[[which(names(env.vars)%in%gdm.rast$predictors)]]
  rastTrans <- gdm::gdm.transform(gdm.rast, env.vars.red)
  map.gbd.results <- MapGDMLight(rastTrans)
  return(list(occ.table=occ.table, gdm.res=gdm.rast, gdm.rasters=rastTrans, gdm.map=map.gbd.results))
}
