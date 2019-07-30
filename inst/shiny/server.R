options(shiny.maxRequestSize=5000*1024^2)
logInit <- function() {
  intro <- '***Welcome to WhereNext***'
  brk <- paste(rep('------', 14), collapse='')
  expl <- 'Please find messages for the user in this log window.'
  return(c(paste(intro, brk, expl, brk, "\n",sep='\n')))
}

tableStyle <- "<style>
table {
font-family: arial, sans-serif;
border-collapse: collapse;
width: 100%;
}

th {
background-color: #333333;
color:white;
}

td, th {
border: 1px solid #dddddd;
text-align: left;
padding: 2px;
}
</style>"

server <- function(input, output, session) {
  rv <- reactiveValues(logs = logInit())
  output$log <- renderText({rv$logs})
  ###########################
  # Module 1: Enter GBIF credentials
  ###########################
  observeEvent(input$login, {
    showModal(modalDialog(
      title = "Enter GBIF credentials",
      textInput("user",label="Username", value="",placeholder = "sample",width="40%"),
      passwordInput("password", label="Password", value="", placeholder = "GBIF password",width="40%"),
      textInput("user.email",label="Email", value="", placeholder = "sample@email.com"),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("ok", "OK")
      )
    ))
  })

  observeEvent(input$ok,
               removeModal()
  )

  ###########################
  # Module 1: GBIF download
  ###########################
  observeEvent(input$download.gbif,{
    print("Running GBIF download module")
    rv$isGBIF <- TRUE
    if(is.null(input$ok)){
      showModal(modalDialog("Please provide your GBIF credentials login info download data\n.
                            To create a free GBIF account go to: https://www.gbif.org"))
      return()
    }
    tryCatch(gbif.key <- rgbif::name_backbone(name = input$grp.selection),
             error = function(e){
               rv$logs <- paste(rv$logs, e, "\n")
             })
    if(!exists("gbif.key")){
      rv$logs <- paste(rv$logs, "An error occurred getting the taxon key. Check your internet connection and try again.\n")
      return()
    }
    if(gbif.key$matchType == "EXACT"){
      key <- gbif.key[input$grp.type]
      if(key=="NULL"){
        rv$logs <- paste(rv$logs, input$grp.selection, "name exists but does not match selected taxonomic hierarchy\n")
        rv$logs <- paste(rv$logs, "Check taxonomic hierarchy and try again\n")
        return()
      }
    } else {
      rv$logs<-paste(rv$logs, input$grp.selection, "not found. Please check spelling and try again\n")
      return()
    }

    if(input$aoiSrc=="aoi.ctr"){
      ind.aoi <- which(ISOcodes::ISO_3166_1$Name == input$country.sel)
      tryCatch(rv$aoi <- raster::getData(name="GADM",country=ISOcodes::ISO_3166_1$Alpha_3[ind.aoi], level=0),
               error = function(e) {
                 rv$logs <- paste(rv$logs, print(e),"\n")
               })
      if(is.null(rv$aoi)){
        rv$logs <- paste(rv$logs, "An error ocurred downloading country shapefile.\nCheck your internet connection, search parameters and/or GBIF credentials and try again\n")
        return()
      }
      tryCatch(gbif.res <- DownloadGBIF(key, input$user, input$user.email, input$password, ISOcodes::ISO_3166_1$Alpha_2[ind.aoi]),
               error = function(e){
                 rv$logs <- paste(rv$logs, print(e),"\n")
               })
      if(!exists("gbif.res")){
        rv$logs <- paste(rv$logs, "An error ocurred downloading data from GBIF. \nCheck your internet connection, search parameters and/or GBIF credentials and try again\n")
        return()
      }

    } else if(input$aoiSrc=="usr.shp") {
      index <- sapply(gregexpr("\\/", input$aoi.shp$datapath), tail, 1)
      shp.dir <- substr(input$aoi.shp$datapath, 1, index)
      zip::unzip(input$aoi.shp$datapath, exdir=shp.dir, junkpaths =T)
      shp.name <- list.files(shp.dir,"*.shp$")[1]
      shp.name <- substr(shp.name, 0, nchar(shp.name)-4)
      test.shp <- anyNA(file.exists(paste0(shp.dir,shp.name,c(".shp",".dbf",".shx"))))
      if(test.shp){
        rv$logs <- paste(rv$logs,"ERROR: Zipped shapefile must have .shp, .dbf, and .shx files.\n
                         Fix your files and try again")
        return()
      }
      rv$aoi <- rgdal::readOGR(shp.dir, shp.name)
      tryCatch(gbif.res <- DownloadGBIF(key, input$user, input$user.email, input$password, NULL, rv$aoi),
               error = function(e){
                 rv$logs <- paste(rv$logs, print(e),"\n")
               })
      if(!exists("gbif.res")){
        rv$logs <- paste(rv$logs, "An error ocurred downloading data from GBIF. \nCheck your internet connection, search parameters and/or GBIF credentials and try again\n")
        return()
      }
    }

    output$citation <- gbif.res$citation
    occ.table <- gbif.res$occ.table
    if(nrow(occ.table)==0){
      rv$logs <- paste(rv$logs, "No records found for selected group in selected area. Change query parameters and try again\n")
      return()
    }
    occ.table <- unique(subset(occ.table, taxonRank=="SPECIES"&(coordinateUncertaintyInMeters < 5000| is.na(coordinateUncertaintyInMeters)),
                               select=c("gbifID", "species", "decimalLongitude" ,"decimalLatitude", "eventDate","countryCode","locality","recordedBy")))
    occ.table$individualCount <- 1
    occ.table$eventDate <- lubridate::as_date(occ.table$eventDate)
    occ.table <- occ.table[!is.na(occ.table$eventDate), ]
    rv$logs <-paste(rv$logs, nrow(occ.table), "records remain after removing records without valid date\n")
    
    row.names(occ.table) <- 1:nrow(occ.table)
    occ.table$countryCode <- ISOcodes::ISO_3166_1$Alpha_3[match(occ.table$countryCode, ISOcodes::ISO_3166_1$Alpha_2)] #Change country code to ISO3 for verification in CoordinateCleaner

    #Refine records based on study area shapefile
    ind.coords <- sf::st_as_sf(data.frame(id=1:nrow(occ.table), decimalLongitude=occ.table$decimalLongitude, decimalLatitude=occ.table$decimalLatitude),
                           coords=c("decimalLongitude","decimalLatitude"))

    aoi.sf <- as(rv$aoi,"sf")
    if(is.na(sf::st_crs(ind.coords))){
      sf::st_crs(ind.coords) <- sf::st_crs(aoi.sf)
    }
    ind.coords <- as.data.frame(sf::st_join(ind.coords, aoi.sf, join = sf::st_intersects))
    rec.sel <- which(!is.na(ind.coords[, 2]))

    occ.table <- occ.table[rec.sel, ]
    rv$sp.data <- occ.table
    rv$sp.data.orig <- occ.table
    rv$logs <- paste(rv$logs,"After excluding records outside study area, ", nrow(occ.table), " records remain\n")

    #Map records
    n.max <- min(nrow(rv$sp.data), 1e5)
    disp.order <- sample(1:nrow(rv$sp.data), n.max)
    labs <- lapply(disp.order, function(i) {
      paste0( tableStyle,
              '<table style="width:100%">',
              '<tr><th>Attribute</th><th>Value</th></tr>',
              '<tr><td> ID','</td><td><a href="https://www.gbif.org/occurrence/', rv$sp.data[i, "gbifID"],'" target="_blank">See on GBIF</a></td></tr>',
              '<tr><td> Species','</td><td>', rv$sp.data[i, "species"],'</td></tr>',
              '<tr><td> Date','</td><td>', rv$sp.data[i, "eventDate"],'</td></tr>',
              '</table>')
    })
    map %>% clearGroup(group="Occ data") %>%
      clearControls() %>%
      addLayersControl(baseGroups=c("Basemap"), overlayGroups = c("Occ data"),
                       options = layersControlOptions(collapsed = FALSE, autoZIndex = TRUE)) %>%
      fitBounds(lng1=min(rv$sp.data$decimalLongitude, na.rm=T),
                lng2=max(rv$sp.data$decimalLongitude, na.rm=T),
                lat1=min(rv$sp.data$decimalLatitude, na.rm=T),
                lat2=max(rv$sp.data$decimalLatitude, na.rm=T)) %>%
      addCircleMarkers(lng = ~decimalLongitude, lat = ~decimalLatitude, data = rv$sp.data[disp.order, ], group = "Occ data",
                       fillColor = 'dodgerblue', fillOpacity = 0.6, weight = 2, radius = 5,
                       popup = lapply(labs, htmltools::HTML),
                       clusterOptions = markerClusterOptions())
    output$downloadOccTable <- downloadHandler("occData.csv",
                                               content = function(file){
                                                 write.csv(rv$sp.data, file, row.names=FALSE)
                                               })
    })

  m <- leaflet(options=leafletOptions(preferCanvas=TRUE)) %>% setView(0, 0, zoom = 2) %>%
    addTiles(group = 'BaseMap') %>%
    addLayersControl(baseGroups=c("Basemap"))

  output$map <- renderLeaflet(m)
  map <- leafletProxy("map")


  options <- list(autoWidth = TRUE,
                  columnDefs = list(list(width = '200px', targets = c(1,6,7)),
                                    list(className = 'dt-center', targets = 0:8),
                                    list(
                                      targets = c(1,6,7),
                                      render = JS(
                                        "function(data, type, row, meta) {",
                                        "return type === 'display' && data.length > 25 ?",
                                        "'<span title=\"' + data + '\">' + data.substr(0, 24) + '...</span>' : data;",
                                        "}"))),
                  scrollX=TRUE, scrollY=400)
  output$occ.table <- DT::renderDataTable(rv$sp.data, filter = 'top', rownames = FALSE, options = options)

  ##############################
  #Module 1B: User file upload
  ##############################

  observeEvent(input$user.occs.go,{
    print("Running user occurrence upload module")
    #Read text file
    rv$isGBIF <- FALSE
    if(input$fileType == "txt"){
      tryCatch(n.recs <- bigreadr::nlines(input$user.occs),
               error=function(e){
                 rv$logs <- paste(rv$logs, "Error:", e, "\n")

               })
      if(!exists("n.recs")){
        rv$logs <- paste(rv$logs, "An error occurred reading file:", input$user.occs, ".\n Check file format and use forward slash as path separator.\n")
        return()
      }
      if(n.recs > 1e6){
        rv$logs <- paste(rv$logs, "Your file has", n.recs,"records. WhereNext will attempt to split and read the file and individualCount will be ignored\n")
        tryCatch(user.file <- bigreadr::big_fread1(input$user.occs, 
                                every_nlines=1000000, 
                                .transform = function(x) {
                                  res<-dplyr::select(x, c("gbifID", "species", "decimalLongitude" ,"decimalLatitude", "eventDate","countryCode","locality","recordedBy"))
                                  res<-dplyr::distinct(res,species, decimalLongitude ,decimalLatitude, eventDate,countryCode, .keep_all=T)
                                  return(res)
                                }),
                 error=function(e){
                   rv$logs <- paste(rv$logs, "Error:", e, "\n")
                 })
      } else {
        tryCatch(user.file <- read.delim(input$user.occs, sep=input$sep, dec=input$dec, header=TRUE,
                                         stringsAsFactors=F, row.names = NULL, quote="", fill=FALSE, as.is=TRUE),
                 error=function(e){
                   rv$logs<-paste(rv$logs,"Error:", e, "\n")
                   return()
                   })
      }
    } else if(input$fileType == "xls"){
      tryCatch(user.file <- openxlsx::read.xlsx(input$user.occs, sheet = 1, colNames = TRUE),
               error=function(e){
                 rv$logs<-paste(rv$logs,"Error:", e, "\n")
                 return()}
      )
    }
    if(!exists("user.file")){
      rv$logs <- paste(rv$logs,"An error occurred reading file. Check file format.\n")
      return()
    }
    
    #Checking columns
    if(is.null(user.file$individualCount)){
      rv$logs <- paste(rv$logs,"individualCount column not found. Assuming data is presence-only\n")
      user.file$individualCount <- 1
    } else {
      user.file$individualCount[is.na(user.file$individualCount)] <- 1
    }
    required.cols <- c("gbifID","species", "decimalLongitude" ,"decimalLatitude", "eventDate","countryCode","locality","recordedBy","individualCount")
    check.cols <- required.cols%in%colnames(user.file)
    if(all(check.cols)){
      user.file <- unique(subset(user.file, select=required.cols))
      user.file$eventDate <- lubridate::as_date(user.file$eventDate)
      if(length(which(is.na(user.file$eventDate))) == nrow(user.file)){
        rv$logs <-paste(rv$logs, "All dates failed to parse. Check that your dates are in ymd or ymd HM format and try again\n")
        return()
      } else {
        user.file <- user.file[!is.na(user.file$eventDate), ]
        rv$logs <-paste(rv$logs, nrow(rv$sp.data), "records remain after removing records without date\n")
      }
    } else {
      notfound <- required.cols[which(!check.cols)]
      rv$logs <- paste(rv$logs, "Column(s):", paste(notfound, collapse = ", "), "not found. File not imported.\n")
      return()
    }
    rv$logs <- paste(rv$logs, "Imported occurrence file with", nrow(user.file), "records.\n")
    
    #Remove rows with empty species cell
    user.file <- user.file[!is.na(user.file$species) & user.file$species != "", ]
    rv$logs <-paste(rv$logs, nrow(user.file), "records remain after removing records with empty species cell\n")

    #Get study area shapefile
    if(input$aoiSrc=="aoi.ctr"){
      ind.aoi <- which(ISOcodes::ISO_3166_1$Name == input$country.sel)
      tryCatch(rv$aoi <- raster::getData(name="GADM",country=ISOcodes::ISO_3166_1$Alpha_3[ind.aoi], level=0),
               error = function(e) {
               rv$logs <- paste(rv$logs, print(e),"\n")
      })
      if(is.null(rv$aoi)){
        rv$logs <- paste(rv$logs, "An error ocurred downloading country shapefile. \nCheck your internet connection, search parameters and/or GBIF credentials and try again\n")
        return()
      }
    } else if(input$aoiSrc=="usr.shp") {
      index <- sapply(gregexpr("\\/", input$aoi.shp$datapath), tail, 1)
      shp.dir <- substr(input$aoi.shp$datapath, 1, index)
      zip::unzip(input$aoi.shp$datapath, exdir=shp.dir, junkpaths =T)
      shp.name <- list.files(shp.dir,"*.shp$")[1]
      shp.name <- substr(shp.name, 0, nchar(shp.name)-4)
      test.shp <- anyNA(file.exists(paste0(shp.dir, shp.name, c(".shp",".dbf",".shx"))))
      if(test.shp){
        rv$logs <- paste(rv$logs,"ERROR: Zipped shapefile must have .shp, .dbf, and .shx files.\n
                         Fix your files and try again")
        return()
      }
      rv$aoi <- rgdal::readOGR(shp.dir, shp.name)
    }

    #Refine records by study area shapefile
    ind.coords <- sf::st_as_sf(data.frame(id=1:nrow(user.file), decimalLongitude=user.file$decimalLongitude, decimalLatitude=user.file$decimalLatitude),
                           coords=c("decimalLongitude","decimalLatitude"))

    aoi.sf <- as(rv$aoi, "sf")
    if(is.na(sf::st_crs(ind.coords))){
      sf::st_crs(ind.coords) <- sf::st_crs(aoi.sf)
    }
    ind.coords <- as.data.frame(sf::st_join(ind.coords, aoi.sf, join = sf::st_intersects))
    rec.sel <-which(!is.na(ind.coords[, 2]))

    if(length(rec.sel)==0){
      rv$logs <- paste(rv$logs, "No records within selected study area\n")
      return()
    } else {
      rv$logs <- paste(rv$logs, "After excluding records outside study area,", length(rec.sel), "records remain.\n")
    }
    user.file <- user.file[rec.sel, ]
    user.file$countryCode <- ISOcodes::ISO_3166_1$Alpha_3[match(user.file$countryCode, ISOcodes::ISO_3166_1$Alpha_2)] #Change country code to ISO3 for verification in CoordinateCleaner
    row.names(user.file) <- 1:nrow(user.file)
    rv$sp.data.orig <- user.file
    rv$sp.data <- user.file
    
    #Map records
    n.max <- min(nrow(rv$sp.data), 1e5)
    disp.order <- sample(1:nrow(rv$sp.data),n.max)
    labs <- lapply(disp.order, function(i) {
      paste0( tableStyle,
              '<table style="width:100%">',
              '<tr><th>Attribute</th><th>Value</th></tr>',
              '<tr><td> ID','</td><td>', rv$sp.data[i, "gbifID"],'</td></tr>',
              '<tr><td> Species','</td><td>', rv$sp.data[i, "species"],'</td></tr>',
              '<tr><td> Date','</td><td>', rv$sp.data[i, "eventDate"],'</td></tr>',
              '</table>')
    })
    map %>% clearGroup(group="Occ data") %>%
      clearControls() %>%
      addLayersControl(baseGroups=c("Basemap"), overlayGroups = c("Occ data"),
                       options = layersControlOptions(collapsed = FALSE, autoZIndex = TRUE)) %>%
      fitBounds(lng1=min(rv$sp.data$decimalLongitude, na.rm=T),
                lng2=max(rv$sp.data$decimalLongitude, na.rm=T),
                lat1=min(rv$sp.data$decimalLatitude, na.rm=T),
                lat2=max(rv$sp.data$decimalLatitude, na.rm=T)) %>%
      addCircleMarkers(lng = ~decimalLongitude, lat = ~decimalLatitude, data = rv$sp.data[disp.order, ], group = "Occ data",
                       fillColor = 'dodgerblue', fillOpacity = 0.6, weight = 2, radius = 5,
                       popup = lapply(labs, htmltools::HTML),
                       clusterOptions = markerClusterOptions())
    output$citation <- renderPrint({"User uploaded data"})
    })

  ###########################
  # Module 1C: Clean data
  ###########################

  observeEvent(input$run.clean,{
    print("Running data cleaning module")

    if(is.null(rv$sp.data)|nrow(rv$sp.data)==0){
      rv$logs <- paste(rv$logs,"Occurrence data is null. Download/upload occurrences before proceeding\n")
      return()
    }

    if(nrow(rv$sp.data) == 0){
      rv$logs <- paste(rv$logs, "No occurrences left\n")
      return()
    }

    #Run coordinate cleaner
    row.names(rv$sp.data) <- 1:nrow(rv$sp.data)
    try(sp.data.clean <- CoordinateCleaner::clean_coordinates(rv$sp.data,
                                             lon="decimalLongitude",
                                             lat="decimalLatitude",
                                             species="species",
                                             countries = "countryCode",
                                             value="clean",
                                             tests=c("countries","capitals","centroids", "equal", "gbif",
                                                     "institutions", "outliers", "seas","zeros")))
    if(exists("sp.data.clean")){
      rv$sp.data <- sp.data.clean
      rv$logs <-paste(rv$logs, nrow(rv$sp.data), "records remain after running CoordinateCleaner\n")
    } else {
      rv$logs <-paste(rv$logs, "CoordinateCleaner failed. Trying now without country test\n")
      tryCatch({sp.data.clean <- CoordinateCleaner::clean_coordinates(rv$sp.data,
                                                                      lon="decimalLongitude",
                                                                      lat="decimalLatitude",
                                                                      species="species",
                                                                      countries = "countryCode",
                                                                      value="clean",
                                                                      tests=c("capitals","centroids", "equal", "gbif","institutions", "seas","zeros"))
               rv$sp.data <- sp.data.clean
               rv$logs <-paste(rv$logs, nrow(rv$sp.data), "records remain after running CoordinateCleaner\n")},
               error = function(e) {
                 rv$logs <-paste(rv$logs, e)
                 rv$logs <-paste(rv$logs, "CoordinateCleaner failed. No data cleaning performed.\n")
               })
    }
             
    #Map cleaned records
    n.max <- min(nrow(rv$sp.data), 1e5)
    disp.order <- sample(1:nrow(rv$sp.data), n.max)
    if(rv$isGBIF){
      labs <- lapply(disp.order, function(i) {
        paste0( tableStyle,
                '<table style="width:100%">',
                '<tr><th>Attribute</th><th>Value</th></tr>',
                '<tr><td> ID','</td><td><a href="https://www.gbif.org/occurrence/', rv$sp.data[i, "gbifID"],'" target="_blank">See on GBIF</a></td></tr>',
                '<tr><td> Species','</td><td>', rv$sp.data[i, "species"],'</td></tr>',
                '<tr><td> Date','</td><td>', rv$sp.data[i, "eventDate"],'</td></tr>',
                '</table>')
      })
    } else {
      labs <- lapply(disp.order, function(i) {
        paste0( tableStyle,
                '<table style="width:100%">',
                '<tr><th>Attribute</th><th>Value</th></tr>',
                '<tr><td> ID','</td><td>', rv$sp.data[i, "id"],'</td></tr>',
                '<tr><td> Species','</td><td>', rv$sp.data[i, "species"],'</td></tr>',
                '<tr><td> Date','</td><td>', rv$sp.data[i, "eventDate"],'</td></tr>',
                '</table>')
      })
    }

    map %>% clearGroup(group="Occ data") %>%
      clearControls() %>%
      addLayersControl(baseGroups=c("Basemap"), overlayGroups = c("Occ data"),
                       options = layersControlOptions(collapsed = FALSE, autoZIndex = TRUE)) %>%
      addCircleMarkers(lng = ~decimalLongitude, lat = ~decimalLatitude, data = rv$sp.data[disp.order, ], group = "Occ data",
                       fillColor = 'dodgerblue', fillOpacity = 0.6, weight = 2, radius = 5,
                       popup = lapply(labs, htmltools::HTML),
                       clusterOptions = markerClusterOptions()) %>%
      fitBounds(lng1=min(rv$sp.data$decimalLongitude,na.rm=T),
                lng2=max(rv$sp.data$decimalLongitude,na.rm=T),
                lat1=min(rv$sp.data$decimalLatitude,na.rm=T),
                lat2=max(rv$sp.data$decimalLatitude,na.rm=T))
    output$downloadCleanTable <- downloadHandler("occCleanData.csv",
                                                 content = function(file){
                                                   write.csv(rv$sp.data, file, row.names=FALSE)
                                                 })
  })
  #############################
  #Reset occurrences
  #############################
  observeEvent(input$reset.occs,{
    if(is.null(rv$sp.data.orig)){
      rv$logs <- paste(rv$logs, "Nothing to reset\n")
    } else {
      rv$sp.data <- rv$sp.data.orig
      rv$logs <- paste(rv$logs, "Occurrence records restored. Total records:", nrow(rv$sp.data), "\n")
      #Map records
      labs <- lapply(seq(nrow(rv$sp.data)), function(i) {
        paste0( tableStyle,
                '<table style="width:100%">',
                '<tr><th>Attribute</th><th>Value</th></tr>',
                '<tr><td> Species','</td><td>', rv$sp.data[i, "species"],'</td></tr>',
                '<tr><td> Date','</td><td>', rv$sp.data[i, "eventDate"],'</td></tr>',
                '</table>')
      })
      #Map records
      n.max <- min(nrow(rv$sp.data), 1e5)
      disp.order <- sample(1:nrow(rv$sp.data),n.max)
      if(rv$isGBIF){
        labs <- lapply(disp.order, function(i) {
          paste0( tableStyle,
                  '<table style="width:100%">',
                  '<tr><th>Attribute</th><th>Value</th></tr>',
                  '<tr><td> ID','</td><td><a href="https://www.gbif.org/occurrence/', rv$sp.data[i, "gbifID"],'" target="_blank">See on GBIF</a></td></tr>',
                  '<tr><td> Species','</td><td>', rv$sp.data[i, "species"],'</td></tr>',
                  '<tr><td> Date','</td><td>', rv$sp.data[i, "eventDate"],'</td></tr>',
                  '</table>')
        })
      } else {
        labs <- lapply(disp.order, function(i) {
          paste0( tableStyle,
                  '<table style="width:100%">',
                  '<tr><th>Attribute</th><th>Value</th></tr>',
                  '<tr><td> ID','</td><td>', rv$sp.data[i, "gbifID"],'</td></tr>',
                  '<tr><td> Species','</td><td>', rv$sp.data[i, "species"],'</td></tr>',
                  '<tr><td> Date','</td><td>', rv$sp.data[i, "eventDate"],'</td></tr>',
                  '</table>')
        })
      }
      map %>% clearGroup(group="Occ data") %>%
        clearControls() %>%
        addLayersControl(baseGroups=c("Basemap"), overlayGroups = c("Occ data"),
                         options = layersControlOptions(collapsed = FALSE, autoZIndex = TRUE)) %>%
        addCircleMarkers(lng = ~decimalLongitude, lat = ~decimalLatitude, data = rv$sp.data[disp.order, ], group = "Occ data",
                         fillColor = 'dodgerblue', fillOpacity = 0.6, weight = 2, radius = 5,
                         popup = lapply(labs, htmltools::HTML),
                         clusterOptions = markerClusterOptions())
    }
  })

  #################################
  #Module 2A: Download environmental data
  #################################
  observeEvent(input$download.clim, {
    print("Running climate download module")

    if(is.null(rv$sp.data)){
      rv$logs <- paste(rv$logs, "Error: first select occurrence data\n")
      return()
    }
    tryCatch(clim <- raster::getData("worldclim", var="bio", res=input$wc.res),
             error= function(e){
               rv$logs <- paste(rv$logs, e, "\n")
             })
    if(!exists("clim")){
      rv$logs <- paste(rv$logs, "Worldclim download failed. Check your internet connection and try again.\n")
      return()
    }
    rv$logs <- paste(rv$logs, "Worldclim data downloaded\n")
    clim.aoi <- raster::crop(clim, rv$aoi)
    clim.aoi <- raster::mask(clim.aoi, rv$aoi)
    clim.aoi <- Normalize(clim.aoi)
    rv$env.vars <- clim.aoi
    rv$env.vars.orig <- clim.aoi
    rv$logs <- paste(rv$logs, "Environmental data cropped, masked and normalized\n")
    cells.wdata <- length(raster::Which(!is.na(rv$env.vars[[1]]),cells=T))
    rv$logs <- paste(rv$logs, "Environmental data has", cells.wdata, "cells with data.\n")
    output$downloadEnvVars <- downloadHandler("envVars.tif",
                                              content = function(file){
                                                writeRaster(rv$env.vars, file, format="GTiff")
                                              })
  })

  #################################
  #Module 2B: Upload environmental data
  #################################
  observeEvent(input$env.files, {
    print("Running environmental data upload module")
    tryCatch(rv$env.vars.orig <- raster::stack(input$env.files$datapath),
             error=function(e){
               rv$logs <- paste(rv$logs, "Error:", e,"\n")
             })
    if(is.null(rv$env.vars.orig)){
      rv$logs <- paste(rv$logs, "An error occurred uploading rasters.\nRasters must be readable by function raster and should have matching extent and resolution.\n")
      return()
    }
    names(rv$env.vars.orig) <- input$env.files$name
    if(!is.null(rv$env.vars.orig)){
      rv$logs <- paste(rv$logs, "Variables: ", paste(names(rv$env.vars.orig),collapse=", "), "uploaded\n")
    }
    rv$env.vars.orig <- raster::crop(rv$env.vars.orig, rv$aoi)
    rv$env.vars.orig <- raster::mask(rv$env.vars.orig, rv$aoi)
    rv$env.vars <- rv$env.vars.orig
    rv$logs <- paste(rv$logs, "Environmental data cropped and masked\n")
    cells.wdata <- length(raster::Which(!is.na(rv$env.vars[[1]]), cells=T))
    rv$logs <- paste(rv$logs, "Environmental data has", cells.wdata, "cells with data.\n")

  })

  #################################
  #Module 2B: Remove colinearity
  #################################

  observeEvent(input$remove.corvars,{
    print("Running variable elimination module")

    if(is.null(rv$env.vars)){
      rv$logs <- paste(rv$logs, "Environmental variables not found\n")
      return()
    }
    if(input$cor.threshold > 1 | input$cor.threshold < 0){
      rv$logs <- paste(rv$logs, "Select correlation value within 0 to 1. Variables not removed\n")
      return()
    }
    rv$env.vars <- RemCorrLayers(rv$env.vars, input$cor.threshold)
    rv$logs <- paste(rv$logs, "After removing correlated variables:", paste(names(rv$env.vars), collapse=","), "remain\n")
  })

  #################################
  #Module 2C: Reset variables (previous to colinearity removal)
  #################################

  observeEvent(input$env.reset,{
    if(is.null(rv$env.vars)){
      rv$logs <- paste(rv$logs, "Environmental variables not found. Nothing to restore.\n")
      return()
    }
    rv$env.vars <- rv$env.vars.orig
    rv$logs <- paste(rv$logs, "Restoring variables:", paste(names(rv$env.vars), collapse=","), "\n")
  })

  #################################
  #Module 2D: Visualize variables
  #################################
  observe({
    if(!is.null(rv$env.vars)){
      updateSelectInput(session, "selVar", choices = names(rv$env.vars))
    }
  })

  observeEvent(input$selVar, {
    if(!is.null(rv$env.vars)){
      map %>%
        clearGroup("Env data") %>%
        clearControls() %>%
        addRasterImage(rv$env.vars[[input$selVar]], opacity =0.8, group = "Env data") %>%
        addLayersControl(baseGroups=c("Basemap"), overlayGroups = c("Occ data","Env data"),
                         options = layersControlOptions(collapsed = FALSE, autoZIndex = TRUE)) %>%
        fitBounds(lng1=raster::extent(rv$env.vars)@xmin,
                  lng2=raster::extent(rv$env.vars)@xmax,
                  lat1=raster::extent(rv$env.vars)@ymin,
                  lat2=raster::extent(rv$env.vars)@ymax)
    }
  })

  #################################
  #Module 3A: Estimate species richness
  #################################

  observeEvent(input$est.rich,{
    print("Running richness estimation module")

    rv$sp.data$cell <-raster::cellFromXY(rv$env.vars, rv$sp.data[,c("decimalLongitude","decimalLatitude")])
    rv$cell.richness <- RichSamp(rv$sp.data, rv$env.vars, c("decimalLongitude" ,"decimalLatitude", "eventDate", "species","individualCount"))
    rv$logs <- paste(rv$logs, "Estimated cell richness, completeness and survey effort per grid cell\n")
    output$cell.richness <- DT::renderDataTable(rv$cell.richness, options = options, rownames= FALSE)
    output$downloadCellStats <- downloadHandler("cellStatsMatrix.csv", content = function(file){
      write.csv(rv$cell.richness, file, row.names=FALSE)})
  })

  #################################
  #Module 3A: Select cells according to sampling criteria
  #################################
  observe({
    if(!is.null(rv$cell.richness)){
      updateSliderInput(session, "completeness", value=input$completeness, min=round(min(rv$cell.richness$C_chao),1), max=1, step=0.05)
      updateSliderInput(session, "richness", value=input$richness, min=min(rv$cell.richness$Species), max=50, step=1)
      updateSliderInput(session, "n.surveys", value=input$n.surveys, min=min(rv$cell.richness$n),max=50, step=1)
    }
  })

  observeEvent({input$richness
      input$completeness
      input$n.surveys
    }, {
    if(!is.null(rv$cell.richness)){
      survey.sites <- rv$cell.richness[which(rv$cell.richness$C_chao >= input$completeness &
                                               rv$cell.richness$Species >= input$richness &
                                               rv$cell.richness$n >= input$n.surveys),]
      rv$selected.cells <- survey.sites$cell
      labs <- lapply(seq(nrow(survey.sites)), function(i) {
        paste0( tableStyle,
                '<table style="width:100%">',
                '<tr><th>Attribute</th><th>Value</th></tr>',
                '<tr><td> Cell','</td><td>', survey.sites[i, "cell"],'</td></tr>',
                '<tr><td> Richness(R)','</td><td>', survey.sites[i, "Species"],'</td></tr>',
                '<tr><td> Surveys','</td><td>', survey.sites[i, "n"],'</td></tr>',
                '<tr><td> Estimated R(chao)','</td><td>', round(survey.sites[i, "chao"], 1),'</td></tr>',
                '</table>')
      })
      map %>%
        clearControls() %>%
        clearGroup("Selected sites") %>%
        addLayersControl(baseGroups=c("Basemap"), overlayGroups = c("Occ data","Env data","Selected sites"),
                         options = layersControlOptions(collapsed = FALSE, autoZIndex = TRUE)) %>%
        addCircleMarkers(lng = ~x, lat = ~y, data = survey.sites, group = "Selected sites",
                         fillColor = 'red', fillOpacity = 0.6, weight = 2, radius = 5,
                         label = lapply(labs, htmltools::HTML))
    }
  })

  observeEvent(input$sel.cells,{
    print("Running cell selection module")

    if(is.null(rv$cell.richness)){
      rv$logs <- paste(rv$logs, "Must first estimate richness, completeness and sampling effort\n")
    } else {
      rv$selected.occs <- rv$sp.data[which(rv$sp.data$cell%in%rv$selected.cells), ]
      rv$logs <- paste(rv$logs, "Selected", nrow(rv$selected.occs), "occurrences from", length(rv$selected.cells),"cells.\n")
      output$downloadFilteredTable <- downloadHandler("selectedOccData.csv",
                                                      content = function(file){
                                                        write.csv(rv$selected.occs, file, row.names = FALSE)
                                                      })
    }
  })


  #################################
  #Module 3C: Run GDM
  #################################

  observeEvent(input$gdm.run, {
    print("Running GDM module")
    if(is.null(rv$selected.occs)){
      rv$logs <- paste(rv$logs, "Must filter occurrences first.\n")
    } else {
      rv$gdm <- RunGDM(rv$selected.occs, rv$env.vars, input$gdm.beta,
                       input$gdm.dist, input$gdm.varsel, c("species", "decimalLongitude", "decimalLatitude"))
      wmcrs <- "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext +no_defs"

      rv$gdm$gdm.map$pcaRast <- projectRaster(rv$gdm$gdm.map$pcaRast, projectExtent(rv$gdm$gdm.map$pcaRast, wmcrs))

      map %>%
        clearControls() %>%
        addRasterImage(rv$gdm$gdm.map$pcaRast[[1]], colors=WhereNext:::rgbPalette(rv$gdm$gdm.map$pcaRast), opacity =0.8, group = "GDM") %>%
        addLayersControl(baseGroups=c("Basemap"), overlayGroups = c("Occ data","Env data","Selected sites","GDM"),
                         options = layersControlOptions(collapsed = FALSE, autoZIndex = TRUE))
      rv$logs <- paste(rv$logs, "Computed GDM from input occurrence and environmental data\n")
    }
    output$downloadGDM <- downloadHandler("gdm.tif", content = function(file){
      writeRaster(rv$gdm$gdm.map$pcaRast, file, format="GTiff")})
    output$gdmSummary <- renderPrint({
      summary(rv$gdm$gdm.res)
    })
  })


  #################################
  #Module 4A: Identify survey priorities from entire study area
  #################################
  observeEvent(input$run.ed, {
    print("Running ED complementarity module")

    if(is.null(rv$gdm)){
      rv$logs <- paste(rv$logs, "Must run GDM first\n")
      return()
    } else {
      if(input$edSel=="ed.all"){
        rv$ed.res <- PreFindNext(rv$gdm$gdm.res,
                                 rv$gdm$occ.table,
                                 rv$gdm$gdm.rasters,
                                 1000,
                                 "1")
        if(class(rv$ed.res)=="object_size"){
          rv$logs <- paste(rv$logs,"Cannot allocate extra", format(rv$ed.res, standard="SI", units="GB"), "in memory.\n Use larger cell size or reduce extent and try again.\n")
        }
      } else {
        tryCatch(candidate.sites <- read.csv(input$ed.sites$datapath),
                 error=function(e){
                   rv$logs <- paste(rv$logs, "Error reading csv file:", e, "\n")
                   return()})
        test.cols <- identical(colnames(candidate.sites), c("decimalLongitude", "decimalLatitude"))
        if(!test.cols){
          rv$logs <- paste(rv$logs, "Candidate site file must be comma separated (.csv) and contain only decimalLongitude and decimalLatitude columns\n")
        } else {
          rv$logs <- paste(rv$logs, "Loaded candidate site file\n")
        }

        rv$ed.res <- PreFindNext(rv$gdm$gdm.res,
                                 rv$gdm$occ.table,
                                 rv$gdm$gdm.rasters,
                                 1000,
                                 "2",
                                 candidate.sites)
#        if(class(rv$ed.res)=="object_size"){
#          rv$logs <- paste(rv$logs,"Cannot allocate extra", format(rv$ed.res, standard="SI", units="GB"), "in memory.\n Use larger cell size or reduce extent and try again.\n")
#        }
      }
      if(nrow(rv$ed.res$selCoords)>1){
        rv$logs <- paste(rv$logs, "More than one cell with the highest complementarity. Selecting first cell.\n")
      }
      rv$ed.table <- data.frame(x=rv$ed.res$selCoords[1, 1], y= rv$ed.res$selCoords[1, 2],initED=rv$ed.res$initED[1], outED=rv$ed.res$outED[1])
      pal <- colorQuantile(c("#ffeda0","#feb24c","#f03b20"), n=5, values(rv$ed.res$out.raster),
                          na.color = "#00000000", alpha=TRUE)
      map %>%
        clearControls() %>%
        addRasterImage(rv$ed.res$out.raster, colors=pal, opacity=0.8, group = "Complementarity") %>%
        addLegend(pal = pal, values = values(rv$ed.res$out.raster), title = "Complementarity", position = "bottomright",
                  labFormat = labelFormat(digits=1), opacity = 0.8, group = "Legend") %>%
        addLayersControl(overlayGroups = c("Occ data","Env data","Selected sites","GDM","ED Complementarity","Suggested"), baseGroups=c("Basemap"), 
                         options = layersControlOptions(collapsed = FALSE, autoZIndex = TRUE)) %>%
        addCircleMarkers(lng = ~x, lat = ~y, data = rv$ed.table, group = "Suggested",
                         fillColor = 'cyan', fillOpacity = 0.6, weight = 2, radius = 5,
                         label = paste("Site: ", 1:nrow(rv$ed.table)))
      rv$logs <- paste(rv$logs, "Computed ED Complementarity\n")
      ed.diff <- rv$ed.table$initED[1] - rv$ed.table$outED[1]
      output$plot <- renderPlot({plot(0:nrow(rv$ed.table),c(rv$ed.table$initED[1],rv$ed.table$outED), type="l",
           xlab="Sites selected", ylab="Total complementarity", lwd=2, col="blue",
           main = paste("Iteration difference =", round(ed.diff, 2)))})
      rv$stop.wn <- FALSE
    }
  })

  #################################
  #Module 4B: Iterative selection of sites
  #################################
  observeEvent(input$ed.go,{
    if(rv$stop.wn){
      rv$logs <- paste(rv$logs, "No more cells to select. WhereNext won't suggest any further sites on current parameters.\n")
      return()
    }
    print("Running iterative site selection module")
    if(is.null(rv$ed.res)){
      rv$logs <- paste(rv$logs, "Must Run ED first\n")
      return()
    }
    if(is.character(rv$ed.res)){
      rv$logs <- paste(rv$logs, "No gain from further iterations of WhereNext\n")
      return()
    }
    if(input$edAction=="ed.add"){
      rv$ed.res <- FindNext(rv$ed.res, "add")
      rv$ed.res$selCoords
      rv$logs <- paste(rv$logs, "Added previous suggestion to list and found next site\n")
    }
    if(input$edAction=="ed.reject"){
      if(nrow(rv$ed.table)==0){
        rv$logs <- paste(rv$logs, "No suggestions to reject\n")
        return()
      }
      rv$logs <- paste(rv$logs, "Added previous suggestion to list and found next site\n")
      rv$ed.res <- FindNext(rv$ed.res, "reject")
      rv$ed.table <- rv$ed.table[-nrow(rv$ed.table), ]
      rv$logs <- paste(rv$logs, "Rejected previous suggestion from list and found next site\n")
    }
    if(input$edAction=="ed.modify"){
      if(nrow(rv$ed.table)==0){
        rv$logs <- paste(rv$logs, "No suggestions to modify\n")
        return()
      }
      customCoords <- data.frame(x=input$customCoords.x, y=input$customCoords.y)
      check.na <- raster::extract(rv$ed.res$out.raster, customCoords)
      if(is.na(check.na)){
        rv$logs <- paste(rv$logs, "No environmental data for provided coordinates\n")
        return()
      }

      rv$ed.res <- FindNext(rv$ed.res, "modify", customCoords)
      rv$ed.table[nrow(rv$ed.table), 1:2] <- customCoords
      rv$ed.table$outED[nrow(rv$ed.table)] <- rv$ed.res$initED[1]
      rv$logs <- paste(rv$logs, "Modified last suggestion from list and found next site\n")
    }

    #Update ed.table
    if(is.character(rv$ed.res)){
      rv$logs <- paste(rv$logs, "No gain from further iterations of WhereNext\n")
      return()
    }
    if(nrow(rv$ed.res$selCoords)>1){
      rv$logs <- paste(rv$logs, "More than 1 cell has the highest complementarity. Selecting first cell.\n")
      rv$stop.wn <- TRUE
      return()
    }
    rv$ed.table <- rbind(rv$ed.table,
                         data.frame(x=rv$ed.res$selCoords[1,1], y=rv$ed.res$selCoords[1,2], initED=rv$ed.res$initED[1], outED=rv$ed.res$outED[1]))
    row.names(rv$ed.table)<-1:nrow(rv$ed.table)
    pal <- colorQuantile(c("#ffeda0","#feb24c","#f03b20"), n=5, values(rv$ed.res$out.raster),
                         na.color = "#00000000", alpha=TRUE)
    map %>%
      clearControls() %>%
      clearGroup(group =c("Complementarity", "Suggested")) %>%
      addRasterImage(rv$ed.res$out.raster, colors=pal, opacity=0.8, group = "Complementarity") %>%
      addLegend(pal = pal, values = values(rv$ed.res$out.raster), title = "Complementarity", position = "bottomright",
                labFormat = labelFormat(digits=1), opacity = 0.8, group = "Legend") %>%
      addLayersControl(overlayGroups = c("Occ data","Env data","Selected sites","GDM", "ED Complementarity","Suggested"),baseGroups=c("Basemap"), 
                       options = layersControlOptions(collapsed = FALSE, autoZIndex = TRUE)) %>%
      addCircleMarkers(lng = ~x, lat = ~y, data = rv$ed.table, group = "Suggested",
                       fillColor = 'cyan', fillOpacity = 0.6, weight = 2, radius = 5,
                       label = paste("Site: ", 1:nrow(rv$ed.table)))

    output$downloadSuggestions <- downloadHandler("surveySuggestions.csv",
                                                  content = function(file){
                                                    write.csv(rv$ed.table, file, row.names = FALSE)
                                                  })
    ed.diff <- rv$ed.table$outED[nrow(rv$ed.table)-1] - rv$ed.table$outED[nrow(rv$ed.table)]
    output$plot <- renderPlot({plot(0:nrow(rv$ed.table),c(rv$ed.table$initED[1],rv$ed.table$outED), type="l",
         xlab="Sites selected", ylab="Total complementarity", lwd=2, col="blue", 
         main = paste("Iteration difference =", round(ed.diff,2)))})
  })

  observeEvent(input$map_click, {
    click <- input$map_click
    text<-paste0("Lat: ", round(click$lat, 2), "- Lon: ", round(click$lng, 2))

    map %>%
      clearPopups() %>%
      addPopups(click$lng, click$lat, text)

    updateTextInput(session, "customCoords.y", value = round(click$lat,8))
    updateTextInput(session, "customCoords.x", value = round(click$lng,8))
  })


  #Output survey suggestions table
  output$ed.table <- DT::renderDataTable(rv$ed.table,
                                         options = list(autoWidth = TRUE, columnDefs = list(list(width = '40%', targets = 4)),
                                                        scrollX=FALSE, scrollY=400), rownames= TRUE)
  output$downloadED <- downloadHandler("edComplementarity.tif", content = function(file){
    writeRaster(rv$ed.res$out.raster, file, format="GTiff")})

  }
