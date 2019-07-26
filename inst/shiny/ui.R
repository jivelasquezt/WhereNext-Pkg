js <- "
$(document).on('shiny:value', function(evt){
if(evt.name == 'log'){
setTimeout(function(){
var objDiv = document.getElementById('log');
objDiv.scrollTop = objDiv.scrollHeight - objDiv.clientHeight;
}, 500); // wait 500 milliseconds
}
});
"

ui <- navbarPage(theme=shinythemes::shinytheme('sandstone'), id='tabs', collapsible=TRUE,
                 title="WhereNext",
                 tabPanel("Intro", value=0),
                 tabPanel("1 Occurrence Data", value=1),
                 tabPanel("2 Environmental Data", value=2),
                 tabPanel("3 Run GDM", value=3),
                 tabPanel("4 Recommend survey", value=4),
                 shinyjs::useShinyjs(),
                 fluidRow(column(3,
                                 conditionalPanel("input.tabs == 1",
                                                  h4(HTML("<b>1.1 Select study area<b>")),
                                                  radioButtons("aoiSrc", NULL, choices=list("By country"='aoi.ctr',"User shape"='usr.shp')),
                                                  conditionalPanel("input.aoiSrc=='aoi.ctr'",
                                                                   selectInput("country.sel", label="Select country:", choices = ISOcodes::ISO_3166_1$Name, width = '100%')
                                                  ),
                                                  conditionalPanel("input.aoiSrc=='usr.shp'",
                                                                   fileInput("aoi.shp", label="Enter full path to zipped shapefile:")
                                                  ),
                                                  HTML('<hr>'),
                                                  h4(HTML("<b>1.2 Select data source<b>")),
                                                  radioButtons("occSel", NULL,
                                                               choices = list("Download from GBIF" = 'db', "User-specified" = 'user'),
                                                               selected = 'db'),
                                                  conditionalPanel("input.occSel == 'db'",
                                                                   radioButtons("grp.type", label = "Select taxonomic rank:",
                                                                                choices = list("Class"='classKey', "Order"='orderKey',"Family"='familyKey'), inline=TRUE),
                                                                   textInput("grp.selection", label="Enter biological group (capitalize first letter):", width = '100%'),
                                                                   actionButton("login", "GBIF login"),
                                                                   actionButton(inputId ="download.gbif", label="Download Data"),
                                                                   downloadButton("downloadOccTable", "")
                                                  ),
                                                  conditionalPanel("input.occSel == 'user'",
                                                                   radioButtons("fileType", "Choose file type:",choices=list("Excel"="xls", "Plain text"="txt"), inline=T),
                                                                   conditionalPanel("input.fileType == 'txt'",
                                                                                    radioButtons("sep", "Field delimiter", choices=list("(,) comma" = ',', "(;) semicolon" = ';', "tab" = '\t'), inline=T),
                                                                                    radioButtons("dec", "Decimal delimiter", choices=list("(.) period" = '.',"(,) comma" = ','), inline=T)
                                                                   ),
                                                                   textInput("user.occs", "Enter path to occurrence file:"),
                                                                   actionButton("user.occs.go", "Go")
                                                  ),
                                                  HTML('<hr>'),
                                                  h4(HTML("<b>1.3 Clean occurrence data (optional)<b>")),
                                                  actionButton(inputId ="run.clean", label="Clean data"),
                                                  actionButton(inputId ="reset.occs", label="Reset clean"),
                                                  downloadButton("downloadCleanTable", "")
                                 ),
                                 conditionalPanel("input.tabs==2",
                                                  h4(HTML("<b>2.1 Select environmental data source<b>")),
                                                  radioButtons("envSource", NULL,
                                                               choices=list("WorldClim (cropped by selected country)"='wc', "User-specified"='user.env')),
                                                  conditionalPanel("input.envSource=='wc'",
                                                                   radioButtons("wc.res", "Select raster resolution:",
                                                                                choices=list("2.5'"=2.5, "5'"=5, "10'"=10), inline=TRUE),
                                                                   actionButton(inputId ="download.clim", label="Download Worldclim data"),
                                                                   downloadButton("downloadEnvVars", ""),
                                                                   helpText("", a(href="https://www.worldclim.org/bioclim", target="_blank", "WorldClim variable definitions"))
                                                  ),
                                                  conditionalPanel("input.envSource=='user.env'",
                                                                   fileInput(inputId ="env.files", label="Select input rasters:",
                                                                             multiple = TRUE, width="100%")
                                                  ),
                                                  HTML('<hr>'),
                                                  h4(HTML("<b>2.2 Reduce colinearity (optional)<b>")),
                                                  h5(HTML("<b>Correlation threshold:<b>")),
                                                  splitLayout(
                                                    numericInput(inputId="cor.threshold", label=NULL, value=0.8, min = 0, max = 1, width="100%"),
                                                    actionButton(inputId ="remove.corvars", "Run", width="100%"),
                                                    actionButton("env.reset","Reset", width="100%"),
                                                    cellWidths=c("40%","30%","30%")
                                                  ),
                                                  HTML('<hr>'),
                                                  h4(HTML("<b>Visualize environmental layer<b>")),
                                                  selectInput("selVar","Select layer:",choices = "")
                                 ),
                                 conditionalPanel("input.tabs == 3",
                                                  h4(HTML("<b>3.1 Estimate cell stats<b>")),
                                                  actionButton(inputId="est.rich", label="Run estimation"),
                                                  downloadButton("downloadCellStats", ""),
                                                  HTML('<hr>'),
                                                  h4(HTML("<b>3.2 Select sites with complete surveys<b>")),
                                                  sliderInput(inputId="n.surveys", label="Surveys", value=0, min=0, max=50),
                                                  sliderInput(inputId="richness", label="Richness", value=0, min=0, max=50),
                                                  sliderInput(inputId="completeness", label="Completeness", value=0, min=0, max=1, round=-1),

                                                  actionButton(inputId="sel.cells", label="Filter occurrences"),
                                                  downloadButton("downloadFilteredTable", ""),
                                                  HTML('<hr>'),
                                                  h4(HTML("<b>3.3 Configure GDM and run<b>")),
                                                  selectInput(inputId="gdm.beta", label="Disimilarity index",
                                                              choices=list("Jaccard"="jaccard","Bray-Curtis"="bray","Beta similarity"="betasim"),
                                                              selected = "jaccard"),
                                                  checkboxInput(inputId="gdm.dist",label="Include distance in model",value=TRUE),
                                                  checkboxInput(inputId="gdm.varsel",label="Do variable selection",value=TRUE),
                                                  actionButton(inputId="gdm.run",label="Run GDM"),
                                                  downloadButton("downloadGDM", "")
                                 ),
                                 conditionalPanel("input.tabs == 4",
                                                  h4(HTML("<b>4.1 Find complementary sites<b>")),
                                                  radioButtons("edSel","Find sampling suggestions",
                                                               choices = list("Anywhere in study area"='ed.all',"From preselected sites"='ed.some')),
                                                  conditionalPanel("input.edSel=='ed.some'",
                                                                   fileInput("ed.sites", "Upload preselected sites table (.csv)")
                                                  ),
                                                  actionButton("run.ed","Compute complementarity"),
                                                  HTML('<hr>'),
                                                  h4(HTML("<b>4.2 Iterative site selection<b>")),
                                                  conditionalPanel("input.edSel=='ed.all'",
                                                                   radioButtons("edAction", "Select action & press go to find the next site",choices=list("Add"='ed.add',"Reject"='ed.reject',"Modify"='ed.modify'), inline =T)),
                                                  conditionalPanel("input.edSel=='ed.some'",
                                                                   radioButtons("edAction", "Select action & press go to find the next site",choices=list("Add"='ed.add',"Reject"='ed.reject'), inline =T)),
                                                  conditionalPanel("input.edAction=='ed.modify'",
                                                                   splitLayout(
                                                                     numericInput("customCoords.y","Latitude", value=0),
                                                                     numericInput("customCoords.x","Longitude", value=0),
                                                                     cellWidths = c("50%", "50%")
                                                                   )
                                                  ),
                                                  actionButton("ed.go","Go"),
                                                  HTML('<hr>'),
                                                  h4(HTML("<b>Download results")),
                                                  downloadLink("downloadSuggestions","Selected sites table"),
                                                  br(),
                                                  downloadLink("downloadED","ED Complementarity raster")
                                 )
                 ),
                 column(9,
                        shinybusy::add_busy_spinner(spin = "fading-circle"),
                        conditionalPanel("input.tabs == 0",
                                         includeMarkdown("Rmd/introMD.md")
                        ),
                        conditionalPanel("input.tabs != 0",
                                         tags$head(
                                           tags$script(HTML(js))
                                         ),
                                         tags$style(HTML("#log {height:80px}")),
                                         br(),
                                         verbatimTextOutput("log")
                        ),
                        br(),
                        conditionalPanel("input.tabs != 0",
                                         tabsetPanel(id = "inTabset",
                                                     tabPanel('Map', leafletOutput("map", width = "100%", height = 500), type = 1),
                                                     tabPanel("Occurrence table", DT::dataTableOutput("occ.table")),
                                                     tabPanel('Results',
                                                              conditionalPanel("input.tabs == 1",
                                                                               h4(HTML("<b>Citation:<b>")),
                                                                               verbatimTextOutput("citation")),
                                                              conditionalPanel("input.tabs == 3", verbatimTextOutput("gdmSummary")),
                                                              conditionalPanel("input.tabs == 4", plotOutput("plot"))
                                            #                  conditionalPanel("input.tabs == 4", DT::dataTableOutput("ed.table"))
                                                     )
                                         )
                        )
                 )
                 )
)
