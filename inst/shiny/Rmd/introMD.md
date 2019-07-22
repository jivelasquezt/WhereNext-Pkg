## User guide

Welcome to *WhereNext*, a recommendation system to identify iteratively sites that complement the most existing biodiversity surveys. *WhereNext* will guide you through the process of developing a Generalized Dissimilarity Model (GDM) for your group and area of interest and identifying iteratively survey sites that maximize your potential to record new species.

To get started, browse in order through tabs and follow the numbered steps from 1 to 8.

1. Select study area: either select a country from the dropdown list or browse to a zipped shapefile in your browser (e.g. C:/Shapes/country.zip). 

2. Select data source: either download occurrence data from GBIF or upload your own occurrence file. To download from gbif, you must create first an account at https://www.gbif.org/ and provide your user credentials using the login button. Next, you can download all occurrence data for a class, order or family within GBIF. Notice that if you input a lower taxonomic rank name and select a taxonomic rank above it, *WhereNext* will download the data for the highest taxonomic rank selected. For example, if you select the taxonomic rank “class” but type “Hylidae” (a family), WhereNext will download all the data for class Amphibia, within the study area.

    __The citation of your GBIF download will appear in the results tab__

    To upload your own occurrences, the file may be either in Excel or plain text format. In the latter case, you must know the field and delimiters of your file. Some files may failed to parse due to quotes or other escape characters within the file. In those cases it may be best to import the text file to Excel and then into *WhereNext*. The occurrence file may have any number of columns, but it must include at least the fields `gbifID`, `species`, `decimalLongitude`, `decimalLatitude`, `eventDate`, `countryCode`, `locality` and `recordedBy`. Bear in mind:
    * Fields `gbifID`, `locality` and `recordedBy` are not used in computations and may be filled with dummy data (e.g. "") if not available. 
    * `eventDate` must parse to a valid date format by package *lubridate* (e.g. YYYY-MM-DD).
    * Field `countryCode` must be the ISO2 country code. Not using the correct ISO2 code will result in removing all your records if you use the clean data option.
    * An `individualCount` column will be interpreted as abundance and may be useful to run GDM using the Bray-Curtis distance. In case this column is omitted it will be filled with 1 by default. 
    * The fields `decimalLongitude` and `decimalLatitude may be used to store x and y coordinates, regardless of the projection.
    * Bear in mind that *WhereNext* will assume that the coordinates of your occurrence data match the projection of the study area shapefile. If they don’t match, you will have to make any necessary transformations outside *WhereNext*.

    Optional: coordinates may be cleaned by pulsing the clean data button. This action will run package *CoordinateCleaner* (Zizka et al. 2019) which will run the `countries`, `capitals`, `centroids`, `equal`, `gbif`, `institutions`, `outliers`, `seas` and `zeros` tests. To learn more about these tests check the *CoordinateCleaner* package documentation (https://cran.r-project.org/web/packages/CoordinateCleaner/index.html).

3. Select environmental data: either download WorldClim data or provide your own environmental raster layers. Bear in mind that the extent and resolution of your raster layers must be the same and that they must be in the same projection as the study area shapefile. 

    Optional: you may specify a correlation threshold and run the reduce collinearity tool to keep only variables with correlations among them below the specified threshold.

4. Estimate cell stats: in *WhereNext*, occurrences are aggregated spatially into cells, according to the environmental data resolution. For each cell, the number of surveys, species richness and sampling completeness is calculated. Sampling completeness varies between 0 (incomplete) and 1 (complete), and is the ratio of observed species richness over expected species richness, estimated by the jackniffe2 non-parametric estimator implemented in package *vegan* (Oksanen et al. 2019).  By default, the number of surveys in each cell is the number of dates in which at least an occurrence was reported. You may change this behavior by manipulating the `eventDate` field of your occurrence file outside *WhereNext* (however, it must parse to a valid date format).

5. Select sites with complete surveys: GDM relies on proper estimation of the dissimilarity between sampled cells. To do so, only cells that have been sufficiently sampled should be considered in GDM development. *WhereNext* allows you to place lower bounds on the cells to include in the analysis, based on the number of surveys, its richness and/or sampling completeness. Using the sliders, explore the effect of selecting different lower bounds based on the knowledge of your study system and data availability. When satisfied, pulse the "select sites" button.

6. Configure and run GDM: enter your choice of dissimilarity index, geographic distance and variable selection and pulse "Run GDM". *WhereNext* will configure your input files to run a GDM using the *gdm* package (Manion et al. 2018). The map tab will show the resulting GDM (to see, desactivate other raster layers). Different colors in this layer represent different biotic communities.

7. Find complementary sites: you have two choices to find complementary survey sites, either to select from new survey areas anywhere in the study area or to select survey areas from a list of sites. The latter option is appropriate to prioritize sampling when there’s a predefined list of sites that can be visited. This list must be a .csv file with only  two columns, `decimalLongitude` and `decimalLatitude`. After clicking on the “Compute ED” button, a layer will appear on the map depicting the most complementary areas to existing survey sites and the first sampling suggestion. 

8. Select sites iteratively: after computing the first complementarity surface, you may either add the site to a list of sites to survey, reject the suggestion or modify the suggestion by entering alternative coordinates. If using the modify option, clicking on the map will populate the latitude and longitude fields. The selected sites to survey are shown in the map and results tab.

### References
Glenn Manion, Matthew Lisk, Simon Ferrier, Diego Nieto-Lugilde, Karel Mokany and Matthew C. Fitzpatrick (2018). gdm: Generalized Dissimilarity Modeling. R package version 1.3.11. https://CRAN.R-project.org/package=gdm

Jari Oksanen, F. Guillaume Blanchet, Michael Friendly, Roeland Kindt, Pierre Legendre, Dan McGlinn, Peter R. Minchin, R. B. O'Hara, Gavin L. Simpson, Peter Solymos, M. Henry H. Stevens, Eduard Szoecs and Helene Wagner (2019). vegan: Community Ecology Package. R package version 2.5-5. https://CRAN.R-project.org/package=vegan

Zizka A, Silvestro D, Andermann T, Azevedo J, Duarte Ritter C, Edler D, Farooq H, Herdean A, Ariza M, Scharn R, Svanteson S, Wengstrom N, Zizka V, Antonelli A (2019). “CoordinateCleaner: standardized cleaning of occurrence records from biological collection databases.” _Methods in Ecology and Evolution_, 0. doi: 10.1111/2041-210X.13152 (URL:http://doi.org/10.1111/2041-210X.13152), R package version 2.0-11, <URL: https://github.com/ropensci/CoordinateCleaner>

### Issues
Post any issues to https://github.com/jivelasquezt/WhereNext-Pkg/issues

### Author
Jorge Velásquez-Tibatá (jorge.velasquez at tnc.org).

