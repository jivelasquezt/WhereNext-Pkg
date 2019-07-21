[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
# WhereNext (1.0.0)
WhereNext is a recommendation system to identify iteratively sites that complement the most existing biodiversity surveys. WhereNext will guide you through the process of developing a Generalized Dissimilarity Model (GDM) for your group and area of interest and identifying iteratively survey sites that maximize your potential to record new species.

Install and run the WhereNext shiny app using the following R code.

```R
install.packages("devtools")
devtools::install_github("jivelasquezt/WhereNext-Pkg")
library(WhereNext)
RunWhereNext()
```