#' @title Run WhereNext
#'
#' @description This function starts the shiny application WhereNext
#'
#' @return NULL
#'
#' @author Jorge Velásquez-Tibatá
#'
#' @note This implementation of WhereNext is limited to occurrence files of up to 2 million rows.
#' If greater, it will sample at random up to 2 million rows in unix based systems or none in Windows.
#' Only up to 100.000 records selected at random will be shown in the map tab.
#'
#' @examples WhereNext()
#'
#' @export

RunWhereNext <- function ()
{
  app_path <- system.file("shiny", package = "WhereNext")
  return(shiny::runApp(app_path, launch.browser = TRUE))
}
