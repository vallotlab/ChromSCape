#' Launch ChromSCape
#'
#' Main function to launch ChromSCape in your favorite browser. You can pass
#' additional parameters that you would pass to shiny::runApp
#' (\code{\link[shiny]{runApp}})
#'
#' @param launch.browser Wether to launch browser or not
#' @param ... Additional parameters passed to \code{\link[shiny]{runApp}}
#'
#' @return Launches the shiny application
#' @import shiny
#'
#' @export
#'
#' @examples
#' \dontrun{
#' launchApp()
#' }
launchApp <- function(launch.browser=TRUE, ...){
    source(file.path(system.file(package="ChromSCape"),
                "Module_preprocessing_filtering_and_reduction.R"),)
    shiny::runApp(system.file(package="ChromSCape"),launch.browser = TRUE, ...)
}
