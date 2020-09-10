#' Launch ChromSCape
#'
#' @export
#'
#' @examples launchApp()
#' 
#' @return Launches the shiny application
#' @import shiny
#'  
launchApp <- function(){
  
  # shinyApp(ui = shinyAppUI, server = shinyAppServer)
  source(file.path(system.file(package="ChromSCape"),
                   "Module_preprocessing_filtering_and_reduction.R"),)
  shiny::runApp(system.file(package="ChromSCape"),launch.browser = TRUE)
}
