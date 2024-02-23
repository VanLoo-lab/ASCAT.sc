run_PloidyModifier <- function(){
  suppressPackageStartupMessages(require(shiny))
   
  appdir <- paste0(system.file(package = "ASCAT.sc"), "/PloidyModifier/")
  shiny::runApp(appdir)
  }
