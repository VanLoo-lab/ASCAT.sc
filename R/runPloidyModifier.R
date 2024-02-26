run_PloidyModifier <- function(){
  suppressPackageStartupMessages(require(shiny))
   
  appdir <- paste0(system.file(package = "ASCAT.sc"), "/inst/PloidyModifier/")
  shiny:::runApp(appdir)
  }
