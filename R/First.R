##  spatstat.univar/R/First.R

.onLoad <- function(...) {
#  reset.spatstat.options()
}

.onAttach <- function(libname, pkgname) {
  vs <- read.dcf(file=system.file("DESCRIPTION", package="spatstat.univar"),
                 fields="Version")
  vs <- as.character(vs)
#  putSpatstatVariable("SpatstatUnivarVersion", vs)
  packageStartupMessage(paste("spatstat.univar", vs))
  return(invisible(NULL))
}

  
