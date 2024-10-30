#' utility functions
#'
#' Copyright (c) 2024 Adrian Baddeley, Tilman Davies and Martin Hazelton

check.bandwidth <- function(bw,
                            descrip=paste("bandwidth", sQuote("bw")),
                            fatal=TRUE) {
  if(!fatal)
    return(is.numeric(bw) && length(bw) == 1 && bw > 0)
  if(!is.numeric(bw)) 
    stop(paste(descrip, "was not numeric"), call.=FALSE)
  if(length(bw) != 1) 
    stop(paste(descrip, "was not a single number"), call.=FALSE)
  if(bw <= 0) 
    stop(paste(descrip, "was not a positive number"), call.=FALSE)
  return(TRUE)
}
