#Taken from package DMwR2
#Author: Luis Torgo ltorgo@dcc.fc.up.pt

nrLinesFile <- function (f) {
  if (.Platform$OS.type == "unix") 
    as.integer(strsplit(trimws(system(paste("wc -l", f), 
                                      intern = TRUE)), " ")[[1]][1])
  else stop("This function requires unix-based systems")
}
