#' Sample rows from a text file
#' @description This function samples either a percentage or a number of rows from a text file at random without loading the entire file into R, thus reducing memory usage. Useful to get a sample of very large files. Modified from package DMwR2 by Luis Torgo ltorgo@@dcc.fc.up.pt.
#'
#' @param file The full path to the file to be read
#' @param percORn the proportion (if less than 1) or the number of rows to be sampled from the file.
#' @param nrLines integer. If known, the number rows in the file.
#' @param sep the field separator character. See ?read.delim for details.
#' @param dec the character used in the file for decimal points. See ?read.delim for details.
#' @param mxPerc maximum proportion (0 to 1) of rows to extract.
#'
#' @return a data frame with a subset of rows of the original file
#' @examples
#' occs <- SampleCSV("largeFile.csv", 10000)

SampleCSV <- function(file, percORn, nrLines, sep="\t", dec=".", mxPerc=1) {
  print("Running random subsampling module")
  if (.Platform$OS.type != "unix")
    stop("This function requires unix-based systems")

  if (missing(nrLines)) nrLines <- bigreadr::nlines(file)

  if (percORn < 1)
    if (percORn > mxPerc)
      stop("This function is not adequate for that big samples.")
  else percORn <- as.integer(percORn*nrLines)
  perc <- min(2*percORn/nrLines, mxPerc)
  print("Initialize random subsampling")
  system(paste0("perl -ne 'print if $. >= 1 && (rand() < ", perc,")' ",file, " > ", file, ".tmp.csv"))

  print("Read file into R")
  header <- read.delim(file, sep=sep, dec=dec, header=TRUE, nrows=1)
  dt <- read.delim(paste0(file,".tmp.csv"), sep=sep, dec=dec, col.names=colnames(header), as.is=TRUE, row.names=NULL,
                   nrows=percORn, quote="", fill=FALSE)
  file.remove(paste0(file,".tmp.csv"))
  if (nrow(dt) != percORn)
    warning(paste("Expecting",percORn,"rows, but got",nrow(dt)))
  return(dt)
}
