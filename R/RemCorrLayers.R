#' Remove correlated layers
#'
#' @description For a pair of correlated raster layers, this function will remove one. The output stack will have layers with a correlation between them below the user specified threshold.
#'
#' @param in.stack a RasterStack object
#' @param t.cor the maximum allowed correlation between layers in the RasterStack. Must be a value between 0 and 1.
#'
#' @return A RasterStack of layers with correlation below the user specified correlation threshold.
#'
#' @examples
#' env.vars.uncorr <- RemCorrLayers(env.vars, 0.8)
#' @export


RemCorrLayers <- function(in.stack, t.cor){
  if(ncell(in.stack[[1]])<10000){
    ssize <- ncell(in.stack[[1]])
  } else {
    ssize <- 10000
  }
  sample.vals<-sampleRandom(in.stack, ssize)
  cor.matrix<-cor(sample.vals)
  i=1
  while (i <= ncol(cor.matrix)){
    rem.ind <- which(abs(cor.matrix[,i])>t.cor & abs(cor.matrix[,i])<1)
    if(length(rem.ind)==0){
      i <- i+1
      next
    } else {
      print(paste("removing ",colnames(cor.matrix)[rem.ind]))
      cor.matrix<-cor.matrix[-rem.ind, -rem.ind]
      i <- i+1
    }
  }
  out.stack<-in.stack[[which(names(in.stack)%in%colnames(cor.matrix))]]
  return(out.stack)
}
