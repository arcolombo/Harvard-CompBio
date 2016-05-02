#' This will download a 20 Mb
#' @param accessNumber a character string to download
#' @importFrom GEOquery getGEO
#' @export
#' @return expression set
accessingGEO<-function(accessNumber){

gse <- getGEO(accessNumber, GSEMatrix = TRUE)
#show(gse)


#filePaths = getGEOSuppFiles("GSE21653")
#filePaths
return(gse)
#library(GEOquery)

}
