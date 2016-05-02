#' function to read in raw affy and access affy library

#' @importFrom affy list.celfiles
#' @importFrom affy ReadAffy
fectchRawAffy<-function(affyPath, normalize=TRUE,findAnnotation=TRUE)

stopifnot(file.exists(affyPath)=="TRUE")

fns<-list.celfiles(path="~/Documents/Harvard-Liu/lab-work/GSE10940/",full.names=TRUE)
stopifnot(all(file.exists(fns)==TRUE))
data.affy<-ReadAffy(filenames=fns)
final.affy<-data.affy

if(normalize==TRUE) {
#robust multiarray average
final.affy<-rma(data.affy)

    if(findAnnotation==TRUE){


    }

}




return(final.affy)
