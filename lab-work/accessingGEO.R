### This will download a 20 Mb
gse <- getGEO("GSE21653", GSEMatrix = TRUE)
show(gse)


filePaths = getGEOSuppFiles("GSE21653")
filePaths

library(GEOquery)
