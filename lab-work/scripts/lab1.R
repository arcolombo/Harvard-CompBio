library(RcompBio)
library(affy)
library(GEOquery)
library(drosophila2.db)


gse<- getGEO("GSE10940",GSEMatrix=TRUE)
 filePaths<-getGEOSuppFiles("GSE10940")

fns<-list.celfiles(path="~/Documents/Harvard-Liu/lab-work/GSE10940/",full.names=TRUE)
data.affy<-ReadAffy(filenames=fns)

 MAplot(data.affy,pairs=TRUE,which=c(1,2,3,4),plot.method="smoothScatter")
image(data.affy)
hist(data.affy) #unnormalized

boxplot(data.affy,col=seq(2,7,by=1))

 data.rma=rma(data.affy) #robust microarray average (normalization)
#how does this work?

 MAplot(data.rma,pairs=TRUE,which=c(1,2,3,4),plot.method="smoothScatter")


expr.rma<-exprs(data.rma)

boxplot(data.frame(expr.rma))

my_frame<-data.frame(expr.rma)

Annot <- data.frame(REFSEQ=sapply(contents(drosophila2REFSEQ), paste, collapse=", "),
SYMBOL=sapply(contents(drosophila2SYMBOL), paste, collapse=", "),
DESC=sapply(contents(drosophila2GENENAME), paste, collapse=", "))

 all<-merge(Annot,my_frame,by.x=0,by.y=0,all=T)


 control<-c(1,2,3,4,5,6)



mutants<-c(7:12)


#simple contrast

foldchange<-apply(expr.rma,1,function(x)mean(x[mutants])-mean(x[control]))

T.p.value=apply(expr.rma,1,function(x) t.test(x[mutants],x[control],var.equal=TRUE)$p.value)  #is the variance between groups really true?


fdr<-p.adjust(T.p.value,method="fdr")


 idx<-(which(all$SYMBOL!="NA"))
genes.refseq<-all$SYMBOL

upIdx<-which(fdr < 0.05 & foldchange>0)
downIdx<-which(fdr < 0.05 & foldchange < 0)

 genes.up<-all$SYMBOL[upIdx]

genes.down<-all$SYMBOL[downIdx]


  design<-c(rep(0,6),rep(1,6))
 design.mat<-model.matrix(~design)
 fit<-lmFit(expr.rma,design.mat)
 fit<-eBayes(fit)
topTable(fit,coef="design")


deGenes<-rownames(topTable(fit,coef="design",number=50))

#run a for loop to find the gene names in deGenes

#two factor experiment contrasting within groups and across


design2<-matrix(data=NA,nrow=12,ncol=4)

colnames(design2)<-c("WT.U","WT.S","Mu.U","Mu.S")


 rownames(design2)<-colnames(gse)


design2[1:3,1]<-1
 design2[4:12,1]<-0
 design2[,2]<-c(rep(0,3),rep(1,3),rep(0,6))
 c(rep(0,3),rep(0,3),rep(1,3),rep(0,3))->design2[,3]
c(rep(0,3),rep(0,3),rep(0,3),rep(1,3))->design2[,4]


 cont.matrix<-makeContrasts(MT.SvsU = WT.S- WT.U,
                            Mu.SvsU = Mu.S - Mu.U,
                         Diff= (Mu.S-Mu.U) - (WT.S-WT.U),
                     levels=design)


fitted<-lmFit(expr.rma,design2)
fitted2<-contrasts.fit(fitted,cont.matrix)
 fitted2<-eBayes(fitted2)


topGenes2<-topTable(fitted2,adjust="BH")




grp1<-var(expr.rma[,1:6])
grp2<-var(expr.rma[,7:12])

mean(grp1)
mean(grp2)
mean(grp1)-mean(grp2)

# on average sort of robust. for rma normalized.      raw values not so good.


