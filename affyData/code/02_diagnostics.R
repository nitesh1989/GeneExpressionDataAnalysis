###################################################
### Luigi Marchionni
### Edited by: Nitesh Turaga
### May 9, 2014
### Goal: from Affymetrix raw data stored in CEL files to differential gene expression


###################################################
### Getting  the working directory, in my case
getwd()


###################################################
### Setting the working directory, in my case
setwd("~/Documents/GeneExpressionDataAnalysis/affyData/")


###################################################
###Check the data
date()


###################################################
### Clear the workspace, note the two "embedded" function
rm(list=ls())


###################################################
### Load libraries needed for the analysis
require(affy)
require(hgu95a.db)


###################################################
### Load data previously save
load("./objs/affyData.rda")


###################################################
### The workspace content
ls()


###################################################
### Check classes
class(dat)
class(dat.rma)


###################################################
### Slots/elements contained in the two classes can be seen using:
### NOTE: SLOT NAMES USUALLY WORK AS ACCESSORS
slotNames(dat)
slotNames(dat.rma)

### This extract normalized and summarized gene expression
dat.expr <- exprs(dat.rma)
dim(dat.expr)

### This extract phenotypic information
pData(dat.rma)

### This extract features information
annotation(dat.rma)

### This extract perfect match probe intensities
dat.pm <- pm(dat)
dim(dat.pm)

### This extract mismatch probe intensities
dat.mm <- mm(dat)
dim(dat.mm)


##################################################
### We can now make an MA-plot
###################################################

###We first define a function to divide the plotting area in rows and columns
mypar <- function (nRow = 1, nCol = 1, ptsExp = 1) {
	par(mar = c(2, 2, 2, 1))
	par(oma = c(2, 1, 1, 1))
	par(mfrow = c(nRow, nCol))
	par(cex = ptsExp)
}


###################################################
### We can compute rows and columns number based on the nuber of samples
### You can do this also manually....
nc <- ceiling(sqrt(ncol(dat.rma)))
nr <- ceiling(ncol(dat.rma)/nc)


###################################################
### We can make a plot for the raw data, before normalization
mypar(nr,nc,0.5)
MAplot(dat[,],pairs=F,plot.method="smoothScatter")


###################################################
###We can save the plot to a file
bitmap(file="./figs/MApl.raw.png",width=20, height=10, res = 400)
mypar(nr,nc,0.5)
par(oma=c(2,1,1,1))
MAplot(dat[,],pairs=F,plot.method="smoothScatter")
dev.off()


###################################################
### MA-plot after normalization
mypar(nr,nc,0.5)
par(oma=c(2,1,1,1))
a <- exprs(dat.rma)
for (i in 1:ncol(a)){
  ma.plot(A=((a[,i]+apply(a,1,median))/2),
    M=(a[,i]-apply(a,1,median)),
    show.statistics=TRUE,cex.main=1,
    span=1/3, family.loess="gaussian", cex = 0.75,
    plot.method="smoothScatter",add.loess = TRUE,
    lwd = 2, lty = 2, loess.col = "red",ylim=c(-6,6),
    main=paste(colnames(a)[i],
	  '\n vs pseudo-median reference chip'))
}


###################################################
### MA-plot after normalization: saved to a file
bitmap(file='./figs/MApl.rma.png',width=20,height=10,res=400)
mypar(nr,nc,0.5)
par(oma=c(2,1,1,1))
a <- exprs(dat.rma)
for (i in 1:ncol(a)){
  ma.plot(A=((a[,i]+apply(a,1,median))/2),
    M=(a[,i]-apply(a,1,median)),
    show.statistics=TRUE,cex.main=1,
    span=1/3, family.loess="gaussian", cex = 0.75,
    plot.method="smoothScatter",add.loess = TRUE,
    lwd = 2, lty = 2, loess.col = "red",ylim=c(-6,6),
    main=paste(colnames(a)[i],
	  'vs pseudo-median reference chip'))
}
dev.off()


###################################################
### Diagnostic Plots: 2D-Image
bitmap(file="./figs/ima2D.png",width=2*nr, height=2*nc, res = 144)
mypar(nr,nc,1) ; par(mar=c(1,1,1,1))
image(dat[,])
dev.off()


###################################################
### Log2 intensities boxplots: before and after normalization
bitmap(file="./figs/exprBoxplot.png",width=15, height=15, res = 144)
mypar(2,1,0.75)
par(mar=c(10,1,1,1))
### Before
boxplot(log2(exprs(dat))~col(exprs(dat)),
	names=pData(dat.rma)[,"Row.names"],
	las=2,col='blue')
### After
boxplot(exprs(dat.rma)~col(exprs(dat.rma)),
	names=pData(dat.rma)[,"Row.names"],
	las=2,col='blue')
### Close the device
dev.off()



###################################################
### Compute RNA degradation
deg <- AffyRNAdeg(dat)

### Make RNA degradation plot
bitmap(file="./figs/degRNA.png",width=5, height=5, res = 1000)
mypar(2,2,0.75)
plotAffyRNAdeg(deg,cols=c(1:ncol(dat.rma)))
plot(density(deg$slope),main="Slope")
boxplot(deg$slope,main="Slope")
### Close the device
dev.off()


### Find the scan date
dateOfScan <- protocolData(dat.rma)@data$ScanDate
dateOfScan <- gsub(" .+", "", dateOfScan)
boxplot(log2(exprs(dat)), col=1+as.numeric(dateOfScan))


##################################################
### Session information
sessionInfo()

### Quit
q("no")
