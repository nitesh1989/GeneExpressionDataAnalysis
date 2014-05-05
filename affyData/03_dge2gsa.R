
## ----getwd---------------------------------------------------------------
getwd()


## ----setwd,eval=TRUE-----------------------------------------------------
setwd("~/Documents/GeneExpressionDataAnalysis/affyData/")


## ----checkDate-----------------------------------------------------------
date()


## ----clearWorkSpace------------------------------------------------------
rm(list=ls())


## ----loadLibraries-------------------------------------------------------
require(affy)
require(limma)
require(hgu95a.db)
require(annotate)


## ----loadData------------------------------------------------------------
# Gene expression data
load("./objs/affyData.rda")

# Load  linear model results
load("./objs/linearModel.rda")


## ----FGS-----------------------------------------------------------------
hgu95a()


## ----extractMapping------------------------------------------------------
kegg <- as.list(hgu95aPATH2PROBE)


## ----geneOntology--------------------------------------------------------
go <- as.list(hgu95aGO2ALLPROBES)


## ----explainList---------------------------------------------------------
class(kegg)
length(kegg)
str(head(kegg))


## ------------------------------------------------------------------------

geneSetTest(tG2$ID%in%kegg[[1]],tG2$t)



## ------------------------------------------------------------------------
geneSetTest(tG2$ID%in%kegg[[1]],tG2$t,alternative="down")


## ------------------------------------------------------------------------
geneSetTest(tG2$ID%in%kegg[[1]],tG2$t,alternative="up")


## ----barCodePlot---------------------------------------------------------
barcodeplot(tG2$ID%in%kegg[[1]],tG2$t)


## ----testAllFGS----------------------------------------------------------
kegg2 <- kegg[!is.na(kegg)]
gse <- vector()
for (i in 1:length(kegg2)){
  gse[i] <- geneSetTest(tG2$ID%in%kegg2[[i]],tG2$t)
	names(gse)[i] <- names(kegg2)[i]
}

### This is an alternative
tmp <- sapply(kegg2, function(gs, stat, ids) {
	geneSetTest(ids%in%gs,stat)
}, stat=tG2$t, ids=tG2$ID)


## ----checkResults--------------------------------------------------------
str(gse)
head(gse)


## ----removeNAInsideLoop--------------------------------------------------
gse <- vector()
for (i in 1:length(kegg)){
  #check if all values are NOT true
	if (all(!tG2$ID%in%kegg[[i]])) {
		gse[i] <- NA
	}else{
		gse[i] <- geneSetTest(tG2$ID%in%kegg[[i]],tG2$t)
	}
	names(gse)[i] <- names(kegg)[i]
}


## ----checkLength---------------------------------------------------------
length(gse)
str(gse)


## ----firstGeneSet--------------------------------------------------------
i=10
barcodeplot(tG2$ID%in%kegg[[i]],tG2$t)


## ----visualizeSmallestPValue---------------------------------------------
i=which.min(gse)
barcodeplot(tG2$ID%in%kegg[[i]],tG2$t)


## ----visualizeLargestPValue----------------------------------------------
i=which.max(gse)
barcodeplot(tG2$ID%in%kegg[[i]],tG2$t)


## ----eval=FALSE----------------------------------------------------------
## ?geneSetTest
## 
## # Or simply print the arguments of the function
## args(geneSetTest)


## ------------------------------------------------------------------------
geneSetTest(tG2$ID%in%kegg[[1]],tG2$t,alternative="up",type="t",
	    ranks.only=FALSE,nsim=1000)

geneSetTest(tG2$ID%in%kegg[[94]],tG2$t,alternative="down",type="t",
	    ranks.only=FALSE,nsim=1000)


## ------------------------------------------------------------------------
mat <- exprs(dat.rma)


## ------------------------------------------------------------------------
sel <- rownames(mat)%in%kegg[[1]]


## ------------------------------------------------------------------------
newDmat <- data.frame(
		      Intercept=1,
		      CellType=dMat2[,1]+dMat2[,2] + (-1* (dMat2[,3]+dMat2[,4]))
		      )


## ------------------------------------------------------------------------
roast(iset= sel, y=mat, design=newDmat, contrast=2)


## ------------------------------------------------------------------------
gse2 <- list()
for (i in 1:5){
	sel <- rownames(mat)%in%kegg[[i]]
	gse2[[i]] <- roast(sel, mat, design=newDmat, contrast=2)
	names(gse2)[i] <- names(kegg)[i]
}


## ------------------------------------------------------------------------
gse.go <- list()
for (i in 1:5){
	sel <- rownames(mat)%in%go[[i]]
	gse.go[[i]] <- roast(sel, mat, design=newDmat, contrast=2)
	names(gse.go)[i] <- names(go)[i]
}


## ------------------------------------------------------------------------
library(KEGG.db)
library(GO.db)


## ------------------------------------------------------------------------
source("http://bioconductor.org/biocLite.R")
biocLite("KEGG.db")


## ------------------------------------------------------------------------
KEGG()
GO()


## ------------------------------------------------------------------------
pathNames <- mget(names(gse), KEGGPATHID2NAME, ifnotfound=NA)

goNames <- mget(names(gse.go),GOTERM,ifnotfound=NA)
goNames <- lapply(goNames,Term)


## ------------------------------------------------------------------------
print(pathNames)

print(goNames)


## ----smallestPval--------------------------------------------------------
pathNames[names(pathNames)==names(gse[which.min(gse)])]


## ------------------------------------------------------------------------
pathNames[which.min(gse)]


## ------------------------------------------------------------------------
pathNames2 <- pathNames[length(pathNames):1]


## ------------------------------------------------------------------------
pathNames2[names(pathNames2)==names(gse[which.min(gse)])]


## ------------------------------------------------------------------------
pathNames2[which.min(gse)]


## ------------------------------------------------------------------------
pathNames[gse<1e-5]


## ------------------------------------------------------------------------
index <- which(gse<1e-2)
myFGS <- kegg[index]
sum(myFGS[[1]]%in%myFGS[[2]])


## ------------------------------------------------------------------------
all <- unique(unlist(myFGS))
table(fgs1=all%in%myFGS[[1]],fgs2=all%in%myFGS[[2]])


## ------------------------------------------------------------------------
vennDiagram(data.frame(fgs1=all%in%myFGS[[1]],
		       fgs2=all%in%myFGS[[2]],fgs3=all%in%myFGS[[3]]))


## ----multipleTestingCorrection-------------------------------------------
library(multtest)

# Correct p-values
adjp <- mt.rawp2adjp(gse)
adjp <- adjp$adjp[order(adjp$index),]

# Bonferroni only
adjp <- mt.rawp2adjp(gse,proc="Bonferroni")
adjp <- adjp$adjp[order(adjp$index),]

# Benjamini-Hochberg only
adjp <- mt.rawp2adjp(gse,proc="BH")
adjp <- adjp$adjp[order(adjp$index),]


## ------------------------------------------------------------------------
hist(adjp[,1], nclass=20)
hist(adjp[,2],nclass=20,add=TRUE,col="blue")


## ----sessionInformation--------------------------------------------------
sessionInfo()


