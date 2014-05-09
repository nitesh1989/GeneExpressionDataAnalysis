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
### Check the data
date()


###################################################
### Clear the workspace, note the two "embedded" function
rm(list=ls())


###################################################
### Load libraries needed for the analysis
require(limma)
require(affy)
require(hgu95a.db)
require(annotate)


###################################################
### Load  gene expression
load("./objs/affyData.rda")

### Load  linear model results
load("./objs/linearModel.rda")


###################################################
### Lets create the Functional Gene Set (FGS) list
### We will get the data from the annotation medata package hgu95a.db
### Use the package name as a function for listing the annotation content
hgu95a()

### Each environment in a metadata package maps identifiers to annotation information
### The annotation information can be gene/identifiers centered, or annotation centered:
### for example "ID -> SYMBOL" is gene centered,
### while "PATHID -> genes" is annotation centered
### we need the second type of mapping: one FGS to all the genes of the FGS


##################################################
### This extract mapping between kegg and affy identifiers for this platform
kegg <- as.list(hgu95aPATH2PROBE)


##################################################
### This for GENE ONTOLOGY
go <- as.list(hgu95aGO2ALLPROBES)


##################################################
### The code above created a list where each element has name
### (the kegg identifiers) and contains a vector of identifiers, see below
length(kegg)
str(head(kegg))


##################################################
### To run a Wilcoxon rank-sum test we can use the geneSetTest() function in limma
### we need to identify the genes in the FGS, and give also a ranking statistics.
### To this end we will use the t-statistics obtained using the topTable() function

### The test on the fisrt FGS: default parameters
geneSetTest(tG2$ID%in%kegg[[1]],tG2$t)

### The test on the fisrt FGS: shift toward down-reguated genes
geneSetTest(tG2$ID%in%kegg[[1]],tG2$t,alternative="down")

### The test on the fisrt FGS: shift toward up-reguated genes
geneSetTest(tG2$ID%in%kegg[[1]],tG2$t,alternative="up")

### Visual representation
barcodeplot(tG2$t,tG2$ID%in%kegg[[1]])

### The test on all KEGG  FGS (after removing NA
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



### See the results
str(gse)
head(gse)

### The test on all KEGG  FGS (we can also test for NA inside the loop)
gse <- vector()
for (i in 1:length(kegg)){
	if (all(!tG2$ID%in%kegg[[i]])) {
		gse[i] <- NA
	}else{
		gse[i] <- geneSetTest(tG2$ID%in%kegg[[i]],tG2$t)
	}
	names(gse)[i] <- names(kegg)[i]
}

### Check how the gse vector
length(gse)
str(gse)


### This is the way we can visualize this: the first gene set
i=10
barcodeplot(tG2$t,tG2$ID%in%kegg[[i]])

### This is the way we can visualize this: the smallest pvalue
i=which.min(gse)
barcodeplot(tG2$t,tG2$ID%in%kegg[[i]])

### This is the way we can visualize this: the larest pvalue
i=which.max(gse)
barcodeplot(tG2$t,tG2$ID%in%kegg[[i]])


##################################################
### We can also use the geneSetTest() function differently
### You can learn about this function by typing
?geneSetTest

##################################################
### Or simply print the arguments of the function
args(geneSetTest)


##################################################
### We can compare our FGS to randomly generated lists of genes
### of the same size of our FGS of interest and compute the p-value
### as the proportion of times our FGS returns a smaller p-value
### compared to the randomly generated lists
geneSetTest(tG2$ID%in%kegg[[1]],tG2$t,alternative="up",type="t",
	    ranks.only=FALSE,nsim=1000)

geneSetTest(tG2$ID%in%kegg[[94]],tG2$t,alternative="down",type="t",
	    ranks.only=FALSE,nsim=1000)


##################################################
### limma enables you to run also self-contained tests as opposed to competitive tests
### To do this you can use the roast() function

### Extract the expression values
mat <- exprs(dat.rma)

### Select gene expression for the genes in the FGS
sel <- rownames(mat)%in%kegg[[1]]

##################################################
### For this analysis the design matrix must already account for the contrasts
newDmat <- data.frame(
		      Intercept=1,
		      CellType=dMat2[,1]+dMat2[,2] + (-1* (dMat2[,3]+dMat2[,4]))
		      )

### Run the test on one FGS: the contrast is the column index from the design matrix
roast( y=mat,index= sel, design=newDmat, contrast=2)


### The test on the first five FGS for KEGG
gse2 <- list()
for (i in 1:5){
	sel <- rownames(mat)%in%kegg[[i]]
	gse2[[i]] <- roast( mat,sel, design=newDmat, contrast=2)
	names(gse2)[i] <- names(kegg)[i]
}


### The test on the first five FGS for GO
gse.go <- list()
for (i in 1:5){
	sel <- rownames(mat)%in%go[[i]]
	gse.go[[i]] <- roast(mat,sel, design=newDmat, contrast=2)
	names(gse.go)[i] <- names(go)[i]
}


###################################################
### Now how we do find out which patwhays are enriched?

### Load the KEGG library
library(KEGG.db)
library(GO.db)

### If you do not have it install it from bioconductor
source("http://bioconductor.org/biocLite.R")
biocLite("KEGG.db")

### Check the content of the KEGG metadata package
KEGG()
GO()

### It appears that KEGGPATHID2NAME contains what we want
### so you can use mget() to obtain the patways names
pathNames <- mget(names(gse), KEGGPATHID2NAME, ifnotfound=NA)

goNames <- mget(names(gse.go),GOTERM,ifnotfound=NA)
goNames <- lapply(goNames,Term)

### Here we go
print(pathNames)

### Here we go
print(goNames)

### Retrieve the smallest p-values
pathNames[names(pathNames)==names(gse[which.min(gse)])]

### If they are in the same order
pathNames[which.min(gse)]

### Retrieve the smallest p-values
pathNames2 <- pathNames[length(pathNames):1]

### Always correct also if not in the same order
pathNames2[names(pathNames2)==names(gse[which.min(gse)])]

### Wrong if they are in the same order
pathNames2[which.min(gse)]

### Retrieve fgs with p-value less than 0.001
pathNames[gse<1e-5]

### Common genes
index <- which(gse<1e-2)
myFGS <- kegg[index]
sum(myFGS[[1]]%in%myFGS[[2]])

### Cross tabulation
all <- unique(unlist(myFGS))
table(fgs1=all%in%myFGS[[1]],fgs2=all%in%myFGS[[2]])

### Venn diagram
vennDiagram(data.frame(fgs1=all%in%myFGS[[1]],
		       fgs2=all%in%myFGS[[2]],fgs3=all%in%myFGS[[3]]))

#### Multiple testing correction
library(multtest)

### Correct p-values
adjp <- mt.rawp2adjp(gse)
adjp <- adjp$adjp[order(adjp$index),]

### Bonferroni only
adjp <- mt.rawp2adjp(gse,proc="Bonferroni")
adjp <- adjp$adjp[order(adjp$index),]

### Benjamini-Hochberg only
adjp <- mt.rawp2adjp(gse,proc="BH")
adjp <- adjp$adjp[order(adjp$index),]

### Make an histogram of the p-values
hist(adjp[,1], nclass=20)
hist(adjp[,2],nclass=20,add=TRUE,col="blue")



##################################################
### Session information
sessionInfo()

### Quit
q("no")

