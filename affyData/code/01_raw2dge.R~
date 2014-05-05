###################################################
### Luigi Marchionni
### May 10, 2013

### Goal: from Affymetrix raw data stored in CEL files to differential gene expression


###################################################
### Getting  the working directory, in my case
getwd()


###################################################
### Setting the working directory, in my case

setwd("~/Documents/Feb142014-backup/Documents/JHMI-Research/Luigi_GeneExpression_TA/affyData")
#setwd("~/3EDUCATION/myRtutorials/affyData")


###################################################
### Check the data
date()


###################################################
### We are going to use the following packages from bioconductor:
### "affy", "limma", "hgu95a.db", "annotate"
### To install the missing packages from Bioconductor use biocLite()

### First you need to get the list of available packages
installedPckgs <- installed.packages()[,"Package"]

### Here we define the list of desired libraries
pckgListBIOC <- c("affy", "limma", "hgu95a.db", "annotate")

### We source the biocLite.R script from the Bioconductor website
source("http://bioconductor.org/biocLite.R")

### Load the packages, or install them from Bioconductor if needed
for (pckg in pckgListBIOC) {
	if (! pckg %in% installedPckgs) biocLite(pckg)
	require(pckg, character.only=TRUE)
}


###################################################
### Clear the workspace, note the two "embedded" functions
rm(list=ls())


###################################################
### Load the libraries needed for the analysis, note the use of "require"
require(affy)
require(limma)
require(hgu95a.db)
require(annotate)


##################################################
##################################################
### Reading and preprocessing the expression data
##################################################
##################################################


###################################################
### Load data from CEL files using the ReadAffy function
dat <- ReadAffy(celfile.path="./data/raw/", compress=TRUE)


###################################################
### Let's have a look at the dat object

### Just print
print(dat)

### Class
class(dat)

### You can find information about this class as usual
?AffyBatch

### Dimensions
dim(dat)

### Structure
str(dat)

### Information about the platform
annotation(dat)

### Information about the features
head(featureData(dat))

### Information about the samples
pData(dat)


###################################################
### Read the phenotype data associated with the experiment
### Note the use of colClasses
pheno <- read.table('./data/targets.txt',sep='\t',header=T,
		    colClasses=c("character", "character",
		    "factor","factor","character"))


###################################################
### Let's have a look into the phenotypic information
str(pheno)

### Let's get summarize the data
summary(pheno)

### Factor are very useful!
table(pheno$CancerType, pheno$CellType)

plot(pheno$CancerType, pheno$CellType)


###################################################
### We can add phenotypic information to the AffyBatch object
if(all(rownames(pData(dat))==pheno$Name)){
	print("You're good, go ahead!")
	pheno <- merge(pData(dat), pheno, by.x=0, by.y='Name', sort=FALSE)
	rownames(pheno) <- pheno[,"Row.names"]
	pData(dat) <- pheno
}else{
	print("Check order of rows in pData(object) and pheno data provided!")
}

### Now pData(dat) is the merged object
pData(dat)


###################################################
### We can now normalize the data, for instance the  rma() function
### Alternatives you want to look into are frma() and gcrma()
dat.rma <- rma(dat)


###################################################
### Check the new object after normalization
class(dat.rma)

### Nevertheless the following methods still work:
pData(dat.rma)
annotation(dat.rma)
head(featureData(dat.rma))
dim(dat.rma)

###################################################
### Now we can save the two objects before going ahead
save(dat, dat.rma, file="./objs/affyData.rda")



##################################################
##################################################
### Differential Gene Expression Analysis
##################################################
##################################################


###################################################
### To create the design matrix  to fit the linear model we start
### combining all levels of our factors of interest
groups <- factor(paste(pheno$CancerType, pheno$CellType, sep="."))

### Count the number of sample per group
table(groups)

### Create the design matrix
dMat <- model.matrix(~0+groups)
colnames(dMat) <- levels(groups)

### This is the design matrix
dMat

###In reality we could just use the informations in
## the pheno$CancerType factor, which is equivalent
table(pheno$CancerType)
dMat2 <- model.matrix(~0+pheno$CancerType)
colnames(dMat2) <- levels(pheno$CancerType)

### This is the second design matrix
dMat2

###As you can see they are exactly the same
dMat==dMat2


###################################################
### Let's create now a contrast matrix to extract differences between coefficients

#### A contrast matrix using the first design matrix
cMat <- makeContrasts(levels=colnames(dMat),
			 CellType=( ((colon.Epithelial + breast.Epithelial)/2) -
				   ((Burkitt.Lymphocytes + DLBCL.Lymphocytes)/2) ),
			 CancerType.Epithelial=( colon.Epithelial - breast.Epithelial ),
			 CancerType.Lyphocytes=( Burkitt.Lymphocytes-DLBCL.Lymphocytes )
			 )

### All coefficients that were combine contribute equally:
### the columns sum must be equal to 0
cMat
apply(cMat, 2, sum)


### A contrast matrix using the second design matrix
cMat2 <- makeContrasts(levels=colnames(dMat2),
		       CellType=( ((colon + breast)/2) - ((Burkitt + DLBCL)/2) ),
		       CancerType.Epithelial=( colon - breast ),
		       CancerType.Lyphocytes=( Burkitt - DLBCL )
			 )

### Also in this case all the coefficients that were combine must contribute equally:
### the columns sum again must be  equal to 0
cMat2
apply(cMat2, 2, sum)

### The two contrasts matrices are the same
cMat==cMat2


###################################################
### Fit the model using the first design and contrast matrices
### You might want to learn about lmFit() and contrasts.fit()
fit.ls <- lmFit(dat.rma, dMat, method="ls")
fit.ls <- contrasts.fit(fit.ls, cMat)


###################################################
### Fit the model using the second design matrix and contrast matrices
fit2.ls <- lmFit(dat.rma, dMat2, method="ls")
fit2.ls <- contrasts.fit(fit2.ls, cMat2)


###################################################
### Moderation of standard errors using empirical Bayes
eb.ls <- eBayes(fit.ls, proportion=0.01)


###################################################
### Moderation of standard errors using empirical Bayes
eb2.ls <- eBayes(fit2.ls, proportion=0.01)


###################################################
### With the following you can write the complete linear model analysis
### results to a file in a tabular format
write.fit(eb.ls, file="./text/anovaTable.txt")


###################################################
### You can retrieve the top differentially expressed genes for a specific
### contrast in the contrast matrix by name
tG <- topTable(eb.ls, coef="CellType", number=10, resort="logFC")
tG <- tG[order(tG$ID),]


###################################################
### You can retrieve the top differentially expressed genes for a specific
### contrast in the contrast matrix by index
tG2 <- topTable(eb2.ls, coef=1, number=10, resort="logFC")
tG2 <- tG2[order(tG2$ID),]


##################################################
### The same genes were identified
all(tG == tG2)


###################################################
### We can also return ALL the features investigated by using Inf
tG2 <- topTable(eb2.ls, coef=1, number=Inf, resort="logFC")


###################################################
### As you can see we need to add gene annotation!
str(tG2)


##################################################
##################################################
### Retrieve and add the annotation
##################################################
##################################################

###################################################
### Retrieve and add gene SYMBOLS, ENTREZID, and GENENAME
### Using the "old" hard-way method based on mget():
egid <- mget(tG2$ID, hgu95aENTREZID, ifnotfound=NA)

sym <- mget(tG2$ID, hgu95aSYMBOL, ifnotfound=NA)
nm <- mget(tG2$ID, hgu95aGENENAME, ifnotfound=NA)


##################################################
### Add the annnotation to differential gene expression results
tG2ann <- cbind(EGID=unlist(egid), SYMBOL=unlist(sym), GeneName=unlist(nm),
	    tG2, stringsAsFactors=FALSE)

### Let's check the results now
str(tG2ann)


###################################################
### Retrieve and add gene SYMBOLS, ENTREZID, and GENENAME
### Using the "new" and simple method based on select():
ann <- select(x=hgu95a.db, keys=tG2$ID,
	      cols=c("ENTREZID", "SYMBOL", "GENENAME"))

### Check the dimensions of the annotation the results data.frames
dim(ann)
dim(tG2)


##################################################
### Add the annnotation to differential gene expression results
tG2ann2 <- merge(x=ann, y=tG2, by.x="PROBEID", by.y="ID",
		 all.x=TRUE, all.y=FALSE)

### Let's check the results now
str(tG2ann2)


##################################################
### We can save the results as usual for later use
save(tG2, fit2.ls, eb2.ls, dMat2, cMat2, file="objs/linearModel.rda")


##################################################
### We can also filter the genes and create a report page
#tmp <- tG2[ tG2$B > 10 , ]
tmp <- tG2ann2[ tG2ann2$adj.P.Val < 1e-3 , ]


### Reorder the columns to have ENTREZID as the first column:
### hyperlinks to the NCBI ENTREZ GENE database will be generated
tmp <- tmp[ , c(2, 1, 3:ncol(tmp)) ]


### Write the results to an html file
htmlpage(as.data.frame(tmp[ , 1]),
	 filename="html/result.html",
	 othernames=as.data.frame(tmp[ , 2:ncol(tmp) ]),
	 table.head=colnames(tmp), digits=4)



##################################################
##################################################
### A few useful plots
##################################################
##################################################

##################################################
### We can make a volcano plot
volcanoplot(eb2.ls)


##################################################
### We can Retrieve the original gene expression for the interesting genes

### All the expression values
mat <- exprs(dat.rma)
dim(mat)

### Filter to the interesting ones (differentially expressed genes)
mat <- mat[rownames(mat)%in%tmp$PROBEID,]
dim(mat)


##################################################
### Generate a heatmap()
heatmap(mat, scale='none', na.rm=TRUE, margins=c(15,7),
	distfun=function(x) {dist(x, method='euclidian')},
	hclustfun=function(x) {hclust(x, method='average')},
	ColSideColors=as.character(2+as.numeric(pheno$CellType)))


### Add symbol as rownames using mget()
colnames(mat) <- pheno$CellLine
mySym <- mget(rownames(mat), hgu95aSYMBOL, ifnotfound=NA)
rownames(mat) <- unlist(mySym)


### Generated the annotated heatmap with heatmap()
heatmap(mat, scale='row', na.rm=TRUE, margins=c(15,7),
	distfun=function(x) {dist(x, method='euclidian')},
	hclustfun=function(x) {hclust(x, method='average')},
	ColSideColors=as.character(2+as.numeric(pheno$CellType)))



##################################################
### Session information
sessionInfo()

### Quit
q("no")
