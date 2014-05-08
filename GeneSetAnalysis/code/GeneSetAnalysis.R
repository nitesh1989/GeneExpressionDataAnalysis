### R code from vignette source 'GeneSetAnalysis.Rnw'
### Encoding: UTF-8

setwd("~/Documents/GitHub/GeneExpressionDataAnalysis/GeneSetAnalysis/")

###################################################
### code chunk number 1: start
###################################################
options(width=50)
rm(list=ls())


###################################################
### code chunk number 2: getPackagesBioc
###################################################
###Source the biocLite.R script from Bioconductor
source("http://bioconductor.org/biocLite.R")
###Get the list of available packages
installedPckgs <- installed.packages()[,"Package"]
###Define the list of desired libraries
pckgListBIOC <- c("BiocGenerics", "Biobase", "limma", "RTopper", 
		  "org.Hs.eg.db", "AnnotationDbi", "annotate", 
		  "multtest", "KEGG.db", "GO.db")
###Load the packages, install them from Bioconductor if needed
for (pckg in pckgListBIOC) {
	if (! pckg %in% installedPckgs) {
		biocLite(pckg, suppressUpdates=TRUE, ask=FALSE)
		}
	require(pckg, character.only=TRUE)
}


###################################################
### code chunk number 3: loadData
###################################################
require(RTopper)
data(sepScores)
ls()
class(sepScores)
names(sepScores)
sapply(sepScores, class)
sapply(sepScores,dim)


###################################################
### code chunk number 4: metaData
###################################################
require(org.Hs.eg.db)
org.Hs.eg()


###################################################
### code chunk number 5: listFGS
###################################################
kegg <- as.list(org.Hs.egPATH2EG)
length(kegg)
str(kegg[1:5])
names(kegg)[1:5]
go <- as.list(org.Hs.egGO2ALLEGS)
length(go)
str(go[1:5])
names(go)[1:5]


###################################################
### code chunk number 6: convertIDs
###################################################
numberOfFGSkegg <- 200
kegg <- lapply(kegg[sample(1:length(kegg),numberOfFGSkegg)],
               function(x) unique(unlist(mget(x,org.Hs.egSYMBOL))))
str(kegg[1:5])
### Process GO: keep only Biological Process
length(go)
go <- go[ Ontology(names(go)) == "BP" ]
length(go)
numberOfFGSgo <- 200
go <- lapply(go[sample(1:length(go),numberOfFGSgo)],
             function(x) unique(unlist(mget(x,org.Hs.egSYMBOL))))
str(go[1:5])


###################################################
### code chunk number 7: annotateFGS
###################################################
require(KEGG.db)
KEGG()
names(kegg) <- paste(names(kegg),unlist(mget(names(kegg),KEGGPATHID2NAME)),sep=".")
head(names(kegg), n=10)


###################################################
### code chunk number 8: listFGS
###################################################
require(GO.db)
GO()
names(go) <- paste(names(go),Term(names(go)),sep=".")
head(names(go), n=10)


###################################################
### code chunk number 9: listFGS
###################################################
fgsList <- list(go=go,kegg=kegg)
str(fgsList$go[1:5])
str(fgsList$kegg[1:5])


###################################################
### code chunk number 10: runGSEbatchArgs
###################################################
args(runBatchGSE)


###################################################
### code chunk number 11: runBatchGSE.separate
###################################################
gseABS.sep <- runBatchGSE(dataList=sepScores, fgsList=fgsList)
gseABS.sep <- runBatchGSE(dataList=sepScores, fgsList=fgsList,
				 absolute=TRUE, type="f", alternative="mixed")


###################################################
### code chunk number 12: runBatchGSE.separate2
###################################################
gseUP.sep <- runBatchGSE(dataList=sepScores, fgsList=fgsList,
				 absolute=FALSE, type="t", alternative="up")
gseDW.sep <- runBatchGSE(dataList=sepScores, fgsList=fgsList,
				 absolute=FALSE, type="t", alternative="down")
gseBOTH.sep <- runBatchGSE(dataList=sepScores, fgsList=fgsList,
				 absolute=FALSE, type="t", alternative="either")


###################################################
### code chunk number 13: runBatchGSE.int3
###################################################
gseABSsim.sep <- runBatchGSE(dataList=sepScores, fgsList=fgsList,
				    absolute=TRUE, type="f", alternative="mixed",
				    ranks.only=FALSE, nsim=1000)
gseUPsim.sep <- runBatchGSE(dataList=sepScores, fgsList=fgsList,
				    absolute=FALSE, type="t", alternative="up",
				    ranks.only=FALSE, nsim=1000)


###################################################
### code chunk number 14: runBatchGSE.format1
###################################################
str(gseUP.sep[1:5])
head(gseABSsim.sep$dat.affy$go)
head(gseABSsim.sep$dat.affy$kegg)


###################################################
### code chunk number 15: runBatchGSE.altFunc
###################################################
require(limma)
gseUP.sep.2 <- runBatchGSE(dataList=sepScores, fgsList=fgsList,
				 absolute=FALSE, gseFunc=wilcoxGST, alternative="up")


###################################################
### code chunk number 16: runBatchGSE.format2
###################################################
str(gseUP.sep.2)
all(gseUP.sep.2$go==gseUP.sep$go)


###################################################
### code chunk number 17: runBatchGSE.altFunc2
###################################################
gseFunc <- function (selected, statistics, threshold) {
	diffExpGenes <- statistics > threshold
	tab <- table(diffExpGenes, selected)
	pVal <- fisher.test(tab)[["p.value"]]
	}
gseUP.sep.3 <- runBatchGSE(dataList=sepScores, fgsList=fgsList,
				 absolute=FALSE, gseFunc=gseFunc, threshold=7.5)


###################################################
### code chunk number 18: runBatchGSE.format3
###################################################
str(gseUP.sep.3)
head(data.frame(fisher=gseUP.sep.3$dat.affy$kegg,
                wilcoxon=gseUP.sep$dat.affy$kegg))


###################################################
### code chunk number 19: adjustP
###################################################
gseABS.sep.BH <- adjustPvalGSE(gseABS.sep)
gseABS.sep.holm <- adjustPvalGSE(gseABS.sep, proc = "Holm")


###################################################
### code chunk number 20: adjusted.format
###################################################
names(gseABS.sep.BH)
names(gseABS.sep.holm)
str(gseABS.sep.BH$dat.affy)
str(gseABS.sep.holm$dat.affy)
head(gseABS.sep.BH$dat.affy$go, n=10)
head(gseABS.sep.holm$dat.affy$go, n=10)


###################################################
### code chunk number 21: figure004
###################################################
hist(gseABS.sep.BH$dat.affy$go[ , "rawp"], nclass=10, col="orange",
     main="P-values")


###################################################
### code chunk number 22: figure005
###################################################
hist(gseABS.sep.BH$dat.affy$go[ , "BH"], nclass=10, col="orange",
     main="Q-values")


###################################################
### code chunk number 23: readAndProcessData
###################################################
autismGenes <- read.table("./data/101symbols_v02.txt", sep="\t", header=TRUE,
                          colClasses="character")
str(autismGenes)


###################################################
### code chunk number 24: newFunctionForLists
###################################################
formatGeneList <- function (geneList, allGenes) {
  stats <- 1* allGenes %in% geneList
  names(stats) <- allGenes
  out <- list(GeneList = stats)
	}


###################################################
### code chunk number 25: formatTheGeneLsit
###################################################
allGenes <- unique(unlist(as.list(org.Hs.egSYMBOL)))
str(allGenes)
geneMembership <- formatGeneList(autismGenes$SYMBOL, allGenes)
str(geneMembership)
table(geneMembership$GeneList)


###################################################
### code chunk number 26: runBatchGSE.altFunc3
###################################################
gseMyList <- runBatchGSE(dataList=geneMembership, fgsList=fgsList,
				 absolute=FALSE, gseFunc=gseFunc, threshold=0.5)


###################################################
### code chunk number 27: showTheResults
###################################################
str(gseMyList)
str(gseMyList$GeneList$kegg)
summary(gseMyList$GeneList$kegg)
table(is.na(gseMyList$GeneList$kegg))
table(gseMyList$GeneList$kegg < 0.01)
table(gseMyList$GeneList$go < 0.01)


###################################################
### code chunk number 28: adjustP.2
###################################################
gseMyList.BH <- adjustPvalGSE(gseMyList)
gseMyList.BH$GeneList$go[ gseMyList.BH$GeneList$go[, "BH"] < 0.01 , ]
gseMyList.BH$GeneList$kegg[ gseMyList.BH$GeneList$kegg[, "BH"] < 0.01 , ]


###################################################
### code chunk number 29: sessioInfo
###################################################
sessionInfo()


