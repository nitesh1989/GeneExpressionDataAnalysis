__Author: Luigi Marchionni, May 6th 2014__


Goal: Gene Set Analysis
========================================================


### Some Basics

Getting  the working directory, in my case


```r
getwd()
```

```
## [1] "/Users/niteshturaga/Documents/GeneExpressionDataAnalysis/affyData"
```


Setting the working directory, in my case:

```r
setwd("~/Documents/GeneExpressionDataAnalysis/GeneSetAnalysis/")
```


Check the date


```r
date()
```

```
## [1] "Mon May  5 12:08:59 2014"
```


Clear the workspace, note the two "embedded" functions

```r
rm(list = ls())
```


Bioconductor
-------------

The following R packages were used to perform our analyses 

1. **BiocGenerics**: this library implements several generic classes and methods to work with {\Bioc} packages and can be obtained from Bioconductor.

2. **Biobase**: this library implements classes and methods to work with genomic data and can be obtained from {\Bioc};

3. **limma**: this library provide advanced methods to perform gene expression analysis from raw data to gene set analysis and can be obtained from Bioconductor.

4. **RTopper**: this library implements methods to perform integrated gene set analysis across datasets and platforms and can be obtained from Bioconductor.

5. **multtest**: this library implements methods to perform multiple test correction and can be obtained from Bioconductor.

6. **org.Hs.eg.db**: this library contains human genes annotation and can be obtained from Bioconductor.

7. **KEGG.db**: this library contains KEGG pathway annotation and can be obtained from Bioconductor.

8. **reactome.db**: this library contains Reactome pathway annotation and can be obtained from Bioconductor.

9. **GO.db**: this library contains Gene Ontoogy annotation and can be obtained from Bioconductor.

10. **AnnotationDBI** and **annotate**: these libraries contain classes and methods to access and manipulated annotation packages and can be obtained from Bioconductor.


Source the biocLite.R script from Bioconductor

```r
source("http://bioconductor.org/biocLite.R")
```

```
## Bioconductor version 2.13 (BiocInstaller 1.12.1), ?biocLite for help
## A newer version of Bioconductor is available after installing a new
##   version of R, ?BiocUpgrade for help
```

```r
### Get the list of available packages
installedPckgs <- installed.packages()[, "Package"]
### Define the list of desired libraries
pckgListBIOC <- c("BiocGenerics", "Biobase", "limma", "RTopper", "org.Hs.eg.db", 
    "AnnotationDbi", "annotate", "multtest", "KEGG.db", "GO.db")
### Load the packages, install them from Bioconductor if needed
for (pckg in pckgListBIOC) {
    if (!pckg %in% installedPckgs) {
        biocLite(pckg, suppressUpdates = TRUE, ask = FALSE)
    }
    require(pckg, character.only = TRUE)
}
```

```
## Loading required package: BiocGenerics
## Loading required package: parallel
## 
## Attaching package: 'BiocGenerics'
## 
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
## 
## The following object is masked from 'package:stats':
## 
##     xtabs
## 
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, as.vector, cbind,
##     colnames, duplicated, eval, evalq, Filter, Find, get,
##     intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rep.int, rownames, sapply, setdiff, sort,
##     table, tapply, union, unique, unlist
## 
## Loading required package: Biobase
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
## 
## Loading required package: limma
## 
## Attaching package: 'limma'
## 
## The following object is masked from 'package:BiocGenerics':
## 
##     plotMA
## 
## Loading required package: RTopper
## Loading required package: org.Hs.eg.db
## Loading required package: AnnotationDbi
## Loading required package: DBI
## 
## Loading required package: annotate
## Loading required package: multtest
## Loading required package: KEGG.db
## 
## KEGG.db contains mappings based on older data because the original
##   resource was removed from the the public domain before the most
##   recent update was produced. This package should now be
##   considered deprecated and future versions of Bioconductor may
##   not have it available.  Users who want more current data are
##   encouraged to look at the KEGGREST or reactome.db packages
## 
## Loading required package: GO.db
```


### code chunk number 3: Load Data and packages required

```r
require(RTopper)
data(sepScores)
ls()
```

```
## [1] "installedPckgs" "pckg"           "pckgListBIOC"   "sepScores"
```

```r
class(sepScores)
```

```
## [1] "list"
```

```r
names(sepScores)
```

```
## [1] "dat.affy"       "dat.agilent"    "dat.cnvHarvard" "dat.cnvMskcc"
```

```r
sapply(sepScores, class)
```

```
##       dat.affy    dat.agilent dat.cnvHarvard   dat.cnvMskcc 
##      "numeric"      "numeric"      "numeric"      "numeric"
```

```r
sapply(sepScores, dim)
```

```
## $dat.affy
## NULL
## 
## $dat.agilent
## NULL
## 
## $dat.cnvHarvard
## NULL
## 
## $dat.cnvMskcc
## NULL
```


### code chunk number 4: metaData

```r
require(org.Hs.eg.db)
org.Hs.eg()
```

```
## Quality control information for org.Hs.eg:
## 
## 
## This package has the following mappings:
## 
## org.Hs.egACCNUM has 31277 mapped keys (of 46265 keys)
## org.Hs.egACCNUM2EG has 629908 mapped keys (of 629908 keys)
## org.Hs.egALIAS2EG has 99696 mapped keys (of 99696 keys)
## org.Hs.egCHR has 45772 mapped keys (of 46265 keys)
## org.Hs.egCHRLENGTHS has 93 mapped keys (of 93 keys)
## org.Hs.egCHRLOC has 23648 mapped keys (of 46265 keys)
## org.Hs.egCHRLOCEND has 23648 mapped keys (of 46265 keys)
## org.Hs.egENSEMBL has 25125 mapped keys (of 46265 keys)
## org.Hs.egENSEMBL2EG has 27388 mapped keys (of 27388 keys)
## org.Hs.egENSEMBLPROT has 19731 mapped keys (of 46265 keys)
## org.Hs.egENSEMBLPROT2EG has 101504 mapped keys (of 101504 keys)
## org.Hs.egENSEMBLTRANS has 21903 mapped keys (of 46265 keys)
## org.Hs.egENSEMBLTRANS2EG has 164611 mapped keys (of 164611 keys)
## org.Hs.egENZYME has 2230 mapped keys (of 46265 keys)
## org.Hs.egENZYME2EG has 975 mapped keys (of 975 keys)
## org.Hs.egGENENAME has 46265 mapped keys (of 46265 keys)
## org.Hs.egGO has 18105 mapped keys (of 46265 keys)
## org.Hs.egGO2ALLEGS has 17608 mapped keys (of 17608 keys)
## org.Hs.egGO2EG has 13737 mapped keys (of 13737 keys)
## org.Hs.egMAP has 34976 mapped keys (of 46265 keys)
## org.Hs.egMAP2EG has 2498 mapped keys (of 2498 keys)
## org.Hs.egOMIM has 15648 mapped keys (of 46265 keys)
## org.Hs.egOMIM2EG has 19401 mapped keys (of 19401 keys)
## org.Hs.egPATH has 5870 mapped keys (of 46265 keys)
## org.Hs.egPATH2EG has 229 mapped keys (of 229 keys)
## org.Hs.egPFAM has 19052 mapped keys (of 46265 keys)
## org.Hs.egPMID has 31534 mapped keys (of 46265 keys)
## org.Hs.egPMID2EG has 394407 mapped keys (of 394407 keys)
## org.Hs.egPROSITE has 19052 mapped keys (of 46265 keys)
## org.Hs.egREFSEQ has 30052 mapped keys (of 46265 keys)
## org.Hs.egREFSEQ2EG has 167407 mapped keys (of 167407 keys)
## org.Hs.egSYMBOL has 46265 mapped keys (of 46265 keys)
## org.Hs.egSYMBOL2EG has 46257 mapped keys (of 46257 keys)
## org.Hs.egUCSCKG has 22914 mapped keys (of 46265 keys)
## org.Hs.egUNIGENE has 24059 mapped keys (of 46265 keys)
## org.Hs.egUNIGENE2EG has 25484 mapped keys (of 25484 keys)
## org.Hs.egUNIPROT has 19092 mapped keys (of 46265 keys)
## 
## 
## Additional Information about this package:
## 
## DB schema: HUMAN_DB
## DB schema version: 2.1
## Organism: Homo sapiens
## Date for NCBI data: 2013-Sep12
## Date for GO data: 20130907
## Date for KEGG data: 2011-Mar15
## Date for Golden Path data: 2010-Mar22
## Date for Ensembl data: 2013-Sep3
```



### code chunk number 5: list FGS (Functional Gene Set)

```r
kegg <- as.list(org.Hs.egPATH2EG)
length(kegg)
```

```
## [1] 229
```

```r
str(kegg[1:5])
```

```
## List of 5
##  $ 04610: chr [1:69] "2" "462" "623" "624" ...
##  $ 00232: chr [1:7] "9" "10" "1544" "1548" ...
##  $ 00983: chr [1:52] "9" "10" "978" "1066" ...
##  $ 01100: chr [1:1130] "9" "10" "15" "18" ...
##  $ 00380: chr [1:42] "15" "26" "38" "39" ...
```

```r
names(kegg)[1:5]
```

```
## [1] "04610" "00232" "00983" "01100" "00380"
```

```r
go <- as.list(org.Hs.egGO2ALLEGS)
length(go)
```

```
## [1] 17608
```

```r
str(go[1:5])
```

```
## List of 5
##  $ GO:0000002: Named chr [1:20] "291" "1763" "1890" "3980" ...
##   ..- attr(*, "names")= chr [1:20] "TAS" "IDA" "TAS" "IEA" ...
##  $ GO:0000003: Named chr [1:1328] "18" "49" "49" "49" ...
##   ..- attr(*, "names")= chr [1:1328] "IEA" "IEA" "IMP" "ISS" ...
##  $ GO:0000012: Named chr [1:9] "3981" "7141" "7515" "23411" ...
##   ..- attr(*, "names")= chr [1:9] "IDA" "IDA" "IEA" "IMP" ...
##  $ GO:0000018: Named chr [1:52] "604" "641" "641" "940" ...
##   ..- attr(*, "names")= chr [1:52] "IEA" "IEA" "IMP" "IEA" ...
##  $ GO:0000019: Named chr [1:4] "641" "4292" "4361" "10111"
##   ..- attr(*, "names")= chr [1:4] "IEA" "IEA" "TAS" "IDA"
```

```r
names(go)[1:5]
```

```
## [1] "GO:0000002" "GO:0000003" "GO:0000012" "GO:0000018" "GO:0000019"
```


### code chunk number 6: convertIDs

```r
numberOfFGSkegg <- 200
kegg <- lapply(kegg[sample(1:length(kegg), numberOfFGSkegg)], function(x) unique(unlist(mget(x, 
    org.Hs.egSYMBOL))))
str(kegg[1:5])
```

```
## List of 5
##  $ 05145: chr [1:132] "AKT1" "AKT2" "ALOX5" "BIRC2" ...
##  $ 05340: chr [1:35] "ADA" "AIRE" "BTK" "CD3D" ...
##  $ 00232: chr [1:7] "NAT1" "NAT2" "CYP1A2" "CYP2A6" ...
##  $ 00460: chr [1:7] "GGT1" "GGT7" "GGT5" "SHMT1" ...
##  $ 00563: chr [1:25] "GPLD1" "PIGA" "PIGC" "PIGF" ...
```

```r
### Process GO: keep only Biological Process
length(go)
```

```
## [1] 17608
```

```r
go <- go[Ontology(names(go)) == "BP"]
length(go)
```

```
## [1] 12288
```

```r
numberOfFGSgo <- 200
go <- lapply(go[sample(1:length(go), numberOfFGSgo)], function(x) unique(unlist(mget(x, 
    org.Hs.egSYMBOL))))
str(go[1:5])
```

```
## List of 5
##  $ GO:1901993: chr [1:2] "OVOL1" "STRA8"
##  $ GO:0021917: chr "PAX6"
##  $ GO:0045682: chr [1:43] "AQP3" "ATOH1" "CYP27B1" "EZH2" ...
##  $ GO:0002698: chr [1:59] "A2M" "ABR" "ARRB2" "BCL6" ...
##  $ GO:0010955: chr [1:16] "A2M" "SERPING1" "CD59" "CR1" ...
```



### code chunk number 7: annotateFGS

```r
require(KEGG.db)
KEGG()
```

```
## Quality control information for KEGG:
## 
## 
## This package has the following mappings:
## 
## KEGGENZYMEID2GO has 3999 mapped keys (of 3999 keys)
## KEGGEXTID2PATHID has 75100 mapped keys (of 75100 keys)
## KEGGGO2ENZYMEID has 4129 mapped keys (of 4129 keys)
## KEGGPATHID2EXTID has 3152 mapped keys (of 3152 keys)
## KEGGPATHID2NAME has 390 mapped keys (of 390 keys)
## KEGGPATHNAME2ID has 390 mapped keys (of 390 keys)
## 
## 
## Additional Information about this package:
## 
## DB schema: KEGG_DB
## DB schema version: 2.1
## Date for KEGG data: 2011-Mar15
```

```r
names(kegg) <- paste(names(kegg), unlist(mget(names(kegg), KEGGPATHID2NAME)), 
    sep = ".")
head(names(kegg), n = 10)
```

```
##  [1] "05145.Toxoplasmosis"                                        
##  [2] "05340.Primary immunodeficiency"                             
##  [3] "00232.Caffeine metabolism"                                  
##  [4] "00460.Cyanoamino acid metabolism"                           
##  [5] "00563.Glycosylphosphatidylinositol(GPI)-anchor biosynthesis"
##  [6] "05150.Staphylococcus aureus infection"                      
##  [7] "05200.Pathways in cancer"                                   
##  [8] "00472.D-Arginine and D-ornithine metabolism"                
##  [9] "04710.Circadian rhythm - mammal"                            
## [10] "04660.T cell receptor signaling pathway"
```


### code chunk number 8: listFGS

```r
require(GO.db)
GO()
```

```
## Quality control information for GO:
## 
## 
## This package has the following mappings:
## 
## GOBPANCESTOR has 25193 mapped keys (of 25193 keys)
## GOBPCHILDREN has 14497 mapped keys (of 25193 keys)
## GOBPOFFSPRING has 14497 mapped keys (of 25193 keys)
## GOBPPARENTS has 25193 mapped keys (of 25193 keys)
## GOCCANCESTOR has 3232 mapped keys (of 3232 keys)
## GOCCCHILDREN has 1070 mapped keys (of 3232 keys)
## GOCCOFFSPRING has 1070 mapped keys (of 3232 keys)
## GOCCPARENTS has 3232 mapped keys (of 3232 keys)
## GOMFANCESTOR has 9602 mapped keys (of 9602 keys)
## GOMFCHILDREN has 1930 mapped keys (of 9602 keys)
## GOMFOFFSPRING has 1930 mapped keys (of 9602 keys)
## GOMFPARENTS has 9602 mapped keys (of 9602 keys)
## GOOBSOLETE has 1838 mapped keys (of 1838 keys)
## GOTERM has 38028 mapped keys (of 38028 keys)
## 
## 
## Additional Information about this package:
## 
## DB schema: GO_DB
## DB schema version: 2.1
## Date for GO data: 20130907
```

```r
names(go) <- paste(names(go), Term(names(go)), sep = ".")
head(names(go), n = 10)
```

```
##  [1] "GO:1901993.regulation of meiotic cell cycle phase transition"
##  [2] "GO:0021917.somatic motor neuron fate commitment"             
##  [3] "GO:0045682.regulation of epidermis development"              
##  [4] "GO:0002698.negative regulation of immune effector process"   
##  [5] "GO:0010955.negative regulation of protein processing"        
##  [6] "GO:0050852.T cell receptor signaling pathway"                
##  [7] "GO:0048631.regulation of skeletal muscle tissue growth"      
##  [8] "GO:0001711.endodermal cell fate commitment"                  
##  [9] "GO:0072602.interleukin-4 secretion"                          
## [10] "GO:0001697.histamine-induced gastric acid secretion"
```


### code chunk number 9: listFGS

```r
fgsList <- list(go = go, kegg = kegg)
str(fgsList$go[1:5])
```

```
## List of 5
##  $ GO:1901993.regulation of meiotic cell cycle phase transition: chr [1:2] "OVOL1" "STRA8"
##  $ GO:0021917.somatic motor neuron fate commitment             : chr "PAX6"
##  $ GO:0045682.regulation of epidermis development              : chr [1:43] "AQP3" "ATOH1" "CYP27B1" "EZH2" ...
##  $ GO:0002698.negative regulation of immune effector process   : chr [1:59] "A2M" "ABR" "ARRB2" "BCL6" ...
##  $ GO:0010955.negative regulation of protein processing        : chr [1:16] "A2M" "SERPING1" "CD59" "CR1" ...
```

```r
str(fgsList$kegg[1:5])
```

```
## List of 5
##  $ 05145.Toxoplasmosis                                        : chr [1:132] "AKT1" "AKT2" "ALOX5" "BIRC2" ...
##  $ 05340.Primary immunodeficiency                             : chr [1:35] "ADA" "AIRE" "BTK" "CD3D" ...
##  $ 00232.Caffeine metabolism                                  : chr [1:7] "NAT1" "NAT2" "CYP1A2" "CYP2A6" ...
##  $ 00460.Cyanoamino acid metabolism                           : chr [1:7] "GGT1" "GGT7" "GGT5" "SHMT1" ...
##  $ 00563.Glycosylphosphatidylinositol(GPI)-anchor biosynthesis: chr [1:25] "GPLD1" "PIGA" "PIGC" "PIGF" ...
```


### code chunk number 10: runGSEbatchArgs

```r
args(runBatchGSE)
```

```
## function (dataList, fgsList, ...) 
## NULL
```



### code chunk number 11: runBatchGSE.separate

```r
gseABS.sep <- runBatchGSE(dataList = sepScores, fgsList = fgsList)
gseABS.sep <- runBatchGSE(dataList = sepScores, fgsList = fgsList, absolute = TRUE, 
    type = "f", alternative = "mixed")
```


### code chunk number 12: runBatchGSE.separate2

```r
gseUP.sep <- runBatchGSE(dataList = sepScores, fgsList = fgsList, absolute = FALSE, 
    type = "t", alternative = "up")
gseDW.sep <- runBatchGSE(dataList = sepScores, fgsList = fgsList, absolute = FALSE, 
    type = "t", alternative = "down")
gseBOTH.sep <- runBatchGSE(dataList = sepScores, fgsList = fgsList, absolute = FALSE, 
    type = "t", alternative = "either")
```


### code chunk number 13: runBatchGSE.int3

```r
gseABSsim.sep <- runBatchGSE(dataList = sepScores, fgsList = fgsList, absolute = TRUE, 
    type = "f", alternative = "mixed", ranks.only = FALSE, nsim = 1000)
gseUPsim.sep <- runBatchGSE(dataList = sepScores, fgsList = fgsList, absolute = FALSE, 
    type = "t", alternative = "up", ranks.only = FALSE, nsim = 1000)
```


### code chunk number 14: runBatchGSE.format1

```r
str(gseUP.sep[1:5])
```

```
## List of 5
##  $ dat.affy      :List of 2
##   ..$ go  : Named num [1:200] NA 0.89 0.0906 0.0347 0.6763 ...
##   .. ..- attr(*, "names")= chr [1:200] "GO:1901993.regulation of meiotic cell cycle phase transition" "GO:0021917.somatic motor neuron fate commitment" "GO:0045682.regulation of epidermis development" "GO:0002698.negative regulation of immune effector process" ...
##   ..$ kegg: Named num [1:200] 0.193 0.35 NA NA 0.101 ...
##   .. ..- attr(*, "names")= chr [1:200] "05145.Toxoplasmosis" "05340.Primary immunodeficiency" "00232.Caffeine metabolism" "00460.Cyanoamino acid metabolism" ...
##  $ dat.agilent   :List of 2
##   ..$ go  : Named num [1:200] NA 0.933 0.128 0.282 0.842 ...
##   .. ..- attr(*, "names")= chr [1:200] "GO:1901993.regulation of meiotic cell cycle phase transition" "GO:0021917.somatic motor neuron fate commitment" "GO:0045682.regulation of epidermis development" "GO:0002698.negative regulation of immune effector process" ...
##   ..$ kegg: Named num [1:200] 0.0227 0.3399 NA NA 0.0701 ...
##   .. ..- attr(*, "names")= chr [1:200] "05145.Toxoplasmosis" "05340.Primary immunodeficiency" "00232.Caffeine metabolism" "00460.Cyanoamino acid metabolism" ...
##  $ dat.cnvHarvard:List of 2
##   ..$ go  : Named num [1:200] NA 0.751 0.708 0.476 0.29 ...
##   .. ..- attr(*, "names")= chr [1:200] "GO:1901993.regulation of meiotic cell cycle phase transition" "GO:0021917.somatic motor neuron fate commitment" "GO:0045682.regulation of epidermis development" "GO:0002698.negative regulation of immune effector process" ...
##   ..$ kegg: Named num [1:200] 0.493 0.821 NA NA 0.62 ...
##   .. ..- attr(*, "names")= chr [1:200] "05145.Toxoplasmosis" "05340.Primary immunodeficiency" "00232.Caffeine metabolism" "00460.Cyanoamino acid metabolism" ...
##  $ dat.cnvMskcc  :List of 2
##   ..$ go  : Named num [1:200] NA 0.257 0.671 0.121 0.938 ...
##   .. ..- attr(*, "names")= chr [1:200] "GO:1901993.regulation of meiotic cell cycle phase transition" "GO:0021917.somatic motor neuron fate commitment" "GO:0045682.regulation of epidermis development" "GO:0002698.negative regulation of immune effector process" ...
##   ..$ kegg: Named num [1:200] 0.722 0.688 NA NA 0.673 ...
##   .. ..- attr(*, "names")= chr [1:200] "05145.Toxoplasmosis" "05340.Primary immunodeficiency" "00232.Caffeine metabolism" "00460.Cyanoamino acid metabolism" ...
##  $ NA            : NULL
```

```r
head(gseABSsim.sep$dat.affy$go)
```

```
## GO:1901993.regulation of meiotic cell cycle phase transition 
##                                                           NA 
##              GO:0021917.somatic motor neuron fate commitment 
##                                                       0.1688 
##               GO:0045682.regulation of epidermis development 
##                                                       0.9770 
##    GO:0002698.negative regulation of immune effector process 
##                                                       0.9940 
##         GO:0010955.negative regulation of protein processing 
##                                                       0.3556 
##                 GO:0050852.T cell receptor signaling pathway 
##                                                       0.1169
```

```r
head(gseABSsim.sep$dat.affy$kegg)
```

```
##                                         05145.Toxoplasmosis 
##                                                      0.6484 
##                              05340.Primary immunodeficiency 
##                                                      0.4745 
##                                   00232.Caffeine metabolism 
##                                                          NA 
##                            00460.Cyanoamino acid metabolism 
##                                                          NA 
## 00563.Glycosylphosphatidylinositol(GPI)-anchor biosynthesis 
##                                                      0.8901 
##                       05150.Staphylococcus aureus infection 
##                                                          NA
```


### code chunk number 15: runBatchGSE.altFunc

```r
require(limma)
gseUP.sep.2 <- runBatchGSE(dataList = sepScores, fgsList = fgsList, absolute = FALSE, 
    gseFunc = wilcoxGST, alternative = "up")
```


### code chunk number 16: runBatchGSE.format2

```r
str(gseUP.sep.2)
```

```
## List of 4
##  $ dat.affy      :List of 2
##   ..$ go  : Named num [1:200] NA 0.89 0.0906 0.0347 0.6763 ...
##   .. ..- attr(*, "names")= chr [1:200] "GO:1901993.regulation of meiotic cell cycle phase transition" "GO:0021917.somatic motor neuron fate commitment" "GO:0045682.regulation of epidermis development" "GO:0002698.negative regulation of immune effector process" ...
##   ..$ kegg: Named num [1:200] 0.193 0.35 NA NA 0.101 ...
##   .. ..- attr(*, "names")= chr [1:200] "05145.Toxoplasmosis" "05340.Primary immunodeficiency" "00232.Caffeine metabolism" "00460.Cyanoamino acid metabolism" ...
##  $ dat.agilent   :List of 2
##   ..$ go  : Named num [1:200] NA 0.933 0.128 0.282 0.842 ...
##   .. ..- attr(*, "names")= chr [1:200] "GO:1901993.regulation of meiotic cell cycle phase transition" "GO:0021917.somatic motor neuron fate commitment" "GO:0045682.regulation of epidermis development" "GO:0002698.negative regulation of immune effector process" ...
##   ..$ kegg: Named num [1:200] 0.0227 0.3399 NA NA 0.0701 ...
##   .. ..- attr(*, "names")= chr [1:200] "05145.Toxoplasmosis" "05340.Primary immunodeficiency" "00232.Caffeine metabolism" "00460.Cyanoamino acid metabolism" ...
##  $ dat.cnvHarvard:List of 2
##   ..$ go  : Named num [1:200] NA 0.751 0.708 0.476 0.29 ...
##   .. ..- attr(*, "names")= chr [1:200] "GO:1901993.regulation of meiotic cell cycle phase transition" "GO:0021917.somatic motor neuron fate commitment" "GO:0045682.regulation of epidermis development" "GO:0002698.negative regulation of immune effector process" ...
##   ..$ kegg: Named num [1:200] 0.493 0.821 NA NA 0.62 ...
##   .. ..- attr(*, "names")= chr [1:200] "05145.Toxoplasmosis" "05340.Primary immunodeficiency" "00232.Caffeine metabolism" "00460.Cyanoamino acid metabolism" ...
##  $ dat.cnvMskcc  :List of 2
##   ..$ go  : Named num [1:200] NA 0.257 0.671 0.121 0.938 ...
##   .. ..- attr(*, "names")= chr [1:200] "GO:1901993.regulation of meiotic cell cycle phase transition" "GO:0021917.somatic motor neuron fate commitment" "GO:0045682.regulation of epidermis development" "GO:0002698.negative regulation of immune effector process" ...
##   ..$ kegg: Named num [1:200] 0.722 0.688 NA NA 0.673 ...
##   .. ..- attr(*, "names")= chr [1:200] "05145.Toxoplasmosis" "05340.Primary immunodeficiency" "00232.Caffeine metabolism" "00460.Cyanoamino acid metabolism" ...
```

```r
all(gseUP.sep.2$go == gseUP.sep$go)
```

```
## [1] TRUE
```


### code chunk number 17: runBatchGSE.altFunc2

```r
gseFunc <- function(selected, statistics, threshold) {
    diffExpGenes <- statistics > threshold
    tab <- table(diffExpGenes, selected)
    pVal <- fisher.test(tab)[["p.value"]]
}
gseUP.sep.3 <- runBatchGSE(dataList = sepScores, fgsList = fgsList, absolute = FALSE, 
    gseFunc = gseFunc, threshold = 7.5)
```


### code chunk number 18: runBatchGSE.format3

```r
str(gseUP.sep.3)
```

```
## List of 4
##  $ dat.affy      :List of 2
##   ..$ go  : Named num [1:200] NA 1 1 1 1 1 NA NA NA NA ...
##   .. ..- attr(*, "names")= chr [1:200] "GO:1901993.regulation of meiotic cell cycle phase transition" "GO:0021917.somatic motor neuron fate commitment" "GO:0045682.regulation of epidermis development" "GO:0002698.negative regulation of immune effector process" ...
##   ..$ kegg: Named num [1:200] 1 1 NA NA 1 NA 1 NA 1 1 ...
##   .. ..- attr(*, "names")= chr [1:200] "05145.Toxoplasmosis" "05340.Primary immunodeficiency" "00232.Caffeine metabolism" "00460.Cyanoamino acid metabolism" ...
##  $ dat.agilent   :List of 2
##   ..$ go  : Named num [1:200] NA 1 1 1 1 1 NA NA NA NA ...
##   .. ..- attr(*, "names")= chr [1:200] "GO:1901993.regulation of meiotic cell cycle phase transition" "GO:0021917.somatic motor neuron fate commitment" "GO:0045682.regulation of epidermis development" "GO:0002698.negative regulation of immune effector process" ...
##   ..$ kegg: Named num [1:200] 1 1 NA NA 1 NA 1 NA 1 1 ...
##   .. ..- attr(*, "names")= chr [1:200] "05145.Toxoplasmosis" "05340.Primary immunodeficiency" "00232.Caffeine metabolism" "00460.Cyanoamino acid metabolism" ...
##  $ dat.cnvHarvard:List of 2
##   ..$ go  : Named num [1:200] NA 1 1 1 1 1 NA NA NA NA ...
##   .. ..- attr(*, "names")= chr [1:200] "GO:1901993.regulation of meiotic cell cycle phase transition" "GO:0021917.somatic motor neuron fate commitment" "GO:0045682.regulation of epidermis development" "GO:0002698.negative regulation of immune effector process" ...
##   ..$ kegg: Named num [1:200] 1 1 NA NA 1 NA 1 NA 1 1 ...
##   .. ..- attr(*, "names")= chr [1:200] "05145.Toxoplasmosis" "05340.Primary immunodeficiency" "00232.Caffeine metabolism" "00460.Cyanoamino acid metabolism" ...
##  $ dat.cnvMskcc  :List of 2
##   ..$ go  : Named num [1:200] NA 1 1 1 1 1 NA NA NA NA ...
##   .. ..- attr(*, "names")= chr [1:200] "GO:1901993.regulation of meiotic cell cycle phase transition" "GO:0021917.somatic motor neuron fate commitment" "GO:0045682.regulation of epidermis development" "GO:0002698.negative regulation of immune effector process" ...
##   ..$ kegg: Named num [1:200] 1 1 NA NA 1 NA 1 NA 1 1 ...
##   .. ..- attr(*, "names")= chr [1:200] "05145.Toxoplasmosis" "05340.Primary immunodeficiency" "00232.Caffeine metabolism" "00460.Cyanoamino acid metabolism" ...
```

```r
head(data.frame(fisher = gseUP.sep.3$dat.affy$kegg, wilcoxon = gseUP.sep$dat.affy$kegg))
```

```
##                                                             fisher
## 05145.Toxoplasmosis                                              1
## 05340.Primary immunodeficiency                                   1
## 00232.Caffeine metabolism                                       NA
## 00460.Cyanoamino acid metabolism                                NA
## 00563.Glycosylphosphatidylinositol(GPI)-anchor biosynthesis      1
## 05150.Staphylococcus aureus infection                           NA
##                                                             wilcoxon
## 05145.Toxoplasmosis                                           0.1927
## 05340.Primary immunodeficiency                                0.3502
## 00232.Caffeine metabolism                                         NA
## 00460.Cyanoamino acid metabolism                                  NA
## 00563.Glycosylphosphatidylinositol(GPI)-anchor biosynthesis   0.1012
## 05150.Staphylococcus aureus infection                             NA
```


### code chunk number 19: adjustP

```r
gseABS.sep.BH <- adjustPvalGSE(gseABS.sep)
gseABS.sep.holm <- adjustPvalGSE(gseABS.sep, proc = "Holm")
```


### code chunk number 20: adjusted.format

```r
names(gseABS.sep.BH)
```

```
## [1] "dat.affy"       "dat.agilent"    "dat.cnvHarvard" "dat.cnvMskcc"
```

```r
names(gseABS.sep.holm)
```

```
## [1] "dat.affy"       "dat.agilent"    "dat.cnvHarvard" "dat.cnvMskcc"
```

```r
str(gseABS.sep.BH$dat.affy)
```

```
## List of 2
##  $ go  : num [1:200, 1:2] NA 0.124 0.952 0.987 0.349 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : chr [1:200] "GO:1901993.regulation of meiotic cell cycle phase transition" "GO:0021917.somatic motor neuron fate commitment" "GO:0045682.regulation of epidermis development" "GO:0002698.negative regulation of immune effector process" ...
##   .. ..$ : chr [1:2] "rawp" "BH"
##  $ kegg: num [1:200, 1:2] 0.855 0.702 NA NA 0.918 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : chr [1:200] "05145.Toxoplasmosis" "05340.Primary immunodeficiency" "00232.Caffeine metabolism" "00460.Cyanoamino acid metabolism" ...
##   .. ..$ : chr [1:2] "rawp" "BH"
```

```r
str(gseABS.sep.holm$dat.affy)
```

```
## List of 2
##  $ go  : num [1:200, 1:2] NA 0.124 0.952 0.987 0.349 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : chr [1:200] "GO:1901993.regulation of meiotic cell cycle phase transition" "GO:0021917.somatic motor neuron fate commitment" "GO:0045682.regulation of epidermis development" "GO:0002698.negative regulation of immune effector process" ...
##   .. ..$ : chr [1:2] "rawp" "Holm"
##  $ kegg: num [1:200, 1:2] 0.855 0.702 NA NA 0.918 ...
##   ..- attr(*, "dimnames")=List of 2
##   .. ..$ : chr [1:200] "05145.Toxoplasmosis" "05340.Primary immunodeficiency" "00232.Caffeine metabolism" "00460.Cyanoamino acid metabolism" ...
##   .. ..$ : chr [1:2] "rawp" "Holm"
```

```r
head(gseABS.sep.BH$dat.affy$go, n = 10)
```

```
##                                                                rawp BH
## GO:1901993.regulation of meiotic cell cycle phase transition     NA NA
## GO:0021917.somatic motor neuron fate commitment              0.1236  1
## GO:0045682.regulation of epidermis development               0.9518  1
## GO:0002698.negative regulation of immune effector process    0.9873  1
## GO:0010955.negative regulation of protein processing         0.3490  1
## GO:0050852.T cell receptor signaling pathway                 0.1350  1
## GO:0048631.regulation of skeletal muscle tissue growth           NA NA
## GO:0001711.endodermal cell fate commitment                       NA NA
## GO:0072602.interleukin-4 secretion                               NA NA
## GO:0001697.histamine-induced gastric acid secretion              NA NA
```

```r
head(gseABS.sep.holm$dat.affy$go, n = 10)
```

```
##                                                                rawp Holm
## GO:1901993.regulation of meiotic cell cycle phase transition     NA   NA
## GO:0021917.somatic motor neuron fate commitment              0.1236    1
## GO:0045682.regulation of epidermis development               0.9518    1
## GO:0002698.negative regulation of immune effector process    0.9873    1
## GO:0010955.negative regulation of protein processing         0.3490    1
## GO:0050852.T cell receptor signaling pathway                 0.1350    1
## GO:0048631.regulation of skeletal muscle tissue growth           NA   NA
## GO:0001711.endodermal cell fate commitment                       NA   NA
## GO:0072602.interleukin-4 secretion                               NA   NA
## GO:0001697.histamine-induced gastric acid secretion              NA   NA
```


### code chunk number 21: figure004

```r
hist(gseABS.sep.BH$dat.affy$go[, "rawp"], nclass = 10, col = "orange", main = "P-values")
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15.png) 


### code chunk number 22: figure005

```r
hist(gseABS.sep.BH$dat.affy$go[, "BH"], nclass = 10, col = "orange", main = "Q-values")
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16.png) 



### code chunk number 23: readAndProcessData

```r
autismGenes <- read.table("myData/101symbols_v02.txt", sep = "\t", header = TRUE, 
    colClasses = "character")
```

```
## Warning: cannot open file 'myData/101symbols_v02.txt': No such file or
## directory
```

```
## Error: cannot open the connection
```

```r
str(autismGenes)
```

```
## Error: object 'autismGenes' not found
```


### code chunk number 24: newFunctionForLists

```r
formatGeneList <- function(geneList, allGenes) {
    stats <- 1 * allGenes %in% geneList
    names(stats) <- allGenes
    out <- list(GeneList = stats)
}
```


### code chunk number 25: formatTheGeneLsit

```r
allGenes <- unique(unlist(as.list(org.Hs.egSYMBOL)))
str(allGenes)
```

```
##  chr [1:46257] "A1BG" "NAT2" "ADA" "CDH2" "AKT3" "GAGE12F" ...
```

```r
geneMembership <- formatGeneList(autismGenes$SYMBOL, allGenes)
```

```
## Error: object 'autismGenes' not found
```

```r
str(geneMembership)
```

```
## Error: object 'geneMembership' not found
```

```r
table(geneMembership$GeneList)
```

```
## Error: object 'geneMembership' not found
```


### code chunk number 26: runBatchGSE.altFunc3

```r
gseMyList <- runBatchGSE(dataList = geneMembership, fgsList = fgsList, absolute = FALSE, 
    gseFunc = gseFunc, threshold = 0.5)
```

```
## Error: object 'geneMembership' not found
```


### code chunk number 27: showTheResults

```r
str(gseMyList)
```

```
## Error: object 'gseMyList' not found
```

```r
str(gseMyList$GeneList$kegg)
```

```
## Error: object 'gseMyList' not found
```

```r
summary(gseMyList$GeneList$kegg)
```

```
## Error: error in evaluating the argument 'object' in selecting a method for function 'summary': Error: object 'gseMyList' not found
```

```r
table(is.na(gseMyList$GeneList$kegg))
```

```
## Error: object 'gseMyList' not found
```

```r
table(gseMyList$GeneList$kegg < 0.01)
```

```
## Error: object 'gseMyList' not found
```

```r
table(gseMyList$GeneList$go < 0.01)
```

```
## Error: object 'gseMyList' not found
```


### code chunk number 28: adjustP.2

```r
gseMyList.BH <- adjustPvalGSE(gseMyList)
```

```
## Error: object 'gseMyList' not found
```

```r
gseMyList.BH$GeneList$go[gseMyList.BH$GeneList$go[, "BH"] < 0.01, ]
```

```
## Error: object 'gseMyList.BH' not found
```

```r
gseMyList.BH$GeneList$kegg[gseMyList.BH$GeneList$kegg[, "BH"] < 0.01, ]
```

```
## Error: object 'gseMyList.BH' not found
```


### code chunk number 29: sessioInfo

```r
sessionInfo()
```

```
## R version 3.0.3 (2014-03-06)
## Platform: x86_64-apple-darwin13.1.0 (64-bit)
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] GO.db_2.10.1         KEGG.db_2.10.1       multtest_2.18.0     
##  [4] annotate_1.40.1      org.Hs.eg.db_2.10.1  RSQLite_0.11.4      
##  [7] DBI_0.2-7            AnnotationDbi_1.24.0 RTopper_1.8.0       
## [10] limma_3.18.13        Biobase_2.22.0       BiocGenerics_0.8.0  
## [13] BiocInstaller_1.12.1 knitr_1.5           
## 
## loaded via a namespace (and not attached):
##  [1] evaluate_0.5.5  formatR_0.10    IRanges_1.20.7  MASS_7.3-33    
##  [5] splines_3.0.3   stats4_3.0.3    stringr_0.6.2   survival_2.37-7
##  [9] tools_3.0.3     XML_3.98-1.1    xtable_1.7-3
```


