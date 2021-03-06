Basic R and Other Resources
========================================================

Goal: To introduce a few useful commands and tools in R.
--------------------------------------------------------

**We recommend you use RStudio for this class, on all operating systems. You can use the following link** [http://www.rstudio.com/ide/download/](http://www.rstudio.com/ide/download/) **or good old fashioned google.**

If you want to extract the R code only from a markdown file, use the following command in the library **knitr** called *purl*.



```r
library(knitr)
knitr::purl("01_raw2dge.Rmd")
```



Useful resources for the class and life after
-------------------------------

1. The affy vignette is a very good source, and is highly recommended. You can get this by installing and then just running **vignette("affy")** .

2. The very comprehensive tutorial section by UC-Riverside [http://manuals.bioinformatics.ucr.edu/](http://manuals.bioinformatics.ucr.edu/) .

3. R-Markdown is a very good tool to document your code, and figure out exactly where you are going wrong. The R Package **knitr**, does a fantastic job and is very well documented. All these class materials have been produced using **knitr**.

4. **limma** vignette which is available on Bioconductor. [http://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf](http://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf)

5. More information on how to make robust design matrices for experiments can be found in, **Statistics for Biology and Health - Springer**, **Chapter V, Page 399-417**.[http://www.statsci.org/smyth/pubs/limma-biocbook-reprint.pdf](http://www.statsci.org/smyth/pubs/limma-biocbook-reprint.pdf)

6. RTopper vignette, [http://www.bioconductor.org/packages/release/bioc/vignettes/RTopper/inst/doc/RTopper.pdf](http://www.bioconductor.org/packages/release/bioc/vignettes/RTopper/inst/doc/RTopper.pdf).
