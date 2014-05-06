### make the mean over the levels of a factor
### newDat is a data.frame, we want to make the mean by column 
### the unique ids will be levels of the factor
tmp <- sapply(newDat, FUN=function(x, y) {
	if(is.numeric(x)){
		tapply(X=x, INDEX=y, FUN=mean )
	} else {
		unique(x)
	}
}, y=newDat$ids)
