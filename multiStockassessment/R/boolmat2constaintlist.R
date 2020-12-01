

boolMat2ConstraintList <- function(mm){
    mvals <- 1:sum(lower.tri(mm))
    mvals[which(mm[lower.tri(mm)])] <- NA
    indx <- which(mm,TRUE)
    indx <- indx[indx[,1] > indx[,2],,drop=FALSE]
    cons <- list()
    for(i in unique(indx[,1]))
        cons[[length(cons)+1]] <- as.numeric(c(i,sort(indx[indx[,1]==i,2],TRUE))-1)
    return(cons)
}

