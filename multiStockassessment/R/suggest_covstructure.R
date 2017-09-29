##' Construct valid (band) correlation structures for multisam.fit
##'
##' @title Suggest correlation structure
##' @param x a samset
##' @param nAgeClose Number of close ages to correlate
##' @param noConnection If TRUE, no correlations
##' @param onlyCloseAge If TRUE, age classes no less than nAgeClose are un-correlated
##' @param noCorInArea If TRUE, age classes in same area/stock are un-correlated
##' @param onlyCorInArea If TRUE, only age classes in same area/stock are correlated
##' @param noCorBetween matrix of area/stock numbers that should not be correlated
##' @return A boolean matrix
##' @author Christoffer Moesgaard Albertsen
##' @export
suggestCorStructure <- function(x,
                                nAgeClose = 1,
                                noConnection = FALSE,
                                onlyCloseAge = TRUE,
                                noCorInArea = TRUE,
                                onlyCorInArea = FALSE,
                                noCorBetween = matrix(0,0,2)
                                ){
    nAreas <- length(x)
    dat <- collect_data(x)    
    ageList <- lapply(dat,function(x)x$minAge:x$maxAge)
    minAgeAll <- min(unlist(ageList))
    maxAgeAll <- max(unlist(ageList))

    
    cons <- list()
    ala <- minAgeAll:maxAgeAll
    nn <- nAreas * length(ala)

    mm <- matrix(FALSE,nn,nn)
    indx <- which(unlist(sapply(1:length(ageList),function(i){
        !(ala %in% ageList[[i]]) ##(ageList[[i]] - minAgeAll))) + (i-1)*length(ala)
    })))

    mm[indx,] <- TRUE
    mm[,indx] <- TRUE
    if(noConnection)
        mm[,] <- TRUE

    if(onlyCloseAge)
        mm[which(abs(matrix(ala,nn,nn) - matrix(ala,nn,nn,byrow=TRUE)) >= nAgeClose)] <- TRUE

    if(noCorInArea){
        for(i in 1:nAreas)
            mm[1:length(ala) + max(length(ala))*(i-1),1:length(ala) + max(length(ala))*(i-1)] <- TRUE
    }

    if(onlyCorInArea){
        for(i in 1:(nAreas-1))
            for(j in (i+1):nAreas)
                mm[1:length(ala) + max(length(ala))*(j-1),1:length(ala) + max(length(ala))*(i-1)] <- mm[1:length(ala) + max(length(ala))*(i-1),1:length(ala) + max(length(ala))*(j-1)] <- TRUE
    }


    if(nrow(noCorBetween) > 0){
        for(k in 1:nrow(noCorBetween)){
            i <- which(noCorBetween[k,1] == useAreaNum)
            j <- which(noCorBetween[k,2] == useAreaNum)
            if(length(i) > 0 & length(j) > 0){
                indx1 <- 1:length(ala) + max(length(ala))*(i-1)
                indx2 <- 1:length(ala) + max(length(ala))*(j-1)
                mm[indx1,indx2] <- TRUE
                mm[indx2,indx1] <- TRUE
            }
        }
    }

    
    return(mm)
}
