splitVectors <- function(x){
    vdim <- attr(x,"vdim")
    unname(split(x,factor(rep(1:length(vdim),times = vdim))))
}

splitMatrices <- function(x){
    cdim <- attr(x,"cdim")
    rdim <- attr(x,"rdim")
    res <- unname(split(x,factor(rep(1:length(cdim),times = rdim * cdim))))
    sapply(1:length(res), function(i) matrix(res[[i]], rdim[i], cdim[i]), simplify = FALSE)
}

combineVectors <- function(x){
    vdim <- sapply(x,length)
    vec <- as.array(unlist(x))
    attr(vec,"vdim") <- vdim
    vec
}

combineMatrices <- function(x){
    cdim <- sapply(x,ncol)
    rdim <- sapply(x,nrow)
    vec <- as.array(unlist(x))
    attr(vec,"cdim") <- cdim
    attr(vec,"rdim") <- rdim
    vec    
}

