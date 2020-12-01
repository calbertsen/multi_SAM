splitVectors <- function(x){
    vdim <- attr(x,"vdim")
    unname(split(x,factor(rep(seq_along(vdim),times = vdim), levels = seq_along(vdim))))
}

splitMatrices <- function(x){
    cdim <- attr(x,"cdim")
    rdim <- attr(x,"rdim")
    res <- unname(split(x,factor(rep(seq_along(cdim),times = rdim * cdim), levels = seq_along(cdim))))
    sapply(1:length(res), function(i) matrix(res[[i]], rdim[i], cdim[i]), simplify = FALSE)
}

combineVectors <- function(x){
    vdim <- sapply(x,length)
    vec <- as.array(unlist(x))
    attr(vec,"vdim") <- as.integer(vdim)
    vec
}

combineMatrices <- function(x){
    cdim <- sapply(x,ncol)
    rdim <- sapply(x,nrow)
    vec <- as.array(unlist(x))
    attr(vec,"cdim") <- as.integer(cdim)
    attr(vec,"rdim") <- as.integer(rdim)
    vec    
}

splitParameter <- function(x){
    if(!is.null(attr(x,"vdim"))){
        return(splitVectors(x))
    }else if(!is.null(attr(x,"cdim")) && !is.null(attr(x,"rdim"))){
        return(splitMatrices(x))
    }
    stop("The object to split must have vdim or cdim and rdim attributes.")
}

combineParameter <- function(x){
    if(is.matrix(x[[1]]))
        return(combineMatrices(x))
    return(combineVectors(x))
}
