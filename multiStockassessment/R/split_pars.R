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
