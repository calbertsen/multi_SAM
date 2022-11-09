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
split3DArrays <- function(x){
    cdim <- attr(x,"cdim")
    rdim <- attr(x,"rdim")
    adim <- attr(x,"adim")
    res <- unname(split(x,factor(rep(seq_along(cdim),times = rdim * cdim * adim), levels = seq_along(cdim))))
    sapply(1:length(res), function(i) array(res[[i]], c(rdim[i], cdim[i], adim[i])), simplify = FALSE)
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

combine3DArrays <- function(x){
    rdim <- sapply(x,function(a) dim(a)[1])
    cdim <- sapply(x,function(a) dim(a)[2])
    adim <- sapply(x,function(a) dim(a)[3])
    vec <- as.array(unlist(x))
    attr(vec,"rdim") <- as.integer(rdim)
    attr(vec,"cdim") <- as.integer(cdim)
    attr(vec,"adim") <- as.integer(adim)
    vec    
}

splitParameter <- function(x){
    if(!is.null(attr(x,"vdim"))){
        return(splitVectors(x))
    }else if(!is.null(attr(x,"cdim")) && !is.null(attr(x,"rdim"))){
        if(is.null(attr(x,"adim"))){
            return(splitMatrices(x))
        }else{
            return(split3DArrays(x))
        }
    }
    stop("The object to split must have vdim or cdim and rdim attributes.")
}

combineParameter <- function(x){
    if(is.matrix(x[[1]]))
        return(combineMatrices(x))
    if(is.array(x[[1]]) && length(dim(x[[1]]))==3)
        return(combine3DArrays(x))
    return(combineVectors(x))
}
