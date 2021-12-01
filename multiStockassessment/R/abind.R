
asplit <- function(x, along = length(dim(x)), useDims = dim(x)[-along]){
    r <- lapply(split(x, slice.index(x, along)), function(x) array(x, dim = useDims, dimnames = dimnames(x)[-along]))
    names(r) <- dimnames(x)[[along]]
    r
}

abind <- function(..., along){
    v <- list(...)
    ## check
    allDims <- lapply(v, dim)
    nDims <- unique(unlist(lapply(allDims, length)))
    if(length(nDims) > 1)
        stop("Number of dimensions does not match")  
    if(missing(along))
        along <- nDims
    useDims <- unique(do.call("rbind",allDims)[,-along,drop=FALSE])[1,]
    if(sum(!duplicated((do.call("rbind",allDims)[,-along,drop=FALSE]))) > 1)
        warning("Dimensions does not match")

    r <- do.call("c",lapply(v, function(x) asplit(x, along = along, useDims = useDims)))
    if(all(useDims == 1)){
        finalDims <- integer(length(useDims)+1)
        finalDims[along] <- length(r)
        finalDims[-along] <- useDims
        return(array(unlist(r), finalDims))
    }
    simplify2array(r)    
}
