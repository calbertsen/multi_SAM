collect_pars <- function(x) {
    if(class(x) != "samset")
        stop("x must be a samset.")
    
    newparam <- list()
    ncase <- length(x)
    nparam <- length(x[[1]]$pl)
    for(i in 1:nparam){ # loop over different parameters
        if(is.vector(x[[1]]$pl[[i]])){            
            newparam[[i]] <- combineVectors(lapply(x,function(yy) yy$pl[[i]]))
        }else if(is.matrix(x[[1]]$pl[[i]])){
            newparam[[i]] <- combineMatrices(lapply(x,function(yy) yy$pl[[i]]))
        }else{
            stop("Type not implemented")
        }
    }
    names(newparam) <- names(x[[1]]$pl)
    return(newparam)
}
