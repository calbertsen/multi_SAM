collect_pars <- function(x) {
    if(class(x) != "samset")
        stop("x must be a samset.")
    
    newparam <- list()
    ncase <- length(x)
    nparam <- length(x[[1]]$pl)
    for(i in 1:nparam){ # loop over different parameters
        if(is.vector(x[[1]]$pl[[i]])){
            vec <- numeric()
            vdim <- numeric(ncase)
            for(j in 1:ncase){
                vdim[j] <- length(x[[j]]$pl[[i]])
                vec <- c(vec,x[[j]]$pl[[i]])
            }
            vec <- as.array(vec)
            attr(vec,"vdim") <- vdim
            newparam[[i]] <- vec

        }else if(is.matrix(x[[1]]$pl[[i]])){
            vec <- numeric()
            cdim <- numeric(ncase)
            rdim <- numeric(ncase)
            for(j in 1:ncase){
                cdim[j] <- ncol(x[[j]]$pl[[i]])
                rdim[j] <- nrow(x[[j]]$pl[[i]])
                vec <- c(vec,as.vector(x[[j]]$pl[[i]]))
            }
            vec <- as.array(vec)
            attr(vec,"cdim") <- cdim
            attr(vec,"rdim") <- rdim
            newparam[[i]] <- vec
        }else{
            stop("Type not implemented")
        }
    }
    names(newparam) <- names(x[[1]]$pl)
    return(newparam)
}
