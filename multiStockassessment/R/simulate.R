##' Simulate from a msam object
##' @param object msam object result from multisam.fit
##' @param nsim Number of simulations
##' @param seed random number seed
##' @param full.data should a full data set for sam.fit be returned?
##' @param ... Other arguments not used
##' @return a list of lists.
##' @author Christoffer Moesgaard Albertsen
simulate.msam <- function(object,
                          nsim = 1,
                          seed = NULL,
                          full.data = TRUE,
                          ...){
    if(!is.null(seed)){
        set.seed(seed)
    }else{
        if(!exists(.Random.seed))
            set.seed(NULL)
    }
    rngSeed <- .Random.seed
    obj <- attr(object,"m_obj")
    par <- unlist(attr(object,"m_pl"))
    if(full.data){
        ret <- replicate(nsim, {
            allDat <- lapply(object,function(x)x$data)
            simObs <- attr(object,"m_obj")$simulate(par)
            for(i in 1:length(allDat)){
                allDat[[i]]$logobs <- simObs[[i]]$logobs
            }
            attr(allDat,"sim_logN") <- simObs$logN
            attr(allDat,"sim_logF") <- simObs$logF
        }, simplify=FALSE)
    }else{
  	ret <- replicate(nsim, attr(object,"m_obj")$simulate(par), simplify=FALSE)
    }
    attr(ret,".Random.seed") <- rngSeed
    ret
}
