##' Simulate from a msam object
##' @param object msam object result from multisam.fit
##' @param nsim Number of simulations
##' @param seed random number seed
##' @param full.data should a full data set for sam.fit be returned?
##' @param ... Other arguments not used
##' @return a list of lists.
##' @importFrom stats simulate
##' @method simulate msam
##' @inheritParams stats::logLik
##' @author Christoffer Moesgaard Albertsen
##' @export
simulate.msam <- function(object,
                          nsim = 1,
                          seed = NULL,
                          full.data = TRUE,
                          ready.to.fit = FALSE,
                          ...){
    if(!is.null(seed)){
        set.seed(seed)
    }else{
        if(!exists(".Random.seed"))
            set.seed(NULL)
    }
    rngSeed <- .Random.seed
    obj <- attr(object,"m_obj")
    par <- obj$env$last.par.best
    dots <- list(...)
    sn <- getStockNames(object)
    if(full.data){
        ret <- replicate(nsim, {
            allDat <- lapply(object,function(x)x$data)
            simObs <- obj$simulate(par)
            for(i in 1:length(allDat)){
                allDat[[i]]$logobs <- simObs$logobs[[i]]
            }
            if(ready.to.fit){
                fa <- formalArgs(stockassessment::sam.fit)
                fa <- fa[fa != "..."]
                res <- vector("list",length(object))
                for(i in 1:length(res)){
                    res[[i]] <- list(data = allDat[[i]],
                                     conf = object[[i]]$conf,
                                     parameters = stockassessment::defpar(allDat[[i]], object[[i]]$conf))
                    indx <- which(names(dots) %in% fa)
                    if(length(indx) > 0)
                        res[[i]] <- c(res[[i]],dots[indx])
                }
                attr(res,"sim_logN") <- simObs$logN
                attr(res,"sim_logF") <- simObs$logF
                names(res) <- sn
                return(res)
            }else{
                attr(allDat,"sim_logN") <- simObs$logN
                attr(allDat,"sim_logF") <- simObs$logF
                names(allDat) <- sn
                return(allDat)
            }
        }, simplify=FALSE)
    }else{
        if(ready.to.fit)
            warning("ready.to.fit is ignored when full.data is false")
  	ret <- replicate(nsim, obj$simulate(par), simplify=FALSE)
    }
    attr(ret,".Random.seed") <- rngSeed
    attr(ret,"m_fit") <- object
    ret
}
