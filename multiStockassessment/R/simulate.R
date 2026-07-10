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
                          simFlag,
                          sim.keepRec,
                          pl,
                          ...){
    if(!is.null(seed)){
        set.seed(seed)
    }else{        
        if(!exists(".Random.seed"))
            set.seed(NULL)
    }
    rngSeed <- .Random.seed
    obj <- attr(object,"m_obj")
    shObs <- attr(object,"m_dat")$sharedObs
    if(!missing(pl)){
        plMap <- attr(object,"m_pl")
        ## Should probably check names and dimensions
        plMap[names(pl)] <- pl
        map <- obj$env$map
        with.map <- intersect(names(plMap), names(map))
        applyMap <- function(par.name) {
            tapply(plMap[[par.name]], map[[par.name]], mean)
        }
        plMap[with.map] <- sapply(with.map, applyMap, simplify = FALSE)         
        par <- unlist(plMap)
    }else{
        par <- obj$env$last.par.best
    }
    dots <- list(...)
    sn <- getStockNames(object)
    if(!missing(simFlag)){
        for(i in seq_along(obj$env$data$sam))
            obj$env$data$sam[[i]]$simFlag <- rep(pmin(1,pmax(as.integer(simFlag),0)),length.out=4)
    }
    if(!missing(sim.keepRec)){
        for(i in seq_along(obj$env$data$sam))
            obj$env$data$sam[[i]]$simKeepRec <- as.numeric(sim.keepRec)
    }
    if(!missing(simFlag) || !missing(sim.keepRec))
        obj$env$retape(FALSE)
    if(full.data){
        ret <- replicate(nsim, {
            allDat <- lapply(object,function(x)x$data)
            simObs <- obj$simulate(par)
            set.seed(NULL)
            for(i in 1:length(allDat)){
                allDat[[i]]$logobs <- simObs$logobs[[i]]
            }
            shObs$logobs <- simObs$shared_logobs
            shObs$predPerStock <- simObs$shared_predPerStock
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
                bionm <- names(simObs)[grepl("^(bio|total)_",names(simObs))]
                for(nn in bionm)
                    attr(res,nn) <- simObs[[nn]]
                attr(res,"shared_data") <- shObs
                attr(res,"sim_logN") <- simObs$logN
                attr(res,"sim_logF") <- simObs$logF
                attr(res,"sim_logfbar") <- simObs$logfbar
                attr(res,"sim_logssb") <- simObs$logssb
                names(res) <- sn
                return(res)
            }else{
                bionm <- names(simObs)[grepl("^(bio|total)_",names(simObs))]
                for(nn in bionm)
                    attr(allDat,nn) <- simObs[[nn]]
                attr(allDat,"shared_data") <- shObs
                attr(allDat,"sim_logN") <- simObs$logN
                attr(allDat,"sim_logF") <- simObs$logF
                attr(allDat,"sim_logfbar") <- simObs$logfbar
                attr(allDat,"sim_logssb") <- simObs$logssb
                names(allDat) <- sn
                return(allDat)
            }
        }, simplify=FALSE)
    }else{
        if(ready.to.fit)
            warning("ready.to.fit is ignored when full.data is false")
  	ret <- replicate(nsim, {a <- obj$simulate(par = par); a}, simplify=FALSE)
    }
    attr(ret,".Random.seed") <- rngSeed
    attr(ret,"m_fit") <- object
    ret
}
