
##' @importFrom stockassessment runwithout reduce saveConf loadConf defpar sam.fit
##' @method runwithout msam
##' @export
runwithout.msam <- function(fit, year = NULL, fleet = NULL, initializePars = TRUE, initializeRandomEffects = FALSE, silent = FALSE, ...){
    ran <- c("logN", "logF", "missing","logSW", "logCW", "logitMO", "logNM", "logP","logitFseason", "shared_logFscale","shared_missingObs", "logGst", "logGtrip")
    fake_runwithout_sam <- function (fit, year = NULL, fleet = NULL, map = fit$obj$env$map, 
                                     ...) {
        data <- stockassessment::reduce(fit$data, year = year, fleet = fleet, conf = fit$conf)
        conf <- attr(data, "conf")
        fakefile <- file()
        sink(fakefile)
        stockassessment::saveConf(conf, file = "")
        sink()
        conf <- stockassessment::loadConf(data, fakefile, patch = TRUE)
        close(fakefile)
        par <- stockassessment::defpar(data, conf)        
        par[!names(par) %in% ran] <- fit$pl[!names(fit$pl) %in% ran]
        ret <- stockassessment::sam.fit(data, conf, par, rm.unidentified = FALSE, map = map, 
                                        lower = fit$low, upper = fit$hig, run = FALSE, ...)
        class(ret) <- "sam"
        return(ret)
    }
    fitRed <- do.call("c",lapply(fit, fake_runwithout_sam, year = year, fleet = fleet, ...))
    args <- as.list(attr(fit,"m_call"))[-1]
    args$x <- fitRed
    ## args <- list(x = fitRed,
    ##              formula = formula(attr(attr(fit,"m_data")$X,"terms")),
    ##              corStructure = attr(fit,"corStructure"),
    ##              usePartialCors=attr(fit,"partialCors"),
    ##              newtonsteps = 0,
    ##              lower = list(RE = -10),
    ##              upper = list(RE = 10))
    ## ## Reduced shared observation
    shObs <- attr(fit,"m_data")$sharedObs
    if(shObs$hasSharedObs){
        shRed <- stockassessment::reduce(shObs, year = year, fleet = fleet)
        args$shared_data <- shRed[names(shObs)]
    }        
    ## ##
    ## if(!is.null(attr(fit,"newtonsteps")))
    ##     args$newtonsteps <- attr(fit,"newtonsteps")
    ## if(!is.null(attr(fit,"lower")))
    ##     args$lower <- attr(fit,"lower")
    ## if(!is.null(attr(fit,"upper")))
    ##     args$upper <- attr(fit,"upper")
    ## if(!is.null(attr(fit,"dotargs")))
    ##     args <- c(args, attr(fit,"dotargs"))
    ## args$rm.unidentified <- FALSE
    oldDat <- attr(fit,"m_data")
    oldPars <- attr(fit,"m_pl")
    ##
    defPars <- multiStockassessment:::collect_pars(fitRed)
    if(initializePars){
        nn <- setdiff(intersect(names(defPars),names(oldPars)),ran)
        defPars[nn] <- oldPars[nn]
    }
    if(initializeRandomEffects){
        vv <- multiStockassessment:::splitParameter(oldPars$logF)
        defPars$logF <- multiStockassessment:::combineParameter(lapply(seq_along(vv), function(i) vv[[i]][,!(oldDat$sam[[i]]$years %in% year),drop=FALSE]))
        vv <- multiStockassessment:::splitParameter(oldPars$logN)
        defPars$logN <- multiStockassessment:::combineParameter(lapply(seq_along(vv), function(i) vv[[i]][,!(oldDat$sam[[i]]$years %in% year),drop=FALSE]))
    }
    ##
    args$parlist <- defPars
    args$silent <- silent

    try({do.call(multisam.fit, args, envir = attr(fit,"m_envir"))})
}


##' @method retro msam
##' @importFrom stockassessment runwithout reduce retro
##' @importFrom parallel detectCores makeCluster stopCluster clusterExport clusterEvalQ
##' @export
retro.msam <- function(fit, year=NULL, ncores=parallel::detectCores()-1, initializePars = TRUE, initializeRandomEffects = FALSE, ...){   
    dd <- attr(fit,"m_data")

    if(length(year) != 1 || any(year > length(dd$minYearAll:dd$maxYearAll)))
        stop("year should be a single number giving the maximum number of years to exclude. Other options from stockassessment are not implemented")
 
    
    yy <- lapply(tail(dd$minYearAll:dd$maxYearAll,year), function(y) seq(y, dd$maxYearAll, 1)) 
    doRun <- function(years){
        runwithout(fit, year = years, initializePars=initializePars, initializeRandomEffects=initializeRandomEffects)
    }
    if(ncores>1){
        if(!requireNamespace("parallel"))
            stop("retro can only be used with multiple cores if the parallel package is available.")
    cl <- parallel::makeCluster(ncores) #set up nodes
    on.exit(parallel::stopCluster(cl)) #shut it down
    lib.ver1 <- dirname(path.package("stockassessment"))
    lib.ver2 <- dirname(path.package("multiStockassessment"))
        parallel::clusterExport(cl, varlist=c("fit","doRun","lib.ver1","lib.ver2"), envir=environment())
               parallel::clusterExport(cl, varlist=ls(attr(fit,"m_envir")), envir=attr(fit,"m_envir"))
    parallel::clusterEvalQ(cl, {library(stockassessment, lib.loc=lib.ver1)})
    parallel::clusterEvalQ(cl, {library(multiStockassessment, lib.loc=lib.ver2)})
    runs <- parallel::parLapply(cl, yy, function(y)doRun(y))
  } else {
    runs <- lapply(yy, function(y)doRun(y))
  }
  ## converg <- unlist(lapply(runs, function(x)x$opt$conv))
  ## if(any(converg!=0)) warning(paste0("retro run(s) ", paste0(which(converg!=0),collapse=",")," did not converge."))
  attr(runs, "fit") <- fit
  class(runs)<-c("msam_retro","msamset")
  runs

    
}



##' @method mohn msam_retro
##' @importFrom stockassessment mohn
##' @export
mohn.msam_retro <- function(fits, what=NULL, lag=0, addTotal = FALSE,...){
    if(is.null(what)){
        what <- function(fit,...){
            ## if(class(fit)=="try-error"){
            ##     r <- matrix(NA_real_, 4)
            ##     colnames(r) <- c("R","SSB","Fbar","Catch")
            ## }
            recList <- rectable(fit,..., returnList = TRUE, addTotal = addTotal)
            ssbList <- ssbtable(fit,..., returnList = TRUE, addTotal = addTotal)
            fbarList <- fbartable(fit,..., returnList = TRUE, addTotal = addTotal)
            catchList <- catchtable(fit,..., returnList = TRUE, addTotal = addTotal)
            add <- 0
            dots <- list(...)
            if(!is.null(dots$lagR)){
                if(dots$lagR == TRUE){
                    add <- 1
                }
            }
            sn <- getStockNames(fit)
            if(length(recList) > length(sn))
                sn <- c(sn,"Total")                 
            ret <- lapply(seq_along(recList), function(i){
                ##r <- cbind(recList[[i]][,1], ssbList[[i]][,1], fbarList[[i]][,1], catchList[[i]][,1])
                r <- cbindYearTables(list(recList[[i]][,1, drop=FALSE],
                                          ssbList[[i]][,1, drop=FALSE],
                                          fbarList[[i]][,1, drop=FALSE],
                                          catchList[[i]][,1, drop=FALSE]))
                if(i <= length(fit)){
                    colnames(r) <- c(paste("R(age ", fit[[i]]$conf$minAge + add, ")", sep = ""), "SSB",
                                     paste("Fbar(", fit[[i]]$conf$fbarRange[1], "-", fit[[i]]$conf$fbarRange[2], ")", sep = ""),
                                     "Catch")
                }else{
                    colnames(r) <- c("R","SSB","Fbar","Catch")
                }
                r
            })
            names(ret) <- sn
            ret
        }
    }
    ref <- what(attr(fits,"fit"),...)
    ret <- lapply(fits, what, ...)
    getBias <- function(i){ b <- lapply(ret, function(z){x <- z[[i]]; y<-rownames(x)[nrow(x)-lag]; (x[rownames(x)==y,]-ref[[i]][rownames(ref[[i]])==y,])/ref[[i]][rownames(ref[[i]])==y,]})
        colMeans(do.call(rbind,b), na.rm=TRUE)
    }
    b <- do.call("cbind",lapply(seq_along(ref), getBias))
    sn <- getStockNames(attr(fits,"fit"))
    if(ncol(b) == length(sn)){
        colnames(b) <- sn
    }else{
        colnames(b) <- c(sn,"Total")
    }
    b    
}


##' Holdout forecast for validation
##'
##' @param fit msam fit
##' @param nYears Number of years to hold out
##' @param ... Arguments passed to modelforecast
##' @return A model forecast
##' @author Christoffer Moesgaard Albertsen
##' @export
holdout <- function(fit, nYears, ...){
        UseMethod("holdout")
}

##' @rdname holdout
##' @method holdout msam
##' @export
holdout.msam <- function(fit, nYears, forecastYears = 1, ncores = 1, ...){
    stop("This function needs to be updated!")
    dd <- attr(fit,"m_data")
    if(forecastYears > nYears)
        stop("forecastYears must be less than or equal to nYears.")

    if(nYears == 0){
        yy <- list(tail(dd$minYearAll:dd$maxYearAll,forecastYears))
    }else{
        yy <- lapply(tail(dd$minYearAll:dd$maxYearAll,nYears),
                     function(y) seq(y, dd$maxYearAll, 1))
    }

    doOne <- function(years){
        fitRed <- do.call("c",lapply(fit, runwithout, year = years))
        args <- list(x = fitRed,
                     formula = formula(attr(attr(fit,"m_data")$X,"terms")),
                     corStructure = attr(fit,"corStructure"),
                     usePartialCors=attr(fit,"partialCors"),
                     newtonsteps = 0,
                     lower = list(RE = -10),
                     upper = list(RE = 10))
        if(!is.null(attr(fit,"newtonsteps")))
            args$newtonsteps <- attr(fit,"newtonsteps")
        if(!is.null(attr(fit,"lower")))
            args$lower <- attr(fit,"lower")
        if(!is.null(attr(fit,"upper")))
            args$upper <- attr(fit,"upper")
        if(!is.null(attr(fit,"dotargs")))
            args <- c(args, attr(fit,"dotargs"))
        mfitRed <- do.call(multisam.fit, args)

        fargs <- c(list(fit = mfitRed),
                   list(...))
        if(!match("rec.years",names(fargs),FALSE))
            fargs$rec.years <- lapply(fit, function(x) numeric(0))
        if(!any(match(c("fscale","catchval","fval","nextssb","landval","findMSY","hcr"), names(fargs), FALSE)))
            fargs$fscale = rep(NA, pmin(length(years),forecastYears))
        
        do.call(multiStockassessment:::modelforecast.msam, fargs)
    }

    obs <- catchtable(fit, obs.show = TRUE, returnList = TRUE)

    if(ncores>1){
        if(!requireNamespace("parallel"))
            stop("holdout can only be used with multiple cores if the parallel package is available.")
        cl <- parallel::makeCluster(ncores) #set up nodes
        on.exit(parallel::stopCluster(cl)) #shut it down
        lib.ver1 <- dirname(path.package("stockassessment"))
        lib.ver2 <- dirname(path.package("multiStockassessment"))
        parallel::clusterExport(cl, varlist=c("fit","doOne","lib.ver1","lib.ver2"), envir=environment())
        parallel::clusterEvalQ(cl, {library(stockassessment, lib.loc=lib.ver1)})
        parallel::clusterEvalQ(cl, {library(multiStockassessment, lib.loc=lib.ver2)})
        runs <- parallel::parLapply(cl, yy, function(y)doOne(y))
    } else {
        runs <- lapply(yy, function(y)doOne(y))
    }
    
    tab <- lapply(seq_along(fit), function(i){
        rCI <- do.call("rbind",lapply(runs, function(z){
            x <- attr(z[[i]],"tab")
            y <- rownames(x)[-1]
            v <- rep(NA_real_, forecastYears)
            v[1:length(y)] <- x[y,"catch:low"] <= obs[[i]][y,"sop.catch"] &
                obs[[i]][y,"sop.catch"] <= x[y,"catch:high"]
            names(v) <- seq_len(forecastYears)
            v
        }))
        r <- do.call("rbind",lapply(runs, function(z){
            x <- attr(z[[i]],"tab")
            y <- rownames(x)[-1]
            v <- rep(NA_real_, forecastYears)
            v[1:length(y)] <- log(x[y,"catch:median"]) - log(obs[[i]][y,"sop.catch"])
            names(v) <- seq_len(forecastYears)
            matrix(v, nrow = 1)
        }))
        data.frame(RMSLPE = apply(r,2, function(x) sqrt(mean(x^2, na.rm = TRUE))),
                   MALPE = apply(r, 2, function(x) mean(abs(x), na.rm = TRUE)),
                   Mean = colMeans(r),
                   inCI = colMeans(rCI),
                   N = apply(r,2,function(x) sum(!is.na(x)))
                   )
    })
    names(tab) <- getStockNames(fit)
    
    attr(runs,"tab") <- tab
    attr(runs,"fit") <- fit
    class(runs) <- c("msam_holdout")
    
    runs
}


as_samset <- function(x){
    x <- x[]
    class(x) <- "samset"
    x
}

##' Jitter runs 
##' @param fit a fitted model object as returned from sam.fit
##' @param nojit a list of vectors. Each element in the list specifies a run where the fleets mentioned are omitted
##' @param par initial values to jitter around. The defaule ones are returned from the defpar function
##' @param sd the standard deviation used to jitter the initial values (most parameters are on a log scale, so similar to cv) 
##' @param ncores the number of cores to attemp to use
##' @return A "samset" object, which is basically a list of sam fits
##' @details ...
##' @importFrom parallel detectCores makeCluster clusterEvalQ parLapply stopCluster
##' @importFrom stats rnorm
##' @importFrom stockassessment jit
##' @export
jit.msam <- function(fit, nojit=10, par=NULL, sd=.25, ncores=parallel::detectCores()-1, silent = TRUE){
    if(is.null(par))
        par <- multiStockassessment:::collect_pars(as_samset(fit))
    pars <- lapply(1:nojit, function(i){
        lapply(par, function(xx) xx + rnorm(length(xx),sd=sd))
    })
    ##relist(parv+rnorm(length(parv),sd=sd), par))
    doRun <- function(p){
        args <- as.list(attr(fit,"m_call"))[-1]
        args$parlist <- p
        args$silent <- silent
        try({do.call(multisam.fit, args, envir = attr(fit,"m_envir"))})
    }
    if(ncores>1){
        cl <- parallel::makeCluster(ncores) #set up nodes
        on.exit(parallel::stopCluster(cl)) #shut it down
        lib.ver1 <- dirname(path.package("stockassessment"))
        lib.ver2 <- dirname(path.package("multiStockassessment"))
        parallel::clusterExport(cl, varlist=c("fit","doRun","lib.ver1","lib.ver2"), envir=environment())
        parallel::clusterExport(cl, varlist=ls(attr(fit,"m_envir")), envir=attr(fit,"m_envir"))
        parallel::clusterEvalQ(cl, {library(stockassessment, lib.loc=lib.ver1)})
        parallel::clusterEvalQ(cl, {library(multiStockassessment, lib.loc=lib.ver2)})
        fits <- parallel::parLapply(cl, pars, doRun) 
    } else {
        fits <- lapply(pars, doRun)
    }
    attr(fits,"fit") <- fit
    attr(fits,"jitflag") <- 1
    class(fits) <- c("samset")
    fits
}
