

##' @method retro msam
##' @importFrom stockassessment runwithout reduce retro
##' @export
retro.msam <- function(fit, year=NULL, ncores=detectCores()-1, ...){   
    dd <- attr(fit,"m_data")

    if(length(year) != 1 || any(year > length(dd$minYearAll:dd$maxYearAll)))
        stop("year should be a single number giving the maximum number of years to exclude. Other options from stockassessment are not implemented")
 
    
    yy <- lapply(tail(dd$minYearAll:dd$maxYearAll,year), function(y) seq(y, dd$maxYearAll, 1))

    doRun <- function(years){
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
            args$newtonsteps <- attr(fit,"lower")
        if(!is.null(attr(fit,"upper")))
            args$newtonsteps <- attr(fit,"upper")
        if(!is.null(attr(fit,"dotargs")))
            args <- c(args, attr(fit,"dotargs"))
        do.call(multisam.fit, args)
    }
    if(ncores>1){
        if(!requireNamespace("parallel"))
            stop("retro can only be used with multiple cores if the parallel package is available.")
    cl <- parallel::makeCluster(ncores) #set up nodes
    on.exit(parallel::stopCluster(cl)) #shut it down
    lib.ver1 <- dirname(path.package("stockassessment"))
    lib.ver2 <- dirname(path.package("multiStockassessment"))
    parallel::clusterExport(cl, varlist=c("fit","doRun","lib.ver1","lib.ver2"), envir=environment())
    parallel::clusterEvalQ(cl, {library(stockassessment, lib.loc=lib.ver1)})
    parallel::clusterEvalQ(cl, {library(multiStockassessment, lib.loc=lib.ver2)})
    runs <- parallel::parLapply(cl, yy, function(y)doRun(y))
  } else {
    runs <- lapply(yy, function(y)doRun(y))
  }
  converg <- unlist(lapply(runs, function(x)x$opt$conv))
  if(any(converg!=0)) warning(paste0("retro run(s) ", paste0(which(converg!=0),collapse=",")," did not converge."))
  attr(runs, "fit") <- fit
  class(runs)<-c("msam_retro","msamset")
  runs

    
}



##' @method mohn msam_retro
##' @importFrom stockassessment mohn
##' @export
mohn.msam_retro <- function(fits, what=NULL, lag=0, ...){
    if(is.null(what)){
        what <- function(fit){
            recList <- rectable(fit,..., returnList = TRUE)
            ssbList <- ssbtable(fit,..., returnList = TRUE)
            fbarList <- fbartable(fit,..., returnList = TRUE)
            catchList <- catchtable(fit,..., returnList = TRUE)
            add <- 0
            dots <- list(...)
            if(!is.null(dots$lagR)){
                if(dots$lagR == TRUE){
                    add <- 1
                }
            }          
            ret <- lapply(seq_along(recList), function(i){
                r <- cbind(recList[[i]][,1], ssbList[[i]][,1], fbarList[[i]][,1], catchList[[i]][,1])
                colnames(r) <- c(paste("R(age ", fit[[i]]$conf$minAge + add, ")", sep = ""), "SSB",
                                 paste("Fbar(", fit[[i]]$conf$fbarRange[1], "-", fit[[i]]$conf$fbarRange[2], ")", sep = ""),
                                 "Catch")
                r
            })
            names(ret) <- getStockNames(fit)
            ret
        }
    }
    ref <- what(attr(fits,"fit"))
    ret <- lapply(fits, what)
    getBias <- function(i){ b <- lapply(ret, function(z){x <- z[[i]]; y<-rownames(x)[nrow(x)-lag]; (x[rownames(x)==y,]-ref[[i]][rownames(ref[[i]])==y,])/ref[[i]][rownames(ref[[i]])==y,]})
        colMeans(do.call(rbind,b))
    }
    b <- do.call("cbind",lapply(seq_along(ref), getBias))
    colnames(b) <- getStockNames(attr(fits,"fit"))
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
            args$newtonsteps <- attr(fit,"lower")
        if(!is.null(attr(fit,"upper")))
            args$newtonsteps <- attr(fit,"upper")
        if(!is.null(attr(fit,"dotargs")))
            args <- c(args, attr(fit,"dotargs"))
        mfitRed <- do.call(multisam.fit, args)

        fargs <- c(list(fit = mfitRed),
                   list(...))
        if(!match("rec.years",names(fargs),FALSE))
            fargs$rec.years <- lapply(fit, function(x) numeric(0))
        if(!any(match(c("fscale","catchval","fval","nextssb","landval","findMSY","hcr"), names(fargs), FALSE)))
            fargs$fscale = rep(NA, pmin(length(years),forecastYears))
        
        do.call(modelforecast.msam, fargs)
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
    
    attr(r,"tab") <- tab
    attr(r,"fit") <- fit
    class(r) <- c("msam_holdout")
    
    r
}
