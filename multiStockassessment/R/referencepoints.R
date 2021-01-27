##' Calculate jacobian of a function
##'
##' @param func function
##' @param x parameter values
##' @param ... passed to func
##' @return jacobian matrix
##' @author Christoffer Moesgaard Albertsen
jacobian <- function(func, x,
                     h = abs(1e-04 * x) + 1e-04 * (abs(x) < sqrt(.Machine$double.eps/7e-07)),
                     ## 0.1 * 10^floor(log10(abs(x))) + 1e-4,
                     ...){
         r <- .Call("jacobian",
                   function(x)func(x,...),
                   x,
                   globalenv(),
                   30L,
                   h,##abs(1e-4 * x) + 1e-4 * (abs(x) < 1e-8),
                   1e-12,
                   PACKAGE = "stockassessment")
        do.call("cbind",r[-1])
}
##' Calculate gradient of a function
##'
##' @param func function
##' @param x parameter values
##' @param ... passed to func
##' @return gradient vector
##' @author Christoffer Moesgaard Albertsen
grad <- function(func, x,
                 h = abs(1e-04 * x) + 1e-04 * (abs(x) < sqrt(.Machine$double.eps/7e-07)),
                 ##0.1 * 10^floor(log10(abs(x))) + 1e-4,
                 ...){
         r <- .Call("jacobian",
                   function(x)func(x,...),
                   x,
                   globalenv(),
                   30L,
                   h,##abs(1e-4 * x) + 1e-4 * (abs(x) < 1e-8),
                   1e-12,
                   PACKAGE = "stockassessment")
         v <- do.call("cbind",r[-1])
         if(nrow(v) == 1)
             return(as.vector(v))
         return(diag(v))         
}

svd_solve <- function(x){
    ss <- svd(x)
    ss$v %*% diag(1/ss$d, length(ss$d), length(ss$d)) %*% t(ss$u)
}


##' @method forecastMSY msam
##' @inheritParams stockassessment::forecastMSY
##' @importFrom stockassessment forecastMSY
##' @export
forecastMSY.msam <- function(fit,
                             nYears = 100,
                             nlminb.control = list(eval.max = 2000, iter.max = 2000, trace = 1),
                             rec.years = c(),
                             processNoiseF = FALSE,
                             jacobianHScale = 0.5,
                             nCatchAverageYears = 20,
                             ...){
    argsIn <- modelforecast.msam(fit,
                            findMSY = rep(1,nYears),
                            rec.years = rec.years,
                            processNoiseF = processNoiseF,
                            nCatchAverageYears = nCatchAverageYears)
    argsIn$DLL <- "multiStockassessment"

    args <- argsIn
    args$map$logFScaleMSY <- factor(rep(NA,length(fit)))
   
    objForecast <- do.call(TMB::MakeADFun, args)
    objForecast$fn()

    ## Get joint precision
    jointPrecision <- TMB::sdreport(objForecast,
                                    attr(fit,"m_opt")$par,
                                    svd_solve(attr(fit,"m_sdrep")$cov.fixed),
                                    getJointPrecision = TRUE)$jointPrecision

    rp <- c("logFScaleMSY")


    ## Find MSY value
    args <- argsIn
    args$parameters$logFScaleMSY[] <- -1
    args$parameters$implicitFunctionDelta[] <- 0
    map0 <- args$map
    fix <- setdiff(names(args$parameters), args$random)
    args$map <- lapply(args$parameters[fix], function(x)factor(x*NA))
    args$map <- args$map[names(args$map) != rp]
    args$map$implicitFunctionDelta <- NULL

    
    objOptim <- do.call(TMB::MakeADFun, args)

    objOptim$fn(objOptim$par)
    ## obj$gr(obj$par)

    fn <- Vectorize(function(x){
        objOptim$fn(c(x,0))
        theta <- objOptim$env$last.par
        objOptim$env$f(theta, order = 1)[1,which(names(objOptim$env$last.par) == "implicitFunctionDelta")]
    })

    ## Try different values??
    testStart <- matrix(c(-2, -1, -0.5, -0.1, 0, 0.1, 0.25, 0.3, 0.5, 1, 2), nrow = 11, ncol = length(fit))
    fnTestStart <- apply(testStart,2, fn)
    opt <- nlminb(testStart[which(order(fnTestStart) == 1)], fn, control = nlminb.control)

    ## Object to do Delta method (no map, no random, delta = 1)
    args <- argsIn
    args$parameters <- objOptim$env$parList(x = opt$par)
    args$parameters$implicitFunctionDelta <- 1
    args$map$logFScaleMSY <- NULL
    args$map$implicitFunctionDelta <- factor(rep(NA,length(fit)))
    args$random <- NULL

    objDelta <- do.call(TMB::MakeADFun, args)

    ## Get Jacobian
    gridx <- which(names(objDelta$par) %in% rp)
    JacAll <- jacobian(function(x) objDelta$gr(x)[gridx], objDelta$par, h = jacobianHScale * abs(1e-04 * objDelta$par) + 1e-04 * (abs(objDelta$par) < sqrt(.Machine$double.eps/7e-07)))

    ## Implicit function gradient
    dCdTheta <- -svd_solve(JacAll[,gridx,drop=FALSE]) %*% JacAll[,-gridx,drop=FALSE]
    rownames(dCdTheta) <- rp
    colnames(dCdTheta) <- names(objForecast$env$last.par)


    ## Do delta method
    xtra <- diag(1,length(objForecast$env$last.par.best))
    diag(xtra)[objForecast$env$random] <- 0
    dG <- rbind(xtra[diag(xtra) != 0,,drop = FALSE],dCdTheta)
    covAll <- dG %*% svd_solve(jointPrecision) %*% t(dG)
    covAllOld <- covAll

    ## Object to do sdreport (delta = 0)
    args <- argsIn
    ## Remvoe rp from map
    args$map <- args$map[!(names(args$map) %in% rp)]
    ## Use optimized parameters
    args$parameters <- objOptim$env$parList(x = opt$par)
    args$parameters$implicitFunctionDelta[] <- 0

    objSDR <- do.call(TMB::MakeADFun, args)

    sdr <- TMB::sdreport(objSDR, objSDR$par, svd_solve(covAll))
    ssdr <- summary(sdr)


    toCI <- function(what){
        exp(ssdr[rownames(ssdr) %in% what,] %*% cbind(Estimate=c(1,0),Low=c(1,-2),High=c(1,2)))
    }

    forecastTableSingle <- function(i){
        Fmsy <- toCI(sprintf("SAM_%d_logFMSY",i-1))
        ssb <- toCI(sprintf("SAM_%d_logssb",i-1))
        yield <- toCI(sprintf("SAM_%d_logCatch",i-1))
        rec <- toCI(sprintf("SAM_%d_logR",i-1))
        fbar <- toCI(sprintf("SAM_%d_logfbar",i-1))
        tab <- rbind(Fmsy, ssb[nrow(ssb),], yield[nrow(yield),], rec[nrow(rec),])
        rownames(tab) <- c("Fmsy", "Bmsy", "Yield", "Rmsy")
        rownames(ssb) <- rownames(yield) <- rownames(rec) <- rownames(fbar) <- min(fit[[i]]$data$years) + 0:(nrow(ssb)-1)
        res <- list(table = tab,
                    timeseries = list(SSB = ssb,
                                      Catch = yield,
                                      Recruitment = rec,
                                      Fbar = fbar),
                    opt = NA,
                    sdr = NA,
                    fbarlabel = substitute(bar(F)[X - Y], list(X = fit[[i]]$conf$fbarRange[1], Y = fit[[i]]$conf$fbarRange[2])))
        class(res) <- "sam_forecast_msy"
        res
    }

    res <- lapply(as.list(seq_along(fit)), forecastTableSingle)
    names(res) <- getStockNames(fit)
    attr(res,"m_opt") <- opt
    attr(res,"m_sdr") <- sdr 
    class(res) <- "msam_forecast_msy"
    res

}

##' @method referencepoints msam
##' @inheritParams stockassessment::referencepoints
##' @importFrom stockassessment referencepoints
##' @export
referencepoints.msam <- function(fit,
                                 nYears = 100,
                                 Fsequence = seq(0,4, len = 200),
                                 aveYears = lapply(fit,function(x)max(x$data$years)+(-9:0)),
                                 selYears = lapply(fit,function(x)max(x$data$years)),
                                 SPRpercent = c(0.35),
                                 catchType = "catch",
                                 jacobianHScale = 0.5,
                                 MSYreduction = c(0.05),
                                 newtonSteps = 3,
                                 optN = 20,
                                 ...){

    nStocks <- length(fit)

     
    toList <- function(x){
        if(is.null(x))
            return(replicate(nStocks,numeric(0),FALSE))
        if(is.list(x))
            return(x)
        if(is.vector(x)) ## Same target for all stocks
            x <- matrix(x, nStocks, length(x), byrow = TRUE)
        if(is.matrix(x))
            return(split(x, row(x)))
    }


    ## Get joint precision
    jointPrecision <- TMB::sdreport(attr(fit,"m_obj"),
                                    attr(fit,"m_opt")$par,
                                    svd_solve(attr(fit,"m_sdrep")$cov.fixed),
                                    getJointPrecision = TRUE)$jointPrecision

    ## Prepare arguments to calculate reference points (fix parameters and latent variables, delta = 1)
    obj0 <- attr(fit,"m_obj")
    argsIn <- as.list(obj0$env)[methods::formalArgs(TMB::MakeADFun)[methods::formalArgs(TMB::MakeADFun) != "..."]]
    argsIn$parameters <- attr(fit,"m_pl")
    argsIn$random <- unique(names(obj0$env$par[obj0$env$random]))
    ## Add referencepointSet
    catchType <- lapply(toList(catchType), pmatch, table = c("catch","landing","discard"))
    if(any(is.na(unlist(catchType))))
        stop("Invalid catch type")


    nYears <- toList(nYears)
    Fsequence <- toList(Fsequence)
    SPRpercent <- toList(SPRpercent)
    optN <- toList(optN)
    MSYreduction <- toList(MSYreduction)
    MSYfraction <- lapply(MSYreduction, function(x) 1-x)
   

    aveYears <- toList(aveYears)
    aveYears <- lapply(as.list(seq_along(fit)),function(i)match(aveYears[[i]], table =fit[[i]]$data$years) - 1)
    if(any(is.na(unlist(aveYears))))
        stop("aveYears has years without data.")

    selYears <- toList(selYears)
    selYears <- lapply(as.list(seq_along(fit)),function(i)match(selYears[[i]], table =fit[[i]]$data$years) - 1)
    if(any(is.na(unlist(selYears))))
        stop("selYears has years without data.")


    for(i in seq_along(fit))
        argsIn$data$sam[[i]]$referencepoint <- list(nYears = nYears[[i]],
                                                    aveYears = aveYears[[i]],
                                                    selYears = selYears[[i]],
                                                    Fsequence = Fsequence[[i]],
                                                    xPercent = SPRpercent[[i]],
                                                    catchType = catchType[[i]]-1,
                                                    MSYRange = MSYfraction[[i]],
                                                    optN = optN[[i]]
                                                    )

    args <- argsIn
    ## Remove random
    args$random <- NULL
    ## Remove referencepoint parameters from map
    map0 <- args$map
    fix <- names(args$parameters) ## setdiff(names(args$parameters), args$random)
    args$map <- lapply(args$parameters[fix], function(x)factor(x*NA))
    args$map$logScaleFmsyRange <- NULL

    rp <- vector("list",length=nStocks)
    for(i in seq_along(fit)){
        if(fit[[i]]$conf$stockRecruitmentModelCode %in% c(0)){ # RW
            rp[[i]] <- c("logScaleFmax",
                         "logScaleF01",
                         "logScaleFxPercent")
        }else if(fit[[i]]$conf$stockRecruitmentModelCode %in% c(61,63)){ # Hockey-sticks
            rp[[i]] <- c("logScaleFmsy",
                         "logScaleFmax",
                         "logScaleF01",
                         "logScaleFcrash",
                         "logScaleFext",
                         "logScaleFxPercent",
                         "logScaleFlim",
                         "logScaleFmsyRange")
        }else if(fit[[i]]$conf$stockRecruitmentModelCode %in% c(3,62)){ # constant mean, AR
            rp[[i]] <- c("logScaleFmsy",
                         "logScaleFmax",
                         "logScaleF01",
                         "logScaleFxPercent",
                         "logScaleFmsyRange")
        }else if(fit[[i]]$conf$stockRecruitmentModelCode %in% c(64)){ # Pow CMP
            rp[[i]] <- c("logScaleFmsy",
                         "logScaleFmax",
                         "logScaleF01",
                         "logScaleFxPercent",
                         "logScaleFmsyRange"
                         )
        }else if(fit[[i]]$conf$stockRecruitmentModelCode %in% c(65)){ # Pow Non-CMP
            rp[[i]] <- c("logScaleFmsy",
                         "logScaleFmax",
                         "logScaleF01",
                         "logScaleFxPercent",
                         "logScaleFmsyRange"
                         )
        }else if(fit[[i]]$conf$stockRecruitmentModelCode %in% c(68,69) && fit[[i]]$pl$rec_par[3] > 0){ ## depensatory recruitment; Fcrash does not work.
            rp[[i]] <- c("logScaleFmsy",
                         "logScaleFmax",
                         "logScaleF01",
                         "logScaleFext",
                         "logScaleFxPercent",
                         "logScaleFmsyRange"
                         )
        }else if(fit[[i]]$conf$stockRecruitmentModelCode %in% c(90)){
            rp[[i]] <- c("logScaleFmsy",
                         "logScaleFmax",
                         "logScaleF01",
                         "logScaleFcrash",
                         "logScaleFext",
                         "logScaleFxPercent",
                         "logScaleFmsyRange"
                         )
        }else if(fit[[i]]$conf$stockRecruitmentModelCode %in% c(91,92)){
            rp[[i]] <- c("logScaleFmsy",
                         "logScaleFmax",
                         "logScaleF01",
                         ##"logScaleFcrash",
                         "logScaleFext",
                         "logScaleFxPercent",
                         "logScaleFmsyRange"
                         )
        }else{
            rp[[i]] <- c("logScaleFmsy",
                         "logScaleFmax",
                         "logScaleF01",
                         "logScaleFcrash",
                         "logScaleFext",
                         "logScaleFxPercent",
                         "logScaleFmsyRange"
                         )
        }
    }

    rp2map <- function(rp, what){
        mtmp <- rep(NA,length(fit))
        v <- which(sapply(rp,function(x) any(x == what)))
        mtmp[v] <- seq_along(v)
        factor(mtmp)
    }

    rpMap <- list()
    rpMap$logScaleFmsy <- rp2map(rp, "logScaleFmsy")
    rpMap$logScaleF01 <- rp2map(rp, "logScaleF01")
    rpMap$logScaleFmax <- rp2map(rp, "logScaleFmax")
    rpMap$logScaleFcrash <- rp2map(rp, "logScaleFcrash")
    rpMap$logScaleFext <- rp2map(rp, "logScaleFext")
    rpMap$logScaleFlim <- rp2map(rp, "logScaleFlim")    
   
    ## Referencepoints to estimate
    args$map <- c(args$map[-which(names(args$map) %in% unique(unlist(rp)))],
                  rpMap[which(names(rpMap) %in% unique(unlist(rp)))])

    args$parameters$logScaleFmsy[] <- -1
    args$parameters$logScaleF01[] <- -1
    args$parameters$logScaleFmax[] <- -1
    args$parameters$logScaleFcrash[] <- -1
    args$parameters$logScaleFext[] <- -1
    args$parameters$logScaleFxPercent <- combineVectors(lapply(seq_along(fit), function(i) rep(-1, length(SPRpercent[[i]]))))
    if(length(MSYfraction) > 0){
        args$parameters$logScaleFmsyRange <- combineMatrices(lapply(seq_along(fit), function(i){
            if(any(rp[[i]] == "logScaleFmsyRange")){
                rbind(log(-log(MSYfraction[[i]])),
                      log(log(1+(1-MSYfraction[[i]]))))
            }else{
                matrix(0,2,0)
            }
        }
        ))
    }
    args$parameters$logScaleFlim[] <- -1
    args$parameters$implicitFunctionDelta[] <- 1

    objOptim <- do.call(TMB::MakeADFun, args)
    pStart <- vector("list", nStocks)
    ##pStart$logScaleFmsyRange <- replicate(nStocks, matrix(0,2,0), simplify = FALSE)
        
    ## Take inital look at YPR / SPR to determine if Fmax makes sense
    rep <- objOptim$report()
    tryAgain <- FALSE
    stockNames <- multiStockassessment:::getStockNames(fit)
    for(i in seq_along(stockNames)){
        logFbar <- unname(log(tail(fbartable(fit,returnList = TRUE)[[i]][,"Estimate"],1)))
        ## Fmax
        if(which.max(rep$logYPR[[i]][is.finite(rep$logYPR[[i]])]) == length(rep$logYPR[[i]][is.finite(rep$logYPR[[i]])]) && any(rp[[i]] %in% "logScaleFmax")){
            warning(sprintf("%s does not appear to have a well-defined Fmax. Fmax will not be estimated. Increase the upper bound of Fsequence to try again.",stockNames[i]))
            rp[[i]] <- rp[[i]][-which(rp[[i]] %in% "logScaleFmax")]
            args$map$logScaleFmax[i] <- NA
            args$map$logScaleFmax <- factor(args$map$logScaleFmax)
            rpMap$logScaleFmax[i] <- NA
            rpMap$logScaleFmax <- factor(rpMap$logScaleFmax)
            tryAgain <- TRUE
        }else if(any(rp[[i]] %in% "logScaleFmax")){
            indx <- which(is.finite(rep$logYPR[[i]]) & Fsequence[[i]] > 0)
            ypr <- rep$logYPR[[i]][indx]
            ff <- Fsequence[[i]][indx]
            pStart[[i]]$logScaleFmax <- log(ff[which.max(ypr)]) - logFbar
        }

        ## Fcrash
        if(any(rp[[i]] %in% "logScaleFcrash") && min(rep$logSe[[i]][is.finite(rep$logSe[[i]])], na.rm = TRUE) > -4){
            warning(sprintf("%s does not appear to have a well-defined Fcrash. Fcrash will not be estimated. Increase the upper bound of Fsequence to try again.",stockNames[i]))
            rp[[i]] <- rp[[i]][-which(rp[[i]] %in% "logScaleFcrash")]
            args$map$logScaleFcrash[i] <- NA
            args$map$logScaleFcrash <- factor(args$map$logScaleFcrash)
            rpMap$logScaleFcrash[i] <- NA
            rpMap$logScaleFcrash <- factor(rpMap$logScaleFcrash)
            tryAgain <- TRUE
        }
        
        ## Fmsy
        if(any(rp[[i]] %in% "logScaleFmsy") && which.max(rep$logYe[[i]][is.finite(rep$logYe[[i]])]) == length(rep$logYe[is.finite(rep$logYe[[i]])])){
            warning(sprintf("%s does not appear to have a well-defined Fmsy. Fmsy will not be estimated. Increase the upper bound of Fsequence to try again.",stockNames[i]))
            rp[[i]] <- rp[[i]][-which(rp[[i]] %in% "logScaleFmsy")]
            args$map$logScaleFmsy[i] <- NA
            args$map$logScaleFmsy <- factor(args$map$logScaleFmsy)
            rpMap$logScaleFmsy[i] <- NA
            rpMap$logScaleFmsy <- factor(rpMap$logScaleFmsy)
            tryAgain <- TRUE
        }else if(any(rp[[i]] %in% "logScaleFmsy")){
            indx <- which(is.finite(rep$logYe[[i]]) & Fsequence[[i]] > 0)
            ye <- rep$logYe[[i]][indx]
            ff <- Fsequence[[i]][indx]
            pStart[[i]]$logScaleFmsy <- log(ff[which.max(ye)]) - logFbar
            if(any(rp[[i]] %in% "logScaleFmsyRange")){
                indx2 <- which(is.finite(rep$logYe[[i]]) & Fsequence[[i]] > 0 & Fsequence[[i]] < ff[which.max(ye)])
                ye2 <- rep$logYe[[i]][indx2]
                ff2 <- Fsequence[[i]][indx2]
                fmsy <- ff[which.max(ye)]
                FmsyRangeLower <- sapply(MSYfraction[[i]], function(x){
                    fL <- ff2[which.min((ye2 - x * max(ye))^2)]
                    -log(-(log(fL) - log(fmsy)))
                })
                indx3 <- which(is.finite(rep$logYe[[i]]) & Fsequence[[i]] > 0 & Fsequence[[i]] > ff[which.max(ye)])
                ye3 <- rep$logYe[[i]][indx3]
                ff3 <- Fsequence[[i]][indx3]
                fmsy <- ff[which.max(ye)]
                FmsyRangeUpper <- sapply(MSYfraction[[i]], function(x){
                    fL <- ff3[which.min((ye3 - x * max(ye))^2)]
                    log(log(fL) - log(fmsy))
                })
                pStart[[i]]$logScaleFmsyRange <- rbind(FmsyRangeLower,
                                                  FmsyRangeUpper)
            }
        }
        
        ## F01
        if(any(rp[[i]] %in% "logScaleF01")){
            indx <- which(is.finite(rep$logYPR[[i]]) & Fsequence[[i]] > 0)
            ypr <- rep$logYPR[[i]][indx]
            ff <- Fsequence[[i]][indx]
            pStart[[i]]$logScaleF01 <- log(ff[which.min((diff(ypr)/diff(ff) - 0.1 * diff(ypr)[1] / diff(ff)[1])^2)]) - logFbar
        }
        
        ## Fx%
        if(any(rp[[i]] %in% "logScaleFxPercent")){
            indx <- which(is.finite(rep$logSPR[[i]]) & Fsequence[[i]] > 0)
            spr <- rep$logSPR[[i]][indx]
            ff <- Fsequence[[i]][indx]
            pStart[[i]]$logScaleFxPercent <- sapply(SPRpercent[[i]],function(x){
                log(ff[which.min((spr - x * spr[1])^2)]) - logFbar
            })
        }
    }

    if(tryAgain)                        
        objOptim <- do.call(TMB::MakeADFun, args)
    
    p0 <- objOptim$par    
    for(ii in unique(unlist(lapply(pStart,names))))
        p0[names(p0) %in% ii] <- combineParameter(lapply(pStart,function(x)x[[match(ii,names(x))]]))

    
    opt <- nlminb(p0, objOptim$fn, objOptim$gr, objOptim$he)
    for(ii in seq_len(newtonSteps)){
        g <- as.numeric( objOptim$gr(opt$par) )
        h <- objOptim$he(opt$par)
        opt$par <- opt$par - solve(h, g)
        opt$objective <- objOptim$fn(opt$par)
    }

    ## Object to do Delta method (nothing mapped (that's not mapped in fit$obj, nothing random, delta = 1)
    args <- argsIn
    ## Remove random
    args$random <- NULL
    ## Remvoe rp from map
    ## args$map <- args$map[!(names(args$map) %in% unique(unlist(rp)))]
    ## Remember to handle different reference points for different stocks
    args$map <- c(args$map[-which(names(args$map) %in% unique(unlist(rp)))],
                  rpMap)

    ## Use optimized parameters
    args$parameters <- objOptim$env$parList(x = opt$par)

    objDelta <- do.call(TMB::MakeADFun, args)

    ## Get Jacobian
    gridx <- which(names(objDelta$par) %in% unique(unlist(rp)))
    JacAll <- jacobian(function(x) objDelta$gr(x)[gridx], objDelta$par, h = jacobianHScale * abs(1e-04 * objDelta$par) + 1e-04 * (abs(objDelta$par) < sqrt(.Machine$double.eps/7e-07)))

    ## Implicit function gradient
    dCdTheta <- -svd_solve(JacAll[,gridx,drop=FALSE]) %*% JacAll[,-gridx,drop=FALSE]
    rownames(dCdTheta) <- gsub("^logScale","",names(objDelta$par)[gridx])
    colnames(dCdTheta) <- names(fit$obj$env$last.par)

    ## Do delta method
    xtra <- diag(1,length(attr(fit,"m_obj")$env$last.par.best))
    diag(xtra)[attr(fit,"m_obj")$env$random] <- 0
    dG <- rbind(xtra[diag(xtra) != 0,,drop = FALSE],dCdTheta)
    covAll <- dG %*% svd_solve(jointPrecision) %*% t(dG)
    covAllOld <- covAll

    ## Object to do sdreport (delta = 0)
    args <- argsIn
    ## Remvoe rp from map
    ## args$map <- args$map[!(names(args$map) %in% unique(unlist(rp)))]
    args$map <- c(args$map[-which(names(args$map) %in% unique(unlist(rp)))],
                  rpMap)
    ## Use optimized parameters
    args$parameters <- objOptim$env$parList(x = opt$par)
    args$parameters$implicitFunctionDelta[] <- 0

    objSDR <- do.call(TMB::MakeADFun, args)


    sdr <- TMB::sdreport(objSDR, objSDR$par, svd_solve(covAll))
    ssdr <- summary(sdr)

    toCI <- function(what){
        tmp <- ssdr[rownames(ssdr) == what,,drop=FALSE]
        CI <- exp(tmp %*% cbind(CIL=c(1,-2),CIH=c(1,2)))
        cbind(Estimate = exp(tmp[,1]), CI)
    }

    rpTableSingle <- function(i){
        Ftab <- do.call("rbind",sapply(unique(rownames(ssdr)[grepl(sprintf("SAM_%d_referencepoint.logF",i-1),rownames(ssdr))]),
                                       toCI, simplify = FALSE))
        Btab <- do.call("rbind",sapply(unique(rownames(ssdr)[grepl(sprintf("SAM_%d_referencepoint.logB",i-1),rownames(ssdr))]),
                                       toCI, simplify = FALSE))
        Rtab <- do.call("rbind",sapply(unique(rownames(ssdr)[grepl(sprintf("SAM_%d_referencepoint.logR",i-1),rownames(ssdr))]),
                                       toCI, simplify = FALSE))
        Ytab <- do.call("rbind",sapply(unique(rownames(ssdr)[grepl(sprintf("SAM_%d_referencepoint.logY",i-1),rownames(ssdr)) & !grepl("referencepoint.logYPR",rownames(ssdr))]),
                                       toCI, simplify = FALSE))
        SPRtab <- do.call("rbind",sapply(unique(rownames(ssdr)[grepl(sprintf("SAM_%d_referencepoint.logSPR",i-1),rownames(ssdr))]),
                                         toCI, simplify = FALSE))
        YPRtab <- do.call("rbind",sapply(unique(rownames(ssdr)[grepl(sprintf("SAM_%d_referencepoint.logYPR",i-1),rownames(ssdr))]),
                                         toCI, simplify = FALSE))
        colnames(Ftab) <- colnames(Btab) <- colnames(Rtab) <- colnames(Ytab) <- colnames(SPRtab) <- colnames(YPRtab) <- c("Estimate","Low","High")

        toRowNames <- Vectorize(function(x){
            switch(gsub("^SAM_[[:digit:]]+_referencepoint.logF","",x),
                   "sq"="Status quo",
                   "0"="Zero catch",
                   "msy"="MSY",
                   "msyRange"="xMR",
                   "max"="Max",
                   "01"="0.1",
                   "crash"="Crash",
                   "ext"="Ext",
                   "xPercent"="xP",
                   "lim"="lim",
                   x
                   )               
        })
        rn <- toRowNames(rownames(Ftab))
        rn[which(rn == "xP")] <- sapply(SPRpercent[[i]],function(x)sprintf("%s%%",x * 100))
        rn[which(rn == "xMR")] <- sapply(MSYreduction,function(x){
            c(sprintf("MSY %s%% reduction range (Lower)",x * 100),
              sprintf("MSY %s%% reduction range (Upper)",x * 100))
        })

        rownames(Ftab) <- rownames(Btab) <- rownames(Rtab) <- rownames(Ytab) <- rownames(SPRtab) <- rownames(YPRtab) <- unname(rn)

        ## Ftab["Ext",c("Low","High")] <- NA
        ## Btab["Ext",c("Low","High")] <- NA
        ## Rtab["Ext",c("Low","High")] <- NA
        ## Ytab["Ext",c("Low","High")] <- NA
        ## SPRtab["Ext",c("Low","High")] <- NA
        ## YPRtab["Ext",c("Low","High")] <- NA

        YPRseq <- toCI(sprintf("SAM_%d_logYPR",i-1))
        SPRseq <- toCI(sprintf("SAM_%d_logSPR",i-1))
        Yieldseq <- toCI(sprintf("SAM_%d_logYe",i-1))
        Bseq <- toCI(sprintf("SAM_%d_logSe",i-1))
        Rseq <- toCI(sprintf("SAM_%d_logRe",i-1))

        rownames(YPRseq) <- rownames(SPRseq) <- rownames(Yieldseq) <- rownames(Bseq) <- rownames(Rseq) <- argsIn$data$referencepoint$Fsequence

        res <- list(tables = list(F = Ftab,
                                  Yield = Ytab,
                                  YieldPerRecruit = YPRtab,
                                  SpawnersPerRecruit = SPRtab,
                                  Biomass = Btab,
                                  Recruitment = Rtab
                                  ),
                    graphs = list(F = Fsequence[[i]],
                                  Yield = Yieldseq,
                                  YieldPerRecruit = YPRseq,
                                  SpawnersPerRecruit = SPRseq,
                                  Biomass = Bseq,
                                  Recruitment = Rseq),
                    fbarlabel = substitute(bar(F)[X - Y], list(X = fit[[i]]$conf$fbarRange[1], Y = fit[[i]]$conf$fbarRange[2]))## ,
                    ## diagonalCorrection = tv
                    )
        class(res) <- "sam_referencepoints"
        return(res)
    }

    res <- lapply(as.list(seq_along(fit)), rpTableSingle)
    names(res) <- getStockNames(fit)
    attr(res,"m_opt") <- opt
    attr(res,"m_sdr") <- sdr 
    class(res) <- "msam_referencepoints"
    res
    
}



##' @method addRecruitmentCurve msam
##' @inheritParams stockassessment::addRecruitmentCurve
##' @importFrom stockassessment addRecruitmentCurve
##' @export
addRecruitmentCurve.msam <- function(fit,
                    CI = TRUE,
                    col = rgb(0.6,0,0),
                    cicol = rgb(0.6,0,0,0.3),
                    plot = TRUE,
                    ...){
    
}
