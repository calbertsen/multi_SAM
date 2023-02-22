.timeToSeed <- getFromNamespace(".timeToSeed","stockassessment")
.SAM_setPB <- getFromNamespace(".SAM_setPB","stockassessment")
.SAMpb <- getFromNamespace(".SAMpb","stockassessment")
predict.rpscurvefit <- getFromNamespace("predict.rpscurvefit","stockassessment")

.getDoSim <- function (logf1, fit, nYears, aveYears, selYears, pl, logCustomSel = NULL, 
    constraint = "F=%f", deterministicF = TRUE, ...) 
{
    sel <- replicate(length(fit),numeric(0))
    for(ss in seq_along(fit)){
        if (length(logCustomSel[[ss]]) > 0) {
            sel <- exp(logCustomSel[[ss]])
        } else if (length(selYears[[ss]]) == 0) {
            sel <- NULL
        } else {
            .logFtoSel <- getFromNamespace(".logFtoSel","stockassessment")            
            logFA <- splitParameter(pl$logF)
            logF <- logFA[[ss]]
            dd <- attr(fit,"m_data")
            if(dd$shared_F_type %in% c(-1,2,5)){
                logF <- logFA[[1]]
            }else if(dd$shared_F_type %in% c(3,4,6) && ss > 1){
                Xph <- dd$Xph[[ss]]
                phb <- splitParameter(pl$shared_phbeta)[[ss]]
                predPh <- matrix(Xph %*% phb,dd$maxAgeAll-dd$minAgeAll+1)
                lfAdd <- logF*0
                for(f in seq_len(nrow(fit[[1]]$conf$keyLogFsta)))
                    for(a in seq_len(ncol(fit[[1]]$conf$keyLogFsta)))
                        if(fit[[1]]$conf$keyLogFsta[f,a] > (-1))
                            lfAdd[fit[[1]]$conf$keyLogFsta[f,a],] <- lfAdd[fit[[1]]$conf$keyLogFsta[f,a],] + predPh[a,]
                logF <- logFA[[1]] + lfAdd
            }else if(dd$shared_F_type == 1){
                stop("Not implemented yet")
            }
            sel <- exp(.logFtoSel(logF, match(selYears[[ss]], fit[[ss]]$data$years) - 1, fit[[ss]]$conf))
        }
    }
    suppressWarnings(invisible(doSim <- modelforecast(fit, rep(sprintf(constraint, 
        exp(logf1)), nYears), nosim = 1, progress = FALSE, returnObj = 2, 
        ave.years = aveYears, custom_pl = pl, cstomSel = sel, 
        deterministicF = deterministicF, fastFixedF = grepl("[[:space:]]*F[[:space:]]*=[[:space:]]*%f[[:space:]]*", 
            constraint), ...)))
    doSim
}

.perRecruitSR <- function(logf, fit, nYears, aveYears, selYears, pl = attr(fit,"m_pl"), ct=0, logCustomSel = NULL, nTail = 1,incpb=NULL, doSim = NULL, label="", constraint = "F=%f", deterministicF=TRUE, ...){    
    ## pl$missing <- NULL
    ## attr(pl,"what") <- NULL
    ## f2 <- sam.fit(fit$data,fit$conf,pl, run = FALSE)
    ## f2$opt <- list(par = f2$obj$par, objective=NA)
    ## f2$sdrep <- fit$sdrep[c("estY","covY")]
    ## class(f2) <- "sam"
    if(is.null(doSim))
        doSim <- .getDoSim(logf[1], fit, nYears, aveYears, selYears, pl, logCustomSel, deterministicF)
    if(isTRUE(all.equal(pl,fit$pl))){
        re_pl <- NULL
    }else{
        re_pl <- pl
    }

    ## progress <- TRUE
    ## if(progress){
    ##     pb <- .SAMpb(min = 0, max = length(logf))
    ##     incpb <- function() .SAM_setPB(pb, pb$getVal()+1)
    ## }else{
    ##     incpb <- function(){}
    ## }

    do.call("rbind",lapply(logf, function(lf){
        v <- doSim(sprintf(constraint,exp(lf)))
        rr <- lapply(seq_along(fit), function(ss){
            logRe <- tail(v$logN[[ss]][1,],nTail)
            logSe <- tail(na.omit(v$logssb[[ss]]),nTail)
            logSPR <- tail(na.omit(v$logEmpiricalSPR[[ss]]),nTail)
            if(ct == 0){
                logYe <- tail(v$logCatch[[ss]],nTail)
                logYPR <- tail(na.omit(v$logEmpiricalYPR[[ss]]),nTail)
            }else if(ct == 1){
                logYe <- tail(v$logLand[[ss]],nTail)
                logYPR <- tail(na.omit(v$logEmpiricalYPR_L[[ss]]),nTail)
            }else{
                logYe <- tail(log(exp(v$logCatch[[ss]])-(v$logLand[[ss]])),nTail)
                logYPR <- tail(v$logEmpiricalYPR_D[[ss]],nTail)
            }
            logLE <- NA_real_ # tail(v$logLifeExpectancy,nTail)
            logYL <- NA_real_ #tail(v$logYLTF,nTail)

            res <- list(logF = rep(lf,length.out = nTail),
                        logYPR = logYPR,
                        logSPR = logSPR,
                        logSe = logSe,
                        logRe = logRe,
                        logYe = logYe,
                        dSR0  = NA_real_,
                        logLifeExpectancy = logLE,
                        logYearsLost = logYL,
                        logDiscYe = NA_real_,
                        logDiscYPR = NA_real_)
            if(!is.null(incpb))
                incpb(label)
            as.data.frame(res)
        })
        names(rr) <- getStockNames(fit)
        rr
        }))
}




#' @importFrom stats runif predict
.refpointSFitCriteria <- function(rpArgs, pl, MT, fit, nosim, Frange, aveYears, selYears, nYears, catchType, nTail = 1,doSim=NULL,incpb=NULL,label="",constraint="F=%f", deterministicF = TRUE, ...){
    .refpointSCurveFit <- getFromNamespace(".refpointSCurveFit","stockassessment")
    rfv <- function(n,a,b){
        u <- stats::runif(n)
        v1 <- exp(stats::runif(n,log(ifelse(a==0,0.002,a)), log(b)))
        v2 <- stats::runif(n,a,b)
        ifelse(u < 0.25, v1, v2)
    }
    Fvals <- sort(rfv(nosim,Frange[1],Frange[2]))
    PRvals <- .perRecruitSR(log(Fvals),
                            fit=fit,
                            nYears=nYears,
                            aveYears = aveYears,
                            selYears = selYears,
                            pl = pl,
                            ct = catchType,
                            nTail = nTail,
                            incpb = incpb,
                            doSim = doSim,
                            label=sprintf("%s equilibrium simulations",label),
                            constraint = constraint,
                            deterministicF = deterministicF,
                            ...)  
###### Different for different RP's
    getOneRP <- function(rp, ss){
        if(rp$rpType == 1){ ## MSY
            Crit <- exp(PRvals[[ss]]$logYe)
            cutfun <- function(x) x > max(x) * rp$cutoff
            trans <- function(x, report=FALSE, ...){
                v <- exp(x)
                if(report)
                    names(v) <- "MSY"
                v
            }
            fn <- function(x) -predict(CurveFit,trans(x))
            startVals <- function() log(Fseq[which.max(pv)])
        }else if(rp$rpType == 2){ ## MSYRange
            Crit <- exp(PRvals[[ss]]$logYe)
            cutfun <- function(x) x > max(x) * rp$cutoff
            trans <- function(x, report=FALSE, ...){
                dots <- list(...)
                if("keepMSY" %in% names(dots)){
                    keepMSY <- dots$keepMSY
                }else{
                    keepMSY <- FALSE
                }
                xMSY <- exp(utils::head(x,1))
                x2 <- matrix(utils::tail(x,-1),2)
                x3 <- rbind(exp(log(xMSY) - exp(-x2[1,])),
                            exp(log(xMSY) + exp(x2[2,])))
                if(keepMSY)
                    return(c(xMSY, x3))
                if(report){
                    nm <- as.vector(t(outer(paste0(formatC(rp$xVal),"MSYRange"),c("Lower","Upper"),paste)))
                    xOut <- as.vector(x3)
                    names(xOut) <- nm
                    return(xOut)
                }
                x3
            }
            fn <- function(x){
                xx <- trans(x, keepMSY=TRUE)
                p <- predict(CurveFit,as.vector(xx))
                sum((tail(p,-1) - rep(rp$xVal,each=2) * p[1])^2) - p[1]
            }
            startVals <- function(){
                fmsy <- Fseq[which.max(pv)]
                c2 <- sapply(rp$xVal, function(xx) (pv - xx*max(pv))^2)
                f0 <- apply(c2,2,function(cc){
                    fl <- Fseq[Fseq < fmsy][which.min(cc[Fseq < fmsy])]
                    fu <- Fseq[Fseq > fmsy][which.min(cc[Fseq > fmsy])]
                    ##c(log(fmsy-fl),log(fu-fmsy))
                    c(-log(log(fmsy)-log(fl)), log(log(fu)-log(fmsy)))
                })
                c(log(fmsy),f0)
            }
        }else if(rp$rpType == 3){ ## Max
            Crit <- exp(PRvals[[ss]]$logYPR)
            cutfun <- function(x) rep(TRUE,length(x))
            trans <- function(x, report=FALSE, ...){
                v <- exp(x)
                if(report)
                    names(v) <- "Max"
                v
            }
            fn <- function(x) -predict(CurveFit,trans(x))
            startVals <- function() log(Fseq[which.max(pv)])                        
        }else if(rp$rpType == 4){ ## xdYPR
            stop("Reference point type not implemented yet")
            ## Need derivative!
            Crit <- exp(PRvals[[ss]]$logYPR)
            cutfun <- function(x) rep(TRUE,length(x))
            trans <- function(x, report=FALSE, ...){
                v <- exp(x)
                if(report)
                    names(v) <- paste0(rp$xVal,"dYPR")
                v
            }
            fn <- function(x){
                p0 <- predict(CurveFit,1e-4)
                p <- predict(CurveFit,trans(x))
                sum((p - rp$xVal * p0)^2)
            }
            startVals <- function() sapply(rp$xVal, function(xv) log(Fseq[which.min((pv - xv * pv[1])^2)]))
        }else if(rp$rpType == 5){ ## xSPR
            Crit <- exp(PRvals[[ss]]$logSPR)
            cutfun <- function(x) rep(TRUE,length(x))
            trans <- function(x, report=FALSE, ...){
                v <- exp(x)
                if(report)
                    names(v) <- paste0(rp$xVal,"SPR")
                v
            }
            fn <- function(x){
                p0 <- predict(CurveFit,1e-4)
                p <- predict(CurveFit,trans(x))
                sum((p - rp$xVal * p0)^2)
            }
            startVals <- function() sapply(rp$xVal, function(xv) log(Fseq[which.min((pv - xv * pv[1])^2)]))
        }else if(rp$rpType == 6){ ## xB0
            Crit <- exp(PRvals[[ss]]$logSe)
            cutfun <- function(x) rep(TRUE,length(x))
            trans <- function(x, report=FALSE, ...){
                v <- exp(x)
                if(report)
                    names(v) <- paste0(rp$xVal,"B0")
                v
            }
            fn <- function(x){
                p0 <- predict(CurveFit,1e-4)
                p <- predict(CurveFit,trans(x))
                sum((p - rp$xVal * p0)^2)
            }
            startVals <- function() sapply(rp$xVal, function(xv) log(Fseq[which.min((pv - xv * pv[1])^2)]))
        }else if(rp$rpType == 7){ ## MYPYLdiv
            stop("Reference point type not implemented yet")
            ## Arng <- conf$maxAge - conf$minAge + 1
            ## v <- PRvals$logYe - log(1.0 + exp(PRvals$logYearsLost - log(Arng)))
            ## return(exp(v))
        }else if(rp$rpType == 8){ ## MYPYLprod
            stop("Reference point type not implemented yet")
        }else if(rp$rpType == 9){ ## MDY
            stop("Reference point type not implemented yet")
        }else if(rp$rpType == 10){ ## Crash
            stop("Reference point type not implemented yet")
        }else if(rp$rpType == 11){ ## Ext
            return((exp(PRvals[[ss]]$logSe) - 1)^2)
        }else if(rp$rpType == 12){ ## Lim
            stop("Reference point type not implemented yet")
        }else{
            stop("Reference point type not implemented")
        }
        Frng <- range(Fvals[cutfun(Crit)])
        inRng <- function(x,rng) x > rng[1] & x < rng[2]
        indx <- inRng(Fvals,Frng)
        CurveFit <- .refpointSCurveFit(Fvals[indx], Crit[indx], MT)
        Fseq <- seq(Frng[1],Frng[2],len=200)
        pv <- predict(CurveFit,Fseq)
        opt <- nlminb(startVals(), fn)
        res <- trans(opt$par, report = TRUE)
        attr(res,"curve_fit") <- cbind(F=Fseq,Criteria=pv)
        res
    }
    getDerivedValues <- function(f, ss){
        if(!is.function(MT$derivedSummarizer)){        #Fit
            doOneA <- function(what){
                Crit <- exp(PRvals[[ss]][[what]])
                if(all(is.na(Crit)))
                   return(NA)
                cutfun <- function(x) x > max(x) * 0.1
                Frng <- range(Fvals[cutfun(Crit)])
                inRng <- function(x,rng) x > rng[1] & x < rng[2]
                indx <- inRng(Fvals,Frng)
                CurveFit <- .refpointSCurveFit(Fvals[indx], Crit[indx], MT)
                predict(CurveFit, f)
            }
            return(sapply(c("logYPR","logSPR","logSe","logRe","logYe","logLifeExpectancy","logYearsLost"), doOneA))
        }else if(is.function(MT$derivedSummarizer)){  #Simulate
            return(sapply(lapply(as.list(.perRecruitSR(rep(log(f),nosim),
                                                         fit=fit,
                                                         nYears=nYears,
                                                         aveYears = aveYears,
                                                         selYears = selYears,
                                                         pl = pl,
                                                         ct = catchType,
                                                       nTail = nTail,
                                                       incpb=incpb,
                                                       doSim = doSim,
                                                       label = sprintf("%s derived values",label),
                                                       constraint = constraint,
                                                       deterministicF = deterministicF,
                                                       ...)[[ss]])[c("logYPR","logSPR","logSe","logRe","logYe","logLifeExpectancy","logYearsLost")],
                          exp),
                          function(x){ if(all(is.na(x))) return(NA); MT$derivedSummarizer}))
        }else{
            stop("Derived type not implemented")
        }
    }
    doOneStock <- function(ss){
        Forig <- lapply(rpArgs, getOneRP, ss = ss)
        F <- Reduce("c",(Forig))
        Curves <- lapply(Forig, attr, which = "curve_fit")
        ##names(Curves) <- names(rpArgs)
        D <- sapply(F, getDerivedValues, ss = ss)    
        res <- rbind(logF=unname(F),D)
        colnames(res) <- Reduce("c",lapply(Forig,names))
        rownames(res) <- gsub("^log","",rownames(res))
        ## Fseq <- seq(min(Fvals),max(Fvals),len=200)
        ## GraphVals <- rbind(logF=Fseq,sapply(Fseq, getDerivedValues))
        ## colnames(GraphVals) <- Fseq
        ## rownames(GraphVals) <- gsub("^log","",rownames(GraphVals))
        list(Estimates = res,
             ## GraphVals = GraphVals,
             Fvals = Fvals,
             PRvals = PRvals[[ss]],
             Curves = Curves)
    }
    res <- lapply(seq_along(fit), doOneStock)
    names(res) <- getStockNames(fit)
    res
}




##' @method stochasticReferencepoints msam
##' @importFrom stockassessment stochasticReferencepoints
##' @export
stochasticReferencepoints.msam <- function(fit,
                                          referencepoints,
                                          method = "Q0.5",
                                          catchType = "catch",
                                          nYears = 300,
                                          Frange = c(0,2),
                                          nosim = 1000,
                                          aveYears = lapply(fit,function(x)numeric(0)), #max(fit$data$years)+(-9:0),
                                          selYears = lapply(fit, function(x)max(x$data$years)),
                                          newton.control = list(),
                                          seed = .timeToSeed(),
                                          formula = ~ibc(F,5),
                                          nosim_ci = 200,
                                          derivedSummarizer = NA,
                                          nTail = 1,
                                          constraint = "F=%f",
                                          deterministicF = TRUE,
                                          ...){
    if(nosim_ci > 0){
        warning("Confidence intervals are not implemented yet")
        nosim_ci <- 0
    }
    oldSeed <- NULL
    if(exists(".Random.seed"))
        oldSeed <- .Random.seed
    on.exit(set.seed(oldSeed))
    set.seed(seed)

    .refpointSMethodParser <- getFromNamespace(".refpointSMethodParser","stockassessment")
    .refpointMerger <- getFromNamespace(".refpointMerger","stockassessment")
    .refpointParser <- getFromNamespace(".refpointParser","stockassessment")
    .refpointCheckRecruitment <- getFromNamespace(".refpointCheckRecruitment","stockassessment")
    .refpointSFitCriteria <- getFromNamespace(".refpointSFitCriteria","stockassessment")

    MT <- .refpointSMethodParser(method, formula = formula, derivedSummarizer=derivedSummarizer)

    
    catchType <- pmatch(catchType,c("catch","landing","discard"))-1
    if(is.na(catchType))
        stop("Invalid catch type")

    aveYearsIn <- aveYears
    ## aveYears <- match(aveYears, fit$data$years) - 1
    ## if(any(is.na(aveYears)))
    ##     stop("aveYears has years without data.")

    selYearsIn <- selYears
    ## selYears <- match(selYears, fit$data$years) - 1
    ## if(any(is.na(selYears)))
    ##     stop("selYears has years without data.")

    if(!all(Frange >= 0) && Frange[1] < Frange[2] && length(Frange) ==2)
        stop("Wrong Frange")
    if(!nosim > 0)
        stop("nosim must be a positive integer")

    ## Get RPs for best fit
    rpArgs <- Reduce(.refpointMerger,
                     lapply(referencepoints, .refpointParser, cutoff = 0.1),
                     list())
    invisible(lapply(fit, function(ff) lapply(rpArgs,.refpointCheckRecruitment,fit=ff)))

    doSim <- .getDoSim(logf1=tail(log(fbartable(fit)[,1]),1),
                       fit=fit, nYears = nYears, aveYears = aveYears, selYears = selYears, pl = attr(fit,"m_pl"), constraint=constraint,deterministicF=deterministicF,...)


    pb <- .SAMpb(min = 0, max = nosim * (nosim_ci + 1 + is.function(derivedSummarizer)*length(rpArgs)))
    incpb <- function(label="") .SAM_setPB(pb, pb$getVal()+1,label)
    
    v0 <- .refpointSFitCriteria(rpArgs, pl=fit$pl, MT=MT, fit=fit, nosim=nosim, Frange=Frange, aveYears=aveYears, selYears=selYears, nYears=nYears, catchType=catchType, nTail=nTail,incpb=incpb,doSim=doSim,label="Estimation:",constraint=constraint,deterministicF=deterministicF, ...)

    MakeOne <- function(v0){
        resTabs <- lapply(rownames(v0$Estimates), function(nm){
            tab <- cbind(v0$Estimate[nm,], NA, NA)
            colnames(tab) <- c("Estimate","Low","High")
            rownames(tab) <- colnames(v0$Estimate)
            tab
        })
        names(resTabs) <- rownames(v0$Estimate)
        
        Fseq <- seq(0,2,len=200)
        
        ## Make output tables

        res <- list(tables = list(F = resTabs[["F"]],
                                  Yield = resTabs[["Ye"]],
                                  YieldPerRecruit = resTabs[["YPR"]],
                                  SpawnersPerRecruit = resTabs[["SPR"]],
                                  Biomass = resTabs[["Se"]],
                                  Recruitment = resTabs[["Re"]],
                                  LifeExpectancy = resTabs[["LifeExpectancy"]],
                                  LifeYearsLost = resTabs[["YearsLost"]]
                                  ),
                    graphs = list(F = exp(v0$PRvals$logF),
                                  Yield = exp(v0$PRvals$logYe),
                                  YieldPerRecruit = exp(v0$PRvals$logYPR),
                                  SpawnersPerRecruit = exp(v0$PRvals$logSPR),
                                  Biomass = exp(v0$PRvals$logSe),
                                  Recruitment = exp(v0$PRvals$logRe),
                                  YearsLost = exp(v0$PRvals$logYearsLost),
                                  LifeExpectancy = exp(v0$PRvals$logLifeExpectancy)),
                    regression = v0$Curves,
                    ## opt = NA,
                    ## ssdr = sdr,
                    fbarlabel = substitute(bar(F)[X - Y], list(X = fit$conf$fbarRange[1], Y = fit$conf$fbarRange[2])),
                    stochastic = TRUE
                    ## diagonalCorrection = tv
                    )

        attr(res,"aveYears") <-  aveYearsIn
        attr(res,"selYears") <- selYearsIn

        attr(res,"fit") <- NA
        class(res) <- "sam_referencepoints"

        res
    }
    res <- lapply(v0, MakeOne)
    attr(res,"fit") <- fit
    class(res) <- "msam_referencepoints"
    res

}
