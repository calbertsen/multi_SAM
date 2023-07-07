.forecastDefault <- getFromNamespace(".forecastDefault","stockassessment")
.parseRel <- getFromNamespace(".parseRel","stockassessment")
.parseForecast <- getFromNamespace(".parseForecast","stockassessment")
.SAM_replicate <- getFromNamespace(".SAM_replicate","stockassessment")


## Change targets to lists if specified
targetToList <- function(x, nStocks, constructor = numeric){
    if(is.null(x))
        return(replicate(nStocks,constructor(0),FALSE))
    if(is.list(x))
        return(x)
    if(is.vector(x)) ## Same target for all stocks
        x <- matrix(x, nStocks, length(x), byrow = TRUE)
    if(is.matrix(x))
        return(split(x, row(x)))
}



##' Forecast of a multiStockassessment
##' 
##' @param fit 
##' @param fscale 
##' @param catchval 
##' @param fval 
##' @param nextssb 
##' @param landval 
##' @param nosim 
##' @param year.base 
##' @param ave.years 
##' @param rec.years 
##' @param label 
##' @param overwriteSelYears 
##' @param deterministicF 
##' @param processNoiseF 
##' @param customSel 
##' @param lagR 
##' @param splitLD 
##' @param addTSB 
##' @param biasCorrect 
##' @param returnAllYears 
##' @param ...
##' @return msamforecast
##' @importFrom stockassessment modelforecast
##' @method modelforecast msam
##' @export
modelforecast.msam <- function(fit,
                               constraints = NULL,
                          fscale = NULL,
                          catchval = NULL,
                          fval = NULL,
                          nextssb = NULL,
                          landval = NULL,
                          ## findMSY = NULL,
                          ## hcr = NULL,
                          nosim = 1000,
                          year.base = unlist(min(sapply(fit,function(x)max(x$data$years)))),
                          ave.years = lapply(fit,function(x)max(x$data$years)+(-9:0)),
                          overwriteBioModel = FALSE,
                          rec.years = lapply(fit,function(x)numeric(0)),
                          label = NULL,
                          overwriteSelYears = NULL,
                          deterministicF = FALSE,
                          processNoiseF = FALSE,
                          fixedFdeviation = FALSE,
                          useFHessian = FALSE,
                          resampleFirst = !is.null(nosim) && nosim > 0,
                          resampleParameters = FALSE,
                          useModelLastN = TRUE,
                          fixFirstN = FALSE,
                          customSel = NULL,
                          lagR = FALSE,
                          splitLD = FALSE,
                          addTSB = FALSE,
                          biasCorrect = FALSE,
                          returnAllYears = FALSE,
                          ## useUniroot = FALSE,
                          returnObj = FALSE,
                          progress = nosim > 0 && ncores == 1,
                          estimate = median,
                          silent = TRUE,
                          newton_config = NULL,
                          custom_pl = NULL,
                          useNonLinearityCorrection = (nosim > 0 && !deterministicF),
                          ncores = 1,
                          ...){
   
    dots <- list(...)
    findMSY  <-  NULL
    if(!is.na(match("findMSY",names(dots))))
        findMSY <- dots[[match("findMSY",names(dots))]]
    hcr  <-  NULL
    if(!is.na(match("hcr",names(dots))))
        hcr <- dots[[match("hcr",names(dots))]]
    hcrConf = lapply(fit, function(x) numeric(0))
    if(!is.na(match("hcrConf",names(dots))))
        hcrConf <- dots[[match("hcrConf",names(dots))]]
    hcrCurrentSSB = 0
    if(!is.na(match("hcrCurrentSSB",names(dots))))
        hcrCurrentSSB <- dots[[match("hcrCurrentSSB",names(dots))]]
    nCatchAverageYears  <-  rep(1, length(fit))
    if(!is.na(match("nCatchAverageYears",names(dots))))
        nCatchAverageYears <- dots[[match("nCatchAverageYears",names(dots))]]
    fastFixedF <- FALSE
    if(!is.na(match("fastFixedF",names(dots))))
        fastFixedF <- dots[[match("fastFixedF",names(dots))]]
    
   if(!is.null(nosim) && nosim > 0){ 
        estimateLabel <- paste(deparse(substitute(estimate), 500L), collapse = " ")
    }else{
        estimateLabel <- "mostLikelyTrajectory"
    }

    if(is.list(year.base)){
        year.base <- unlist(year.base)
    }
    if(length(year.base) > 1)
        stop("year.base should be a single year to use for all stocks")

    if(year.base > min(sapply(fit, function(x) max(x$data$years))))
        stop("year.base is larger than the last year for some stocks.")
    
    if(progress && !returnObj && !is.null(nosim) && nosim > 0){
        pb <- utils::txtProgressBar(min = 0, max = nosim+3, style = 3)
        incpb <- function() utils::setTxtProgressBar(pb, pb$getVal()+1)
        on.exit(close(pb))
    }else{
        incpb <- function(){ return(invisible(NULL)) }
    }

    
    ## Handle year.base < max(fit$data$years)
    ## if(year.base > max(fit$data$years)){
    ##     stop("")
    ## }else if(year.base < max(fit$data$years)){
    ##     warning("year.base is ignored for now")
    ## }
    
    ## Checks
    if(deterministicF && length(fscale) > 0 && any(!is.na(fscale)) && is.null(customSel))
        warning("Forecasted F values may depend on the last estimated F vector and can therefore have uncertainty. Provide a custom selectivity to circumvent this.")

    if(!is.null(nosim) && nosim > 0){ 
        estimateLabel <- deparse1(substitute(estimate))
    }else{
        estimateLabel <- "mostLikelyTrajectory"
    }
    
    
    nStocks <- length(fit)

    constraints <- targetToList(constraints, nStocks, character)
    fscale <- targetToList(fscale, nStocks)
    catchval <- targetToList(catchval, nStocks)
    fval <- targetToList(fval, nStocks)
    nextssb <- targetToList(nextssb, nStocks)
    landval <- targetToList(landval, nStocks)
    findMSY <- targetToList(findMSY, nStocks)
    hcr <- targetToList(hcr, nStocks)
        
    
    ## Get number of forecast years
    lengthVec <- cbind(sapply(constraints, length),
                       sapply(fscale,length),
                       sapply(catchval,length),
                       sapply(fval, length),
                       sapply(nextssb, length),
                       sapply(landval,length),
                       sapply(findMSY, length),
                       sapply(hcr, length))
    lengthVec <- split(lengthVec, row(lengthVec))
    
    if(any(sapply(lengthVec,function(x) any(x > 0 & x < max(x)))))
        stop("...")    
    nYears <- sapply(lengthVec,max)
    if(any(sapply(seq_along(fit),function(x)max(fit[[x]]$data$years) - year.base >= nYears[[x]])))
        warning("The forecast does not go past the assessment.")

    ## Convert input to an F model code and a target value
    for(i in 1:nStocks){
        if(length(constraints[[i]]) == 0)
            constraints[[i]] <- rep(NA_character_, nYears[i])
        
        if(length(fscale[[i]]) == 0)
            fscale[[i]] <- rep(NA_real_, nYears[i])
        if(any(!is.na(constraints[[i]]) & !is.na(fscale[[i]])))
            warning("fscale specified for years with other constraints. Using previously defined constraints.")
        constraints[[i]][is.na(constraints[[i]]) & !is.na(fscale[[i]])] <- sprintf("F=%f*",fscale[[i]][is.na(constraints[[i]]) & !is.na(fscale[[i]])])
        
        if(length(catchval[[i]]) == 0)
            catchval[[i]] <- rep(NA_real_, nYears[i])
        if(any(!is.na(constraints[[i]]) & !is.na(catchval[[i]])))
        warning("catchval specified for years with other constraints. Using previously defined constraints.")
        constraints[[i]][is.na(constraints[[i]]) & !is.na(catchval[[i]])] <- sprintf("C=%f",catchval[[i]][is.na(constraints[[i]]) & !is.na(catchval[[i]])])
        
        if(length(fval[[i]]) == 0)
            fval[[i]] <- rep(NA_real_, nYears[i])
        if(any(!is.na(constraints[[i]]) & !is.na(fval[[i]])))
        warning("fval specified for years with other constraints. Using previously defined constraints.")
        constraints[[i]][is.na(constraints[[i]]) & !is.na(fval[[i]])] <- sprintf("F=%f",fval[[i]][is.na(constraints[[i]]) & !is.na(fval[[i]])])
        
        if(length(nextssb[[i]]) == 0)
            nextssb[[i]] <- rep(NA_real_, nYears[i])
        if(any(!is.na(constraints[[i]]) & !is.na(nextssb[[i]])))
        warning("nextssb specified for years with other constraints. Using previously defined constraints.")
        constraints[[i]][is.na(constraints[[i]]) & !is.na(nextssb[[i]])] <- sprintf("SSB=%f",nextssb[[i]][is.na(constraints[[i]]) & !is.na(nextssb[[i]])])
        
        if(length(landval[[i]]) == 0)
            landval[[i]] <- rep(NA_real_, nYears[i])
         if(any(!is.na(constraints[[i]]) & !is.na(landval[[i]])))
        warning("landval specified for years with other constraints. Using previously defined constraints.")
        constraints[[i]][is.na(constraints[[i]]) & !is.na(landval[[i]])] <- sprintf("L=%f",landval[[i]][is.na(constraints[[i]]) & !is.na(landval[[i]])])
        
        if(length(findMSY[[i]]) == 0)
            findMSY[[i]] <- rep(NA_real_, nYears[i])
        if(any(!is.na(constraints[[i]]) & !is.na(findMSY[[i]])))
        warning("findMSY specified for years with other constraints. findMSY will overwrite other values.")
        if(length(hcr[[i]]) == 0)
            hcr[[i]] <- rep(NA_real_, nYears[i])
         if(any(!is.na(constraints[[i]]) & !is.na(hcr[[i]])))
        warning("hcr specified for years with other constraints. hcr will overwrite other values.")
    
    }

    
    FModel <- lapply(seq_len(nStocks), function(s){
        v <- rep(0,nYears[s])
        v[!is.na(constraints[[s]])] <- 1
        if(fastFixedF)
            v[!is.na(constraints[[s]])] <- 2
        v[!is.na(findMSY[[s]])] <- 3
        v[!is.na(hcr[[s]])] <- 4
        v
    })

    cstr <- lapply(seq_len(nStocks), function(s){
        v <- replicate(nYears[s], .forecastDefault(), simplify = FALSE)
        v[!is.na(constraints[[s]])] <- .parseForecast(constraints[[s]][!is.na(constraints[[s]])], fit[[s]]$conf$fbarRange, fit[[s]]$data$fleetTypes, c(fit[[s]]$conf$minAge,fit[[s]]$conf$maxAge), useNonLinearityCorrection)
        })
    
    ## Use custom selectivity?
    if(is.null(customSel)){
        customSel <- replicate(nStocks,numeric(0),FALSE)
        FtabList <- faytable(fit, returnList = TRUE)
        if(!is.null(overwriteSelYears)){
            for(i in 1:nStocks){
                fromto <- fit[[i]]$conf$fbarRange-(fit[[i]]$conf$minAge-1)  
                customSel[[i]] <- colMeans(FtabList[[i]][as.integer(rownames(FtabList[[i]])) %in% overwriteSelYears, , drop=FALSE])
                customSel[[i]] <- customSel[[i]]/mean(customSel[[i]][fromto[1]:fromto[2]])
            }
        }
    }else{
        if(!is.null(overwriteSelYears))
            warning("overwriteSelYears is ignored when customSel is given.")
        for(i in 1:nStocks){
            fromto <- fit[[i]]$conf$fbarRange-(fit[[i]]$conf$minAge-1)  
            customSel[[i]] <- customSel[[i]]/mean(customSel[[i]][fromto[1]:fromto[2]])
        }
    }

    ## Get recruitment model
    if(!is.list(rec.years) && is.numeric(rec.years))
        rec.years <- replicate(nStocks,rec.years, simplify=FALSE)
    recList <- replicate(nStocks,list(), FALSE)
    recTabAll <- rectable(fit,returnList = TRUE)
    for(i in 1:nStocks){
        if(length(rec.years[[i]]) == 0){
            recList[[i]]$recModel <- rep(0,nYears[i])
            recList[[i]]$logRecruitmentMedian <- NA_real_
            recList[[i]]$logRecruitmentVar <- NA_real_
        }else{
            rectab <- recTabAll[[i]]
            recpool <- rectab[rownames(rectab)%in%rec.years[[i]],1]
            recList[[i]]$recModel <- rep(1,nYears[i])
            recList[[i]]$logRecruitmentMedian <- median(log(recpool))
            recList[[i]]$logRecruitmentVar <- stats::var(log(recpool))
        }
    }

    ## Get F process time scale model
    ## By default, scale as random walk
    fsdTimeScaleModel <- lapply(1:nStocks,function(i)rep(0,nYears[i]))
    if(deterministicF){ ## 'Zero' variance of F process
        fsdTimeScaleModel <- lapply(1:nStocks,function(i)rep(2,nYears[i]))
    }else if(fixedFdeviation){
        fsdTimeScaleModel <- lapply(1:nStocks,function(i)rep(3,nYears[i]))
    }else if(!processNoiseF){ ## Constant variance of F process
        fsdTimeScaleModel <- lapply(1:nStocks,function(i)rep(1,nYears[i]))
    }
    ## When F model is used, fsdTimeScaleModel should be 1
    for(i in 1:nStocks)
        fsdTimeScaleModel[[i]][FModel[[i]] == 0] <- 1
    
    
    ## Convert average years to indices
    useModelBio <- !overwriteBioModel  ##sapply(ave.years, function(aa) if(length(aa) == 0){ return(TRUE) }else{ return(FALSE)})
    ave.yearsIn <- ave.years
    ## ave.years[useModelBio] <- lapply(fit[useModelBio],function(x) max(x$data$years)+(-4:0))
    ave.years <- lapply(as.list(1:nStocks),function(i)match(ave.years[[i]], fit[[i]]$data$years) - 1)
    if(any(is.na(ave.years)))
        stop("ave.years has years without data.")

    ## Find base year number
    preYears <- sapply(seq_along(fit), function(ii) match(year.base, fit[[ii]]$data$years))
    postYears <- sapply(seq_along(fit), function(ii) nYears[[ii]] - (fit[[ii]]$data$noYears - preYears[[ii]]))

    ## Prepare forecast
    obj0 <- attr(fit,"m_obj")
    args <- as.list(obj0$env)[methods::formalArgs(TMB::MakeADFun)[methods::formalArgs(TMB::MakeADFun) != "..."]]
    args$random <- c("logN", "logF", "missing","logSW", "logCW", "logitMO", "logNM", "logP","logitFseason", "shared_logFscale","shared_missingObs", "logGst", "logGtrip")
    
    pl <- attr(fit,"m_pl")

    if(nosim > 0 && !is.na(match("logF",names(args$map)))){
        ## Remove logF map for simulating
        args$map$logF <- NULL
        pl$logF <- combineParameter(attr(fit,"m_rep")$logFs)
        
    }
    
    logFTmp <- splitMatrices(pl$logF)
    pl$logF <- combineMatrices(lapply(as.list(1:length(logFTmp)),
                                       function(i)cbind(logFTmp[[i]],
                                                        matrix(logFTmp[[i]][,ncol(logFTmp[[i]])],nrow(logFTmp[[i]]),postYears[i]))))
    logNTmp <- splitMatrices(pl$logN)
    pl$logN <- combineMatrices(lapply(as.list(1:length(logNTmp)),
                                       function(i)cbind(logNTmp[[i]],
                                                        matrix(logNTmp[[i]][,ncol(logNTmp[[i]])],nrow(logNTmp[[i]]),postYears[i]))))

    sfsTmp <- splitMatrices(pl$shared_logFscale)
    pl$shared_logFscale <- combineMatrices(lapply(as.list(1:length(sfsTmp)),
                                                  function(i){
                                                      if(ncol(sfsTmp[[i]]) == 0) return(sfsTmp[[i]]);
                                                      cbind(sfsTmp[[i]],
                                                            matrix(sfsTmp[[i]][,ncol(sfsTmp[[i]])],nrow(sfsTmp[[i]]),postYears[i]))
                                                  }))
 
    if(any(useModelBio)){
        splitArray <- function(a){
            nr <- dim(a)[1]; nc <- dim(a)[2]; na <- dim(a)[3]
            lapply(split(a,rep(seq_len(na),each=nr*nc)), matrix, nrow = nr, ncol = nc)
        }
        extendBio <- function(x, useBio) rbind(x,matrix(x[nrow(x)],postYears,ncol(x)))
        makeBio <- function(x, isArray = FALSE){
            if(isArray){
                biox <- split3DArrays(x)
            }else{                
                biox <- splitMatrices(x)
            }
            nr <- sapply(biox,nrow)
            fn <- extendBio
            if(isArray)
                fn <- function(y) simplify2array(lapply(splitArray(y), extendBio))
            biox[nr > 0 & useModelBio] <- lapply(biox[nr > 0 & useModelBio], fn)
            if(isArray)
                return(combine3DArrays(biox))
            return(combineMatrices(biox))
        }
        pl$logitMO <- makeBio(pl$logitMO)
        pl$logNM <- makeBio(pl$logNM)
        pl$logSW <- makeBio(pl$logSW)
        pl$logCW <- makeBio(pl$logCW, TRUE)
    }
    logitFSTmp <- split3DArrays(pl$logitFseason)
    pl$logitFseason <- combine3DArrays(lapply(as.list(1:length(logitFSTmp)),
                                              function(i){
                                                  d0 <- dim(logitFSTmp[[i]]) 
                                                  lfsOld <- logitFSTmp[[i]]
                                                  lfsNew <- array(0, c(d0[1],postYears[i]+d0[2], d0[3]))
                                                  lfsNew[,1:d0[2],] <- lfsOld
                                                  lfsNew
                                                  }))
    
    args$parameters <- pl
    args$random <- unique(names(obj0$env$par[obj0$env$random]))

    if(useFHessian){
        if(all(sapply(fit, function(x) max(x$data$years) == year.base)) ||
           (all(sapply(fit, function(x) max(x$data$years)-1 == year.base)) && useModelLastN)){
            est <- attr(fit,"m_sdrep")$estY
            cov <- attr(fit,"m_sdrep")$covY
        }else if(all(sapply(fit, function(x) max(x$data$years)-1 == year.base))){
            est <- attr(fit,"m_sdrep")$estYm1
            cov <- attr(fit,"m_sdrep")$covYm1
        }else{
            warning("Using last year correlation for useFHessian")
            est <- attr(fit,"m_sdrep")$estY
            cov <- attr(fit,"m_sdrep")$covY
        }
        ## if(year.base==(max(fit$data$years)-1)){
        ##     stop("year.base before last assessment year is not implemented yet")
        ##     est <- attr(fit,"m_sdrep")$estYm1
        ##     cov <- attr(fit,"m_sdrep")$covYm1
        ## }
        nfSplit <- gsub("(^.*[lL]ast)(Log[NF]$)","\\2",names(est))
        stockSplit <- gsub("(^SAM_)([[:digit:]]+)(.+$)","\\2",names(est))
        FEstCov <- vector("list",nStocks)
        for(i in 1:nStocks){
            jj <- nfSplit == "LogF" & as.numeric(stockSplit) == i-1
            FEstCov[[i]] <- cov[jj,jj]
        }
    }else{
        FEstCov <- replicate(nStocks,matrix(0,0,0),simplify=FALSE)
    }

    
    for(i in 1:nStocks)
        args$data$sam[[i]]$forecast <- list(nYears = as.numeric(nYears[i]),
                                            preYears = as.numeric(preYears[i]),
                                            nCatchAverageYears = as.numeric(nCatchAverageYears[i]),
                                            aveYears = as.numeric(ave.years[[i]]),
                                            ##forecastYear = as.numeric(c(rep(0,fit[[i]]$data$noYears),seq(1,nYears[i],length=nYears[i]))),
                                            forecastYear = as.numeric(c(rep(0,preYears[i]),seq(1,nYears[i],length=nYears[i]))),
                                            FModel = as.numeric(FModel[[i]]),
                                            ##target = as.numeric(target[[i]]),
                                            constraints = cstr[[i]],
                                            cfg = newton_config,
                                            selectivity = as.numeric(customSel[[i]]),
                                            recModel = as.numeric(recList[[i]]$recModel),
                                            logRecruitmentMedian = as.numeric(recList[[i]]$logRecruitmentMedian),
                                            logRecruitmentVar = as.numeric(recList[[i]]$logRecruitmentVar),
                                            fsdTimeScaleModel = as.numeric(fsdTimeScaleModel[[i]]),
                                            simFlag = c(0,0),
                                            hcrConf = hcrConf[[i]],
                                            hcrCurrentSSB = hcrCurrentSSB,
                                            Fdeviation = rnorm(nrow(splitMatrices(pl$logF)[[i]])),
                                            FdeviationCov = diag(1,nrow(splitMatrices(pl$logF)[[i]]),nrow(splitMatrices(pl$logF)[[i]])),
                                            FEstCov = FEstCov[[i]],
                                            useModelLastN = useModelLastN
                                            )
    args$data$maxYearAll <- max(unlist(lapply(args$data$sam,function(x)max(x$years) + x$forecast$nYears)))

    withFMSY <- which(sapply(findMSY,function(x)any(!is.na(x))))    
    if(length(withFMSY) > 0){
        mtmp <- rep(NA,length(findMSY))
        mtmp[withFMSY] <- seq_along(withFMSY)
        args$map$logFScaleMSY <- factor(mtmp)
        return(args)
    }

    ## Correct maps for shared parameters
    if(!is.null(attr(fit,"m_call")$shared_selectivity) && attr(fit,"m_call")$shared_selectivity != 0){
        ## lfm0 <- lapply(splitParameter(args$parameters$logF),seq_along)
        ## lfm0[-1] <- lapply(lfm0[-1],function(x) x*NA)
        ## args$map$logF <- factor(unlist(lfm0))
        ## NOT FOR F!
        ## trans_rho should all be the same
        itr0 <- lapply(splitParameter(args$parameters$itrans_rho),seq_along)
        args$map$itrans_rho <- factor(unlist(itr0))
    }
    if(!is.null(attr(fit,"m_call")$shared_seasonality) && attr(fit,"m_call")$shared_seasonality != 0){
        lfm0 <- lapply(splitParameter(args$parameters$logitFseason),seq_along)
        lfm0[-1] <- lapply(lfm0[-1],function(x) x*NA)
        args$map$logitFseason <- factor(unlist(lfm0))
    }

    

    ## Create forecast object
    obj <- do.call(TMB::MakeADFun, args)

    if(as.integer(returnObj)==1)
         return(obj)

    if(!is.null(nosim) && nosim > 0){
        est <- attr(fit,"m_sdrep")$estY
        cov <- attr(fit,"m_sdrep")$covY
        D <- cov2cor(cov)
        Sd <- sqrt(diag(cov))
        yearInsert <- sapply(fit, function(x) max(x$data$years))
        if(diff(range(sapply(fit,function(x)max(x$data$years)))) == 0 &&
           (year.base==max(sapply(fit, function(x) max(x$data$years))) ||
            (year.base==(max(sapply(fit, function(x) max(x$data$years)))-1) && useModelLastN))){
            ## All have same last year and it is equal to year.base
            est <- attr(fit,"m_sdrep")$estY
            cov <- attr(fit,"m_sdrep")$covY
            D <- cov2cor(cov)
            Sd <- sqrt(diag(cov))
            yearInsert <- sapply(fit, function(x) max(x$data$years))
        }else if(diff(range(sapply(fit,function(x)max(x$data$years)))) == 0 &&
                 year.base==(max(sapply(fit, function(x) max(x$data$years)))-1)){
            ## All have sam last year and year.base is year before
            est <- attr(fit,"m_sdrep")$estYm1
            cov <- attr(fit,"m_sdrep")$covYm1
            D <- cov2cor(cov)
            Sd <- sqrt(diag(cov))
            yearInsert <- sapply(fit, function(x) max(x$data$years)-1)
        }else{
            ## Any other scenario
            plsdx <- attr(fitMulti,"m_plsd")
            plx <- attr(fitMulti,"m_pl")
            sdLogN <- splitParameter(plsdx$logN)
            eLogN <- splitParameter(plx$logN)
            sdLogF <- splitParameter(plsdx$logF)
            eLogF <- splitParameter(plx$logF)
            Ngrab <- lapply(seq_len(nStocks), function(s){
                i0 <- match(year.base, fit[[s]]$data$years)
                eLogN[[s]][,i0]
            })
            Fgrab <- lapply(seq_len(nStocks), function(s){
                i0 <- match(year.base, fit[[s]]$data$years)
                eLogF[[s]][,i0]
            })
            sdNgrab <- lapply(seq_len(nStocks), function(s){
                i0 <- match(year.base, fit[[s]]$data$years)
                sdLogN[[s]][,i0]
            })
            sdFgrab <- lapply(seq_len(nStocks), function(s){
                i0 <- match(year.base, fit[[s]]$data$years)
                sdLogF[[s]][,i0]
            })
            est[] <- unname(c(unlist(Ngrab),unlist(Fgrab)))
            Sd[]  <- unname(c(unlist(sdNgrab),unlist(sdFgrab)))
            cov <- diag(Sd^2,length(Sd), length(Sd)) ## Does not use correlation
            yearInsert <- rep(year.base, nStocks)
        }
        ## if(year.base==max(sapply(fit, function(x) max(x$data$years))) ||
        ##    (year.base==(max(sapply(fit, function(x) max(x$data$years)))-1) && useModelLastN)){
        ##     est <- attr(fit,"m_sdrep")$estY
        ##     cov <- attr(fit,"m_sdrep")$covY
        ##     yearInsert <- max(sapply(fit, function(x) max(x$data$years)))
        ## }else if(year.base==(max(sapply(fit, function(x) max(x$data$years)))-1)){
        ##     est <- attr(fit,"m_sdrep")$estYm1
        ##     cov <- attr(fit,"m_sdrep")$covYm1
        ##     yearInsert <- (max(sapply(fit, function(x) max(x$data$years)))-1)
        ## }else{
        ##     stop("")
        ## }
        ## if(year.base==(max(fit$data$years)-1)){
        ##     stop("year.base before last assessment year is not implemented yet")
        ##     est <- attr(fit,"m_sdrep")$estYm1
        ##     cov <- attr(fit,"m_sdrep")$covYm1
        ## }
        nfSplit <- gsub("(^.*[lL]ast)(Log[NF]$)","\\2",names(est))
        stockSplit <- gsub("(^SAM_)([[:digit:]]+)(.+$)","\\2",names(est))
        i0 <- sapply(seq_len(nStocks), function(i) which(fit[[i]]$data$year == yearInsert[i]))
        useMapOnPl <- function(plMap, map){
            with.map <- intersect(names(plMap), names(map))
            applyMap <- function(par.name) {
                tapply(plMap[[par.name]], map[[par.name]], mean)
            }
            plMap[with.map] <- sapply(with.map, applyMap, simplify = FALSE)
            plMap
        }
        map <- obj$env$map ##attr(fit,"m_obj")$env$map
        plMap <- useMapOnPl(pl, map)
        sniii <- 1
        parameterSigma <- svd_solve(attr(fit,"m_opt")$he)
        objX <- attr(fit,"m_obj")
        doSim <- function(re_constraint = NULL, re_pl = NULL){
            obj2 <- obj
            if(!is.null(re_constraint)){
                ## Check length of constraints?
                re_constraint <- targetToList(re_constraint, nStocks, character)
                cstr <- lapply(seq_len(nStocks), function(s){
                    v <- replicate(nYears[s], .forecastDefault(), simplify = FALSE)
                    v[!is.na(re_constraint[[s]])] <- .parseForecast(re_constraint[[s]][!is.na(re_constraint[[s]])], fit[[s]]$conf$fbarRange, fit[[s]]$data$fleetTypes, c(fit[[s]]$conf$minAge,fit[[s]]$conf$maxAge), useNonFALSELinearityCorrection)
                })
                for(i in seq_along(cstr))
                    obj2$env$data$sam[[i]]$forecast$constraints <- cstr[[i]]
            }
            sim0 <- 0*est
            if(resampleFirst)
                sim0 <- rmvnorm(1, mu=0*est, Sigma=cov + diag(1e-6,length(est)))
             ## update N & F before forecast
            dList0 <- split(as.vector(sim0), nfSplit)
            ##estList0 <- split(as.vector(sim0+est), nfSplit)          
            if(!is.null(re_pl)){
                if(resampleParameters)
                    warning("Parameters are not resamples when re_pl is given")
                plMap2 <- useMapOnPl(re_pl, map)
                ## ##map <- obj2$env$map
                ## with.map <- intersect(names(plMap2), names(map))
                ## applyMap <- function(par.name) {
                ##     tapply(plMap2[[par.name]], map[[par.name]], mean)
                ## }
                ## plMap2[with.map] <- sapply(with.map, applyMap, simplify = FALSE)
                p <- unlist(plMap2)
                names(p) <- rep(names(plMap2), times = sapply(plMap2,length))
            }else{
                if(resampleParameters){
                    p0 <- obj2$par
                    pfix <- rmvnorm(1, p0, Sigma = parameterSigma + diag(1e-6,length(p0)))
                    ## Update random effects
                    objX$fn(pfix)
                    ## p <- unlist(plMap)
                    ## names(p) <- rep(names(plMap), times = sapply(plMap,length))
                    ## p[-obj$env$random] <- pfix
                    ## obj2$fn(pfix)
                    p <- objX$env$last.par
                }else{
                    p <- unlist(plMap)
                    names(p) <- rep(names(plMap), times = sapply(plMap,length))
                }
            }
            ## Only works when year.base is last assessment year
            indxN <- local({
                ii <- which(names(p) %in% "logN")
                attr(ii,"cdim") <- attr(obj$env$parameters$logN,"cdim")
                attr(ii,"rdim") <- attr(obj$env$parameters$logN,"rdim")
                indxS <- splitMatrices(ii)
                unlist(sapply(seq_len(nStocks), function(i) indxS[[i]][,i0[i]]))
            })
            indxF <- local({
                ii <- which(names(p) %in% "logF")
                ## if(any(names(map) %in% "logF")){
                ##     ii2 <- attr(obj$env$parameters$logF,"map")+1
                ##     ii2[ii2>0] <- ii[ii2[ii2>0]+1]
                ##     attr(ii2,"cdim") <- attr(attr(obj$env$parameters$logF,"shape"),"cdim")
                ##     attr(ii2,"rdim") <- attr(attr(obj$env$parameters$logF,"shape"),"rdim")
                ##     indxS <- splitMatrices(ii2)
                ## }else{
                attr(ii,"cdim") <- attr(obj$env$parameters$logF,"cdim")
                attr(ii,"rdim") <- attr(obj$env$parameters$logF,"rdim")
                indxS <- splitMatrices(ii)                    
                ## }
                unlist(sapply(seq_len(nStocks), function(i) indxS[[i]][,i0[i]]))
            })
            p[indxN] <- p[indxN] + dList0$LogN
            p[indxF] <- p[indxF] + dList0$LogF[indxF>0]
            fdvAll <- dList0$LogF ## - mean(dList0$LogF)
            for(i in 1:nStocks){
                fdv <- fdvAll[as.numeric(stockSplit[nfSplit == "LogF"]) == (i-1)]
                if(length(fdv) > 0){
                    obj2$env$data$sam[[i]]$forecast$Fdeviation[] <- fdv
                    cindx <- nfSplit == "LogF" & stockSplit == (i-1)
                    obj2$env$data$sam[[i]]$forecast$FdeviationCov <- cov[cindx,cindx]
                }
            }
            v <- obj2$simulate(par = p)
            v$parameterVector <- p
            sniii <<- sniii+1
            incpb()
            return(v)
        }
        eSAM <- getNamespace("stockassessment")
        for(x in ls(eSAM))
            assign(x,get(x, envir = eSAM),envir=environment(doSim))
        eMSAM <- getNamespace("multiStockassessment")
        for(x in ls(eMSAM))
            assign(x,get(x, envir = eMSAM),envir=environment(doSim))
        if(as.integer(returnObj)==2)
            return(doSim)    
        simvals <- .SAM_replicate(nosim, doSim(), simplify = FALSE, ncores = ncores, env = environment(doSim))
        stocksimlist <- vector("list",nStocks)
        for(ss in 1:nStocks){
            fleetHasF <- apply(fit[[ss]]$conf$keyLogFsta>-1,1,any)
            simlist <- vector("list",length(FModel[[ss]]) + 1)
            for(i in 0:(length(FModel[[ss]]))){
                y<-year.base+i
                ii <- which(fit[[ss]]$data$year == year.base) + i
                simlist[[i+1]] <- list(sim = do.call("rbind",lapply(simvals,function(x) c(x$logN[[ss]][,ii], x$logF[[ss]][,ii]))),
                                       fbar = sapply(simvals,function(x) exp(x$logfbar[[ss]][ii])),
                                       catch = sapply(simvals,function(x) exp(x$logCatch[[ss]][ii])),
                                       ssb = sapply(simvals,function(x) exp(x$logssb[[ss]][ii])),
                                       rec = sapply(simvals,function(x) exp(x$logN[[ss]][1,ii])),
                                       cwF = rep(NA_real_, nosim),
                                       catchatage = do.call("cbind",lapply(simvals,function(x) exp(x$logCatchAge[[ss]][,ii]))),
                                       catchbyfleet = do.call("cbind",lapply(simvals,function(x) t(exp(x$logCatchByFleet[[ss]][ii,fleetHasF,drop=FALSE])))),
                                       fbarbyfleet = do.call("cbind",lapply(simvals,function(x) exp(x$logFbarByFleet[[ss]][fleetHasF,ii,drop=FALSE]))),
                                       land = sapply(simvals,function(x) exp(x$logLand[[ss]][ii])),
                                       fbarL = sapply(simvals,function(x) exp(x$logfbarL[[ss]][ii])),
                                       tsb = sapply(simvals,function(x) exp(x$logtsb[[ss]][ii])),
                                       logEmpiricalSPR  = sapply(simvals,function(x) (x$logEmpiricalSPR[[ss]][ii])),
                                       logEmpiricalYPR  = sapply(simvals,function(x) (x$logEmpiricalYPR[[ss]][ii])),
                                       logEmpiricalYPR_L  = sapply(simvals,function(x) (x$logEmpiricalYPR_L[[ss]][ii])),
                                       logEmpiricalYPR_D  = sapply(simvals,function(x) (x$logEmpiricalYPR_D[[ss]][ii])),
                                       year=y)
                rownames(simlist[[i+1]]$catchatage) <- seq(fit[[ss]]$conf$minAge,fit[[ss]]$conf$maxAge,1)
                
            }
            attr(simlist, "fit")<-fit[[ss]]
            ## collect <- function(x){
            ##     quan <- quantile(x, c(.50,.025,.975, na.rm = TRUE))
            ##     c(median=quan[1], low=quan[2], high=quan[3])
            ## }
            collect <- function(x){
                est <- estimate(x, na.rm = TRUE)
                quan <- unname(quantile(x, c(.025,.975), na.rm = TRUE))
                v <- c(estimate=est, low=quan[1], high=quan[2])
                names(v)  <- c(estimateLabel, "low","high")
                v
            }

            fbar <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$fbar))),3)
            fbarL <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$fbarL))),3)  
            rec <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$rec))))
            ssb <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$ssb))))
            tsb <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$tsb))))
            catch <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$catch))))
            land <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$land))))  
            caytable<-round(do.call(rbind, lapply(simlist, function(xx)apply(xx$catchatage,1,collect))))
            cbftable<-round(do.call(rbind, lapply(simlist, function(xx)apply(xx$catchbyfleet,1,collect))))
            fbftable<-round(do.call(rbind, lapply(simlist, function(xx)apply(xx$fbarbyfleet,1,collect))),3)
            colnames(fbftable) <- attr(fit$data,"fleetNames")[fleetHasF]
            colnames(cbftable) <- attr(fit$data,"fleetNames")[fleetHasF]

            tab <- cbind(fbar,rec,ssb,catch)
            if(splitLD){
                tab<-cbind(tab,fbarL,fbar-fbarL,land,catch-land)
            }
            if(addTSB){
                tab<-cbind(tab,tsb)
            }
            ## if(!missing(customWeights)) tab <- cbind(tab,cwF=round(do.call(rbind, lapply(simlist, function(xx)collect(xx$cwF))),3))
            rownames(tab) <- unlist(lapply(simlist, function(xx)xx$year))
            nam <- c(estimateLabel,"low","high")
            basename<-c("fbar:","rec:","ssb:","catch:")
            if(splitLD){
                basename<-c(basename,"fbarL:","fbarD:","Land:","Discard:")    
            }
            if(addTSB){
                basename<-c(basename,"tsb:")    
            }
            ## if(!missing(customWeights)){
            ##     basename<-c(basename,"cwF:")    
            ## }
            colnames(tab)<-paste0(rep(basename, each=length(nam)), nam)
            
            attr(simlist, "tab")<-tab
            shorttab<-t(tab[,grep(estimateLabel,colnames(tab))])
            ##rownames(shorttab)<-sub(sprintf(":%s",estimateLabel),"",paste0(label,if(!is.null(label)){":"},rownames(shorttab)))
            attr(simlist, "shorttab")<-shorttab
            attr(simlist, "label") <- label    
            attr(simlist, "caytable")<-caytable
            attr(simlist, "catchby")<-cbftable
            attr(simlist, "fbarby")<-fbftable
            attr(simlist, "nosim") <- nosim
            class(simlist) <- "samforecast"
            attr(simlist,"estimateLabel") <- estimateLabel
            attr(simlist,"useNonLinearityCorrection") <- useNonLinearityCorrection
            stocksimlist[[ss]] <- simlist
        }
        names(stocksimlist) <- getStockNames(fit)
        attr(stocksimlist,"fit") <- fit
        attr(stocksimlist,"sdreport") <- NA
        class(stocksimlist) <- "msamforecast"
        ## Done with reporting
        incpb()
        if(progress)
            close(pb)        
        return(stocksimlist)
    }else{    
        obj$fn(attr(fit,"m_opt")$par)
        ## Get results
        ADIndex <- obj$env$ADreportIndex()
        IndxKeep <- c("logfbar","logssb","logR","logCatch","logtsb","logLagR","logLand","logDis","loglandfbar","logdisfbar")[c(TRUE,TRUE,TRUE,TRUE,addTSB,lagR,splitLD,splitLD,splitLD,splitLD)]
        stockNames <- paste0("SAM_",0:(nStocks-1))
        ADIndexUse <- list()
        for(i in 1:nStocks){
            tmp <- ADIndex[paste0(stockNames[i],"_",IndxKeep)]
            ADIndexUse <- c(ADIndexUse,lapply(tmp, utils::tail, n = nYears[i] + 1 + (fit[[i]]$data$noYears-1) * as.numeric(returnAllYears)))
        }
        ADNames <- sapply(IndxKeep,function(i) grep(i,names(ADIndex)))
        sdr <- TMB::sdreport(obj, attr(fit,"m_opt")$par, attr(fit,"m_opt")$he,
                             bias.correct= biasCorrect,
                             skip.delta.method = biasCorrect,
                             bias.correct.control = list(sd = TRUE,
                                                         split = ADIndexUse
                                                         )
                             )
        ssdr <- summary(sdr)
        invisible(obj$fn(attr(fit,"m_opt")$par))       
        rp <- obj$report(obj$env$last.par)
        stocksimlist <- vector("list",nStocks)
        for(ss in 1:nStocks){
            simlist <- list()
            for(i in 0:(length(FModel[[ss]]))){
                y<-year.base+i            
                simlist[[i+1]] <- list(sim=NA, fbar=NA, catch=NA, ssb=NA, rec=NA,
                                       cwF=NA, catchatage=NA, land=NA, fbarL=NA, tsb=NA, year=y)
            }
            attr(simlist, "fit")<-fit[[ss]]    

            toCI <- function(x, trans = exp){
                trans(x %*% matrix(c(1,0,1,-2,1,2), nrow = 2, ncol =3))
            }

            indx <- 1:2
            if(biasCorrect)
                indx <- 3:4

            toStockVarName <- function(ss,x)sprintf("SAM_%d_%s",ss-1,x)
            fbar <- toCI(ssdr[rownames(ssdr) %in% toStockVarName(ss,"logfbar"),indx])
            rec <- toCI(ssdr[rownames(ssdr) %in% toStockVarName(ss,"logR"),indx])
            if(lagR)
                rec <- toCI(ssdr[rownames(ssdr) %in% toStockVarName(ss,"logLagR"),indx])
            ssb <- toCI(ssdr[rownames(ssdr) %in% toStockVarName(ss,"logssb"),indx])
            catch <- toCI(ssdr[rownames(ssdr) %in% toStockVarName(ss,"logCatch"),indx])
            ## caytable<-round(do.call(rbind, lapply(simlist, function(xx)apply(xx$catchatage,1,collect))))
            cayvec <- toCI(ssdr[rownames(ssdr) %in% toStockVarName(ss,"logCatchAge"),indx])
            cbfvec <- toCI(ssdr[rownames(ssdr) %in% toStockVarName(ss,"logCatchByFleet"),indx])
            fbfvec <- toCI(ssdr[rownames(ssdr) %in% toStockVarName(ss,"logFbarByFleet"),indx])
                        
            tab <- cbind(fbar,rec,ssb,catch)
            if(splitLD){
                land <- toCI(ssdr[rownames(ssdr) %in% toStockVarName(ss,"logLand"),indx])
                fbarL <- toCI(ssdr[rownames(ssdr) %in% toStockVarName(ss,"loglandfbar"),indx])
                dis <- toCI(ssdr[rownames(ssdr) %in% toStockVarName(ss,"logDis"),indx])
                fbarD <- toCI(ssdr[rownames(ssdr) %in% toStockVarName(ss,"logdisfbar"),indx])
                tab<-cbind(tab,fbarL,fbarD,land,dis)
            }
            if(addTSB){
                tsb <- toCI(ssdr[rownames(ssdr) %in% toStockVarName(ss,"logtsb"),indx])
                tab<-cbind(tab,tsb)
            }
            
            ## Handle cwF !!
            
            ## Row and column names on table
            futureYears <- unlist(lapply(simlist, function(xx)xx$year))

            rownames(tab) <- sort(unique(c(fit[[ss]]$data$years, futureYears)))
            fullTable <- tab
            
            tab <- tab[as.numeric(rownames(tab)) %in% futureYears, , drop = FALSE]
            
            nam <- c(estimateLabel,"low","high")
            basename<-c("fbar:","rec:","ssb:","catch:")
            if(splitLD){
                basename<-c(basename,"fbarL:","fbarD:","Land:","Discard:")    
            }
            if(addTSB){
                basename<-c(basename,"tsb:")    
            }
            ## if(!missing(customWeights)){
            ##     basename<-c(basename,"cwF:")    
            ## }

            if(FALSE){
################# NEEED TO CHECK FOR MULTI-STOCK ############################
                ## Make caycbftable and fbftable
#### cay
                ii <- matrix(seq_along(rp$logCatchAge), nrow(rp$logCatchAge), ncol(rp$logCatchAge))[,as.numeric(rownames(fullTable)) %in% futureYears,drop=FALSE]
                caytable <- do.call("rbind",lapply(seq_len(ncol(ii)), function(yy){
                    v <- t(cayvec[ii[,yy],,drop=FALSE])
                }))
#### cbf
                ii <- matrix(seq_along(rp$logCatchByFleet), nrow(rp$logCatchByFleet), ncol(rp$logCatchByFleet))[as.numeric(rownames(fullTable)) %in% futureYears,fleetHasF,drop=FALSE]
                cbftable <- do.call("rbind",lapply(seq_len(nrow(ii)), function(yy){
                    v <- t(cbfvec[ii[yy,],,drop=FALSE])
                }))
#### fay
                ii <- matrix(seq_along(rp$logFbarByFleet), nrow(rp$logFbarByFleet), ncol(rp$logFbarByFleet))[fleetHasF,as.numeric(rownames(fullTable)) %in% futureYears,drop=FALSE]
                fbftable <- do.call("rbind",lapply(seq_len(ncol(ii)), function(yy){
                    v <- t(fbfvec[ii[,yy],,drop=FALSE])
                }))
#### dimnames
                rownames(caytable) <- rep(c(estimateLabel,"low","high"), length.out = nrow(caytable))
                rownames(cbftable) <- rownames(fbftable) <- rep(c(estimateLabel,"low","high"), length.out = nrow(cbftable))
                colnames(cbftable) <- colnames(fbftable) <- attr(fit$data,"fleetNames")[fleetHasF]
                colnames(caytable) <- seq(fit$conf$minAge,fit$conf$maxAge,1)

            }

            colnames(tab)<-paste0(rep(basename, each=length(nam)), nam)
            if(returnAllYears || !biasCorrect)
                colnames(fullTable) <- colnames(tab)
            attr(simlist, "tab")<-tab
            shorttab<-t(tab[,grep("median",colnames(tab))])
            rownames(shorttab)<-sub(":median","",paste0(label,if(!is.null(label))":",rownames(shorttab)))
            attr(simlist, "shorttab")<-shorttab
            attr(simlist, "label") <- label    
            ## attr(simlist, "caytable")<-caytable
            attr(simlist, "fulltab") <- fullTable
            attr(simlist,"nosim") <- 0    
            class(simlist) <- "samforecast"
            attr(simlist,"estimateLabel") <- estimateLabel
            attr(simlist,"useNonLinearityCorrection") <- useNonLinearityCorrection   
            stocksimlist[[ss]] <- simlist
        }
        names(stocksimlist) <- getStockNames(fit)
        attr(stocksimlist,"fit") <- fit
        attr(stocksimlist,"sdreport") <- sdr
        class(stocksimlist) <- "msamforecast"
        return(stocksimlist)
    }
}
    
