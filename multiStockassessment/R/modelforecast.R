.forecastDefault <- getFromNamespace(".forecastDefault","stockassessment")
.parseRel <- getFromNamespace(".parseRel","stockassessment")
.parseForecast <- getFromNamespace(".parseForecast","stockassessment")


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
                          nosim = 0,
                          year.base = sapply(fit,function(x)max(x$data$years)),
                          ave.years = lapply(fit,function(x)numeric(0)),
                          rec.years = lapply(fit,function(x)numeric(0)),
                          label = NULL,
                          overwriteSelYears = NULL,
                          deterministicF = FALSE,
                          processNoiseF = TRUE,
                          resampleFirst = !is.null(nosim) && nosim > 0,
                          customSel = NULL,
                          lagR = FALSE,
                          splitLD = FALSE,
                          addTSB = FALSE,
                          biasCorrect = FALSE,
                          returnAllYears = FALSE,
                          ## useUniroot = FALSE,
                          returnObj = FALSE,
                          progress = TRUE,
                          estimate = median,
                          silent = TRUE,
                          newton_config = NULL,
                          custom_pl = NULL,
                          useNonLinearityCorrection = (nosim > 0 && !deterministicF),
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
    
   if(!is.null(nosim) && nosim > 0){ 
        estimateLabel <- paste(deparse(substitute(estimate), 500L), collapse = " ")
    }else{
        estimateLabel <- "mostLikelyTrajectory"
    }
 
    
    if(progress && !returnObj && !is.null(nosim) && nosim > 0){
        pb <- utils::txtProgressBar(min = 0, max = nosim+3, style = 3)
        incpb <- function() utils::setTxtProgressBar(pb, pb$getVal()+1)
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
                customSel[[i]] <- colMeans(FtabList[[i]][as.integer(rownames(Ftab)) %in% overwriteSelYears, , drop=FALSE])
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
    }else if(!processNoiseF){ ## Constant variance of F process
        fsdTimeScaleModel <- lapply(1:nStocks,function(i)rep(1,nYears[i]))
    }
    ## When F model is used, fsdTimeScaleModel should be 1
    for(i in 1:nStocks)
        fsdTimeScaleModel[[i]][FModel[[i]] == 0] <- 1

    
    ## Convert average years to indices
    useModelBio <- sapply(ave.years, function(aa) if(length(aa) == 0){ return(TRUE) }else{ return(FALSE)})
    ave.yearsIn <- ave.years
    ave.years[useModelBio] <- lapply(fit[useModelBio],function(x) max(x$data$years)+(-4:0))
    ave.years <- lapply(as.list(1:nStocks),function(i)match(ave.years[[i]], fit[[i]]$data$years) - 1)
    if(any(is.na(ave.years)))
        stop("ave.years has years without data.")
    
    ## Prepare forecast
    obj0 <- attr(fit,"m_obj")
    args <- as.list(obj0$env)[methods::formalArgs(TMB::MakeADFun)[methods::formalArgs(TMB::MakeADFun) != "..."]]
    
    pl <- attr(fit,"m_pl")
    logFTmp <- splitMatrices(pl$logF)
    pl$logF <- combineMatrices(lapply(as.list(1:length(logFTmp)),
                                       function(i)cbind(logFTmp[[i]],
                                                        matrix(logFTmp[[i]][,ncol(logFTmp[[i]])],nrow(logFTmp[[i]]),nYears[i]))))
    logNTmp <- splitMatrices(pl$logN)
    pl$logN <- combineMatrices(lapply(as.list(1:length(logNTmp)),
                                       function(i)cbind(logNTmp[[i]],
                                                        matrix(logNTmp[[i]][,ncol(logNTmp[[i]])],nrow(logNTmp[[i]]),nYears[i]))))
    if(any(useModelBio)){
        splitArray <- function(a){
            nr <- dim(a)[1]; nc <- dim(a)[2]; na <- dim(a)[3]
            lapply(split(a,rep(seq_len(na),each=nr*nc)), matrix, nrow = nr, ncol = nc)
        }
        extendBio <- function(x, useBio) rbind(x,matrix(x[nrow(x)],nYears,ncol(x)))
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
                                                  lfsNew <- array(0, c(d0[1],nYears[i]+d0[2], d0[3]))
                                                  lfsNew[,1:d0[2],] <- lfsOld
                                                  lfsNew
                                                  }))
    
    args$parameters <- pl
    args$random <- unique(names(obj0$env$par[obj0$env$random]))

    for(i in 1:nStocks)
        args$data$sam[[i]]$forecast <- list(nYears = as.numeric(nYears[i]),
                                            nCatchAverageYears = as.numeric(nCatchAverageYears[i]),
                                            aveYears = as.numeric(ave.years[[i]]),
                                            forecastYear = as.numeric(c(rep(0,fit[[i]]$data$noYears),seq(1,nYears[i],length=nYears[i]))),
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
                                            hcrCurrentSSB = hcrCurrentSSB)
    args$data$maxYearAll <- max(unlist(lapply(args$data$sam,function(x)max(x$years) + x$forecast$nYears)))

    withFMSY <- which(sapply(findMSY,function(x)any(!is.na(x))))    
    if(length(withFMSY) > 0){
        mtmp <- rep(NA,length(findMSY))
        mtmp[withFMSY] <- seq_along(withFMSY)
        args$map$logFScaleMSY <- factor(mtmp)
        return(args)
    }


    ## Create forecast object
    obj <- do.call(TMB::MakeADFun, args)

    if(as.integer(returnObj)==1)
         return(obj)

    if(!is.null(nosim) && nosim > 0){
        if(all(year.base==sapply(fit, function(x) max(x$data$years)))){
            est <- attr(fit,"m_sdrep")$estY
            cov <- attr(fit,"m_sdrep")$covY
        }else{
            stop("year.base before last assessment year is not implemented yet")
        }        
        ## if(year.base==(max(fit$data$years)-1)){
        ##     stop("year.base before last assessment year is not implemented yet")
        ##     est <- attr(fit,"m_sdrep")$estYm1
        ##     cov <- attr(fit,"m_sdrep")$covYm1
        ## }
        nfSplit <- gsub("(^.*[lL]ast)(Log[NF]$)","\\2",names(est))
        stockSplit <- gsub("(^SAM_)([[:digit:]]+)(.+$)","\\2",names(est))
        i0 <- sapply(seq_len(nStocks), function(i) which(fit[[i]]$data$year == year.base[i]))
        plMap <- pl
        map <- attr(fit,"m_obj")$env$map
        with.map <- intersect(names(plMap), names(map))
        applyMap <- function(par.name) {
            tapply(plMap[[par.name]], map[[par.name]], mean)
        }
        plMap[with.map] <- sapply(with.map, applyMap, simplify = FALSE)
        sniii <- 1
        doSim <- function(re_constraint = NULL, re_pl = NULL){
            obj2 <- obj
            if(!is.null(re_constraint)){
                ## Check length of constraints?
                re_constraint <- targetToList(re_constraint, nStocks, character)
                cstr <- lapply(seq_len(nStocks), function(s){
                    v <- replicate(nYears[s], .forecastDefault(), simplify = FALSE)
                    v[!is.na(re_constraint[[s]])] <- .parseForecast(re_constraint[[s]][!is.na(re_constraint[[s]])], fit[[s]]$conf$fbarRange, fit[[s]]$data$fleetTypes, c(fit[[s]]$conf$minAge,fit[[s]]$conf$maxAge), useNonLinearityCorrection)
                })
                
                obj2$env$data$forecast$constraints <- cstr
            }
            sim0 <- est
            if(resampleFirst)
                sim0 <- rmvnorm(1, mu=est, Sigma=cov)
            estList0 <- split(as.vector(sim0), nfSplit)
            if(!is.null(re_pl)){
                plMap2 <- re_pl
                map <- fit$obj$env$map
                with.map <- intersect(names(plMap2), names(map))
                applyMap <- function(par.name) {
                    tapply(plMap2[[par.name]], map[[par.name]], mean)
                }
                plMap2[with.map] <- sapply(with.map, applyMap, simplify = FALSE)
                p <- unlist(plMap2)
                names(p) <- rep(names(plMap2), times = sapply(plMap2,length))
            }else{
                p <- unlist(plMap)
                names(p) <- rep(names(plMap), times = sapply(plMap,length))
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
                attr(ii,"cdim") <- attr(obj$env$parameters$logF,"cdim")
                attr(ii,"rdim") <- attr(obj$env$parameters$logF,"rdim")
                indxS <- splitMatrices(ii)
                unlist(sapply(seq_len(nStocks), function(i) indxS[[i]][,i0[i]]))
            })
            p[indxN] <- estList0$LogN
            p[indxF] <- estList0$LogF
            v <- obj2$simulate(par = p)
            sniii <<- sniii+1
            incpb()
            return(v)
        }
        if(as.integer(returnObj)==2)
            return(doSim)
        simvals <- replicate(nosim, doSim(), simplify = FALSE)
        stocksimlist <- vector("list",nStocks)
        for(ss in 1:nStocks){
            simlist <- vector("list",length(FModel[[ss]]) + 1)
            for(i in 0:(length(FModel[[ss]]))){
                y<-year.base[ss]+i
                ii <- i0[ss] + i
                simlist[[i+1]] <- list(sim = do.call("rbind",lapply(simvals,function(x) c(x$logN[[ss]][,ii], x$logF[[ss]][,ii]))),
                                       fbar = sapply(simvals,function(x) exp(x$logfbar[[ss]][ii])),
                                       catch = sapply(simvals,function(x) exp(x$logCatch[[ss]][ii])),
                                       ssb = sapply(simvals,function(x) exp(x$logssb[[ss]][ii])),
                                       rec = sapply(simvals,function(x) exp(x$logN[[ss]][1,ii])),
                                       cwF = rep(NA_real_, nosim),
                                       catchatage = do.call("cbind",lapply(simvals,function(x) exp(x$logCatchAge[[ss]][,ii]))),
                                       land = sapply(simvals,function(x) exp(x$logLand[[ss]][ii])),
                                       fbarL = sapply(simvals,function(x) exp(x$logfbarL[[ss]][ii])),
                                       tsb = sapply(simvals,function(x) exp(x$logtsb[[ss]][ii])),
                                       year=y)
                rownames(simlist[[i+1]]$catchatage) <- seq(fit[[ss]]$conf$minAge,fit[[ss]]$conf$maxAge,1)
                
            }
            attr(simlist, "fit")<-fit[[ss]]
            collect <- function(x){
                quan <- quantile(x, c(.50,.025,.975, na.rm = TRUE))
                c(median=quan[1], low=quan[2], high=quan[3])
            }

            fbar <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$fbar))),3)
            fbarL <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$fbarL))),3)  
            rec <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$rec))))
            ssb <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$ssb))))
            tsb <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$tsb))))
            catch <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$catch))))
            land <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$land))))  
            caytable<-round(do.call(rbind, lapply(simlist, function(xx)apply(xx$catchatage,1,collect)))) 
            tab <- cbind(fbar,rec,ssb,catch)
            if(splitLD){
                tab<-cbind(tab,fbarL,fbar-fbarL,land,catch-land)
            }
            if(addTSB){
                tab<-cbind(tab,tsb)
            }
            ## if(!missing(customWeights)) tab <- cbind(tab,cwF=round(do.call(rbind, lapply(simlist, function(xx)collect(xx$cwF))),3))
            rownames(tab) <- unlist(lapply(simlist, function(xx)xx$year))
            nam <- c("median","low","high")
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
            shorttab<-t(tab[,grep("median",colnames(tab))])
            rownames(shorttab)<-sub(":median","",paste0(label,if(!is.null(label))":",rownames(shorttab)))
            attr(simlist, "shorttab")<-shorttab
            attr(simlist, "label") <- label    
            attr(simlist, "caytable")<-caytable
            class(simlist) <- "samforecast"
            stocksimlist[[ss]] <- simlist
        }
        names(stocksimlist) <- getStockNames(fit)
        attr(stocksimlist,"fit") <- fit
        attr(stocksimlist,"sdreport") <- NA
        class(stocksimlist) <- "msamforecast"
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
        stocksimlist <- vector("list",nStocks)
        for(ss in 1:nStocks){
            simlist <- list()
            for(i in 0:(length(FModel[[ss]]))){
                y<-year.base[[ss]]+i            
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
            
            nam <- c("median","low","high")
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
            if(returnAllYears || !biasCorrect)
                colnames(fullTable) <- colnames(tab)
            attr(simlist, "tab")<-tab
            shorttab<-t(tab[,grep("median",colnames(tab))])
            rownames(shorttab)<-sub(":median","",paste0(label,if(!is.null(label))":",rownames(shorttab)))
            attr(simlist, "shorttab")<-shorttab
            attr(simlist, "label") <- label    
            ## attr(simlist, "caytable")<-caytable
            attr(simlist, "fulltab") <- fullTable
            class(simlist) <- "samforecast"
            stocksimlist[[ss]] <- simlist
        }
        names(stocksimlist) <- getStockNames(fit)
        attr(stocksimlist,"fit") <- fit
        attr(stocksimlist,"sdreport") <- sdr
        class(stocksimlist) <- "msamforecast"
        return(stocksimlist)
    }
}
    


