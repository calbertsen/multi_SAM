


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
##' @export
modelforecast.msam <- function(fit,
                          fscale = NULL,
                          catchval = NULL,
                          fval = NULL,
                          nextssb = NULL,
                          landval = NULL,
                          findMSY = NULL,
                          hcr = NULL,
                          nosim = NULL,
                          year.base = lapply(fit,function(x)max(x$data$years)),
                          ave.years = lapply(fit,function(x)max(x$data$years)+(-4:0)),
                          rec.years = lapply(fit,function(x)c()),
                          label = NULL,
                          overwriteSelYears = NULL,
                          deterministicF = FALSE,
                          processNoiseF = TRUE,
                          resampleFirst = !is.null(nosim) && nosim > 0,
                          customSel = NULL,
                          lagR = FALSE,
                          splitLD = FALSE,
                          addTSB = FALSE,
                          biasCorrect = TRUE,
                          returnAllYears = FALSE,
                          useUniroot = FALSE,
                          nCatchAverageYears = rep(1,length(fit)),
                          returnObj = FALSE,
                          hcrConf = lapply(fit, function(x)numeric(0)),
                          ...){
    
    ## Handle year.base < max(fit$data$years)
    ## if(year.base > max(fit$data$years)){
    ##     stop("")
    ## }else if(year.base < max(fit$data$years)){
    ##     warning("year.base is ignored for now")
    ## }
    
    ## Checks
    if(deterministicF && length(fscale) > 0 && any(!is.na(fscale)) && is.null(customSel))
        warning("Forecasted F values may depend on the last estimated F vector and can therefore have uncertainty. Provide a custom selectivity to circumvent this.")

    
    
    nStocks <- length(fit)

    ## Change targets to lists if specified
    targetToList <- function(x){
        if(is.null(x))
            return(replicate(nStocks,numeric(0),FALSE))
        if(is.list(x))
            return(x)
        if(is.vector(x)) ## Same target for all stocks
            x <- matrix(x, nStocks, length(x), byrow = TRUE)
        if(is.matrix(x))
            return(split(x, row(x)))
    }
    
    fscale <- targetToList(fscale)
    catchval <- targetToList(catchval)
    fval <- targetToList(fval)
    nextssb <- targetToList(nextssb)
    landval <- targetToList(landval)
    findMSY <- targetToList(findMSY)
    hcr <- targetToList(hcr)
        
    
    ## Get number of forecast years
    lengthVec <- cbind(sapply(fscale,length),
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
        if(length(fscale[[i]]) == 0)
            fscale[[i]] <- rep(NA_real_, nYears[i])
        if(length(catchval[[i]]) == 0)
            catchval[[i]] <- rep(NA_real_, nYears[i])
        if(length(fval[[i]]) == 0)
            fval[[i]] <- rep(NA_real_, nYears[i])
        if(length(nextssb[[i]]) == 0)
            nextssb[[i]] <- rep(NA_real_, nYears[i])
        if(length(landval[[i]]) == 0)
            landval[[i]] <- rep(NA_real_, nYears[i])
        if(length(findMSY[[i]]) == 0)
            findMSY[[i]] <- rep(NA_real_, nYears[i])
        if(length(hcr[[i]]) == 0)
            hcr[[i]] <- rep(NA_real_, nYears[i])
    }

    
    tabList <- lapply(as.list(1:nStocks),
                      function(i)rbind(fscale[[i]],fval[[i]],catchval[[i]],nextssb[[i]],landval[[i]],findMSY[[i]], hcr[[i]]))
    FModel <- lapply(tabList, function(tab)
        apply(tab,2, function(x){
            y <- which(!is.na(x))
            switch(as.character(length(y)),
                   "0"=0,
                   "1"=y,
                   stop("At most one target can be specified per year."))
        })
        )
    target <- lapply(tabList, function(tab)
        apply(tab,2, function(x){
            y <- which(!is.na(x))
            switch(as.character(length(y)),
                   "0"=NA_real_,
                   "1"=x[y],
                   stop("At most one target can be specified per year."))
        })
        )

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

    args$parameters <- pl
    args$random <- unique(names(obj0$env$par[obj0$env$random]))

    for(i in 1:nStocks)
        args$data$sam[[i]]$forecast <- list(nYears = as.numeric(nYears[i]),
                                            nCatchAverageYears = as.numeric(nCatchAverageYears[i]),
                                            aveYears = as.numeric(ave.years[[i]]),
                                            forecastYear = as.numeric(c(rep(0,fit[[i]]$data$noYears),seq(1,nYears[i],length=nYears[i]))),
                                            FModel = as.numeric(FModel[[i]]),
                                            target = as.numeric(target[[i]]),
                                            selectivity = as.numeric(customSel[[i]]),
                                            recModel = as.numeric(recList[[i]]$recModel),
                                            logRecruitmentMedian = as.numeric(recList[[i]]$logRecruitmentMedian),
                                            logRecruitmentVar = as.numeric(recList[[i]]$logRecruitmentVar),
                                            fsdTimeScaleModel = as.numeric(fsdTimeScaleModel[[i]]),
                                            simFlag = c(0,0),
                                            uniroot = as.numeric(useUniroot),
                                            hcrConf = hcrConf[[i]])
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

    if(returnObj)
        return(obj)

    if(!is.null(nosim) && nosim > 0){


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
        stocksimlist <- list()
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
        attr(stocksimlist,"fit") <- fit
        attr(stocksimlist,"sdreport") <- sdr
        class(stocksimlist) <- "msamforecast"
        return(stocksimlist)
    }
}
    


