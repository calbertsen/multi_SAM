## TODO:
## - [ ] Where should genetic data be added?
## - [ ] Extend with custom simulation model (e.g. Kat/NS example)


addSimulatedYears <- function(fit, constraints, ...){
    UseMethod("addSimulatedYears")
}

## Function to add simulated years to fit for single-stock sam
addSimulatedYears.sam <- function(fit, constraints,resampleFirst=FALSE,trueSel=NULL,refit=FALSE,silent=TRUE, resampleParameters = FALSE, ...){

    if(resampleParameters){
        p0 <- rmvnorm(1, fit$opt$par, solve(fit$opt$he), TRUE)
        fit$obj$fn(p0)
        pl0 <- fit$obj$env$parList(par = fit$obj$env$last.par)        
    }else{
        pl0 <- NULL
    }
    
    doSim <- modelforecast(fit, constraints, nosim=1, returnObj=2,addDataYears=TRUE,resampleFirst=resampleFirst, useModelLastN = FALSE,customSel = trueSel, custom_pl = pl0, ...)
    v <- doSim()
    obj <- environment(doSim)$obj
    names(v) <- gsub("dat\\.","",names(v))
    dat <- fit$data
    dat <- dat[!duplicated(names(dat))]
    nms <- intersect(names(dat), names(v))
    dat[nms] <- v[nms]
    nms1 <- c("aux","auxData","idx1","idx2","idxCor","weight")
    dat[nms1] <- obj$env$data[nms1]
    dat$years <- min(as.numeric(dat$aux[,"year"])):max(as.numeric(dat$aux[,"year"]))
    ## Fix dimensionnames
    dmnm <- list(dat$years, fit$conf$minAge:fit$conf$maxAge, NULL)
    dimnames(dat$propMat) <- dimnames(dat$stockMeanWeight) <- dimnames(dat$natMor) <- dimnames(dat$propM) <- dmnm[1:2]
    dimnames(dat$catchMeanWeight) <- dimnames(dat$landFrac) <- dimnames(dat$disMeanWeight) <- dimnames(dat$landMeanWeight) <- dimnames(dat$propF) <- dmnm
    dat$noYears <- length(dat$years)
    cnf <- fit$conf
    pl <- fit$pl
    pl$missing <- NULL
    attr(pl,"what") <- NULL
    nms2 <- intersect(names(pl), names(v))
    ## WARNING: This might break if parameters are reported!
    pl[nms2] <- v[nms2]    
    map <- obj$env$map
    if(refit){        
        ## mpx <- fit$obj$env$map[sapply(fit$obj$env$map,length) > 0]
        ns <- as.list(attr(fit,"call"))$newtonsteps        
        newFit <- suppressWarnings(sam.fit(dat,cnf,pl,map=map,newtonsteps=ifelse(is.null(ns),3,ns),
                                           run=TRUE, check.parameters=FALSE, silent=silent))
    }else{
        newFit <- suppressWarnings(sam.fit(dat,cnf,pl,map=map,
                                           run=FALSE, check.parameters=FALSE))
        newFit$opt <- list(par = newFit$obj$par,
                           objective = NA,
                           convergence = 0)
        plMap <- newFit$pl
        map <- newFit$obj$env$map
        with.map <- intersect(names(plMap), names(map))
        applyMap <- function(par.name) {
            tapply(plMap[[par.name]], map[[par.name]], mean)
        }
        plMap[with.map] <- sapply(with.map, applyMap, simplify = FALSE)         
        p <- unlist(plMap)
        names(p) <- names(newFit$obj$env$last.par)          
        newFit$obj$env$last.par <- newFit$obj$env$last.par.best <- p
        ## SDREPORT
        obj2 <- TMB::MakeADFun(newFit$obj$env$data, newFit$obj$env$parameters, type = "ADFun", 
                               ADreport = TRUE, DLL = newFit$obj$env$DLL, silent = newFit$obj$env$silent)
        
        newFit$rep <- newFit$obj$report(p)
        sdv <- obj2$fn(p)
        sdrep <- list(value = sdv,
                      sd = rep(0,length(sdv)))
        idx <- c(which(names(sdrep$value) == "lastLogN"), which(names(sdrep$value) == 
                                                                "lastLogF"))
        sdrep$estY <- sdrep$value[idx]
        sdrep$covY <- matrix(0,length(idx),length(idx))
        idx <- c(which(names(sdrep$value) == "beforeLastLogN"), which(names(sdrep$value) == 
                                                                      "beforeLastLogF"))
        sdrep$estYm1 <- sdrep$value[idx]
        sdrep$covYm1 <- matrix(0,length(idx),length(idx))
        newFit$sdrep <- sdrep
        class(newFit) <- "sam"
    }
    newFit
}

## Function to add simulated years to fit
## NOT DONE YET
addSimulatedYears.msam <- function(fit, constraints,resampleFirst=FALSE,trueSel=NULL,refit=FALSE,silent=TRUE, resampleParameters = FALSE, ...){

    if(resampleParameters){
        p0 <- rmvnorm(1, attr(fit,"m_opt")$par, solve(attr(fit,"m_opt")$he), TRUE)
        obj <- attr(fit,"m_obj")
        obj$fn(p0)
        pl0 <- obj$env$parList(par = obj$env$last.par)        
    }else{
        pl0 <- NULL
    }
    
    doSim <- modelforecast(fit, constraints, nosim=1, returnObj=2,addDataYears=TRUE,resampleFirst=resampleFirst, useModelLastN = FALSE,customSel = trueSel, custom_pl = pl0, ...)
    v <- doSim()
    obj <- environment(doSim)$obj    
    names(v) <- gsub("dat\\.","",names(v))
    dat <- obj$env$data ## Already updated dimensions by modelforecast
    ## dat <- dat[!duplicated(names(dat))]
    ## dat$sam <- obj$env$data$sam
    for(i in seq_along(fit)){
        ##Update biology + logobs
        nms <- intersect(names(dat$sam[[i]]), names(v))
        for(nn in nms)
            dat$sam[[i]][[nn]] <- v[[nn]][[i]]
        dat$sam[[i]]$years <- min(as.numeric(dat$sam[[i]]$aux[,"year"])):max(as.numeric(dat$sam[[i]]$aux[,"year"]))
        dat$sam[[i]]$noYears <- length(dat$sam[[i]]$years)
        ## Update names of biology
        dmnm <- list(dat$sam[[i]]$years, fit[[i]]$conf$minAge:fit[[i]]$conf$maxAge, dimnames(fit[[i]]$catchMeanWeight)[3])
        dimnames(dat$sam[[i]]$propMat) <- dimnames(dat$sam[[i]]$stockMeanWeight) <- dimnames(dat$sam[[i]]$natMor) <- dimnames(dat$sam[[i]]$propM) <- dmnm[1:2]
        dimnames(dat$sam[[i]]$catchMeanWeight) <- dimnames(dat$sam[[i]]$landFrac) <- dimnames(dat$sam[[i]]$disMeanWeight) <- dimnames(dat$sam[[i]]$landMeanWeight) <- dimnames(dat$sam[[i]]$propF) <- dmnm
        
        ## nms <- intersect(names(dat$sam[[i]]), names(v))
        ## for(nn in nms){
        ##     dat$sam[[i]][[nn]] <- v[[nn]][[i]]
        ## }
        ## nms1 <- c("aux","auxData","idx1","idx2","idxCor","weight")
        ## dat$sam[[i]][nms1] <- obj$env$data$sam[[i]][nms1]
        ## dat$sam[[i]]$years <- min(as.numeric(dat$sam[[i]]$aux[,"year"])):max(as.numeric(dat$sam[[i]]$aux[,"year"]))
        ## ## Fix dimensionnames
        ## dmnm <- list(dat$sam[[i]]$years, fit[[i]]$conf$minAge:fit[[i]]$conf$maxAge, NULL)
        ## dimnames(dat$sam[[i]]$propMat) <- dimnames(dat$sam[[i]]$stockMeanWeight) <- dimnames(dat$sam[[i]]$natMor) <- dimnames(dat$sam[[i]]$propM) <- dmnm[1:2]
        ## dimnames(dat$sam[[i]]$catchMeanWeight) <- dimnames(dat$sam[[i]]$landFrac) <- dimnames(dat$sam[[i]]$disMeanWeight) <- dimnames(dat$sam[[i]]$landMeanWeight) <- dimnames(dat$sam[[i]]$propF) <- dmnm
        ## dat$sam[[i]]$noYears <- length(dat$sam[[i]]$years)
    }
    ## Update shared data
    if(dat$sharedObs$hasSharedObs){        
        dat$sharedObs$logobs <- v$shared_logobs
    }
    ## Update genetics?
    cnf <- lapply(fit,function(x)x$conf)
    pl <- lapply(seq_along(fit),function(i) defpar(dat$sam[[i]],cnf[[i]]))
    ## pl$missing <- NULL
    ## attr(pl,"what") <- NULL
    ## nms2 <- intersect(names(pl), names(v))
    ## ## WARNING: This might break if parameters are reported!
    ## pl[nms2] <- v[nms2]    
    sdat <- lapply(seq_along(fit), function(i){
        nn <- intersect(setdiff(names(dat$sam[[i]]),names(cnf[[i]])),names(fit[[i]]$data))
        dat$sam[[i]][nn]
    })
    mpl <- obj$env$parList(obj$env$last.par.best)
    ## Replace simulated values
    mpl$logN <- combineParameter(v$logN)
    mpl$logF <- combineParameter(v$logF)
    mpl$logitFseason <- combineParameter(v$logitFseason)
    ## BIOPAR
    map <- lapply(fit,function(x) x$obj$env$map)
    mmap <- obj$env$map
    if(refit){        
        ## mpx <- fit$obj$env$map[sapply(fit$obj$env$map,length) > 0]
        ns <- lapply(fit,function(x) as.list(attr(x,"call"))$newtonsteps)
        newFitS <- do.call(c,lapply(seq_along(fit), function(i) suppressWarnings(sam.fit(sdat[[i]],cnf[[i]],pl[[i]],map=map[[i]],newtonsteps=ifelse(is.null(ns[[i]]),3,ns[[i]]),
                                                                                         run=TRUE, check.parameters=FALSE, silent=silent))))
        names(newFitS) <- getStockNames(fit)
        ## ADD OTHER OPTIONS!
        ee <- attr(fit,"m_envir")
        cc <- as.list(attr(fit,"m_call"))[-1]
        ee$newFitS <- newFitS
        ee$sharedObs <- dat$sharedObs
        cc$x <- newFitS
        cc$shared_data <- sharedObs
        newFit <- do.call(multisam.fit,cc, envir = ee)
    }else{
        newFitS <- do.call("c",lapply(seq_along(sdat), function(i){
            ##dd <- sdat[[snm]]
            ## lo <- dd$logobs
            ## lo[is.na(lo)] <- 0
            ## dd$logobs[] <- NA
            suppressWarnings(f <- sam.fit(sdat[[i]], cnf[[i]], pl[[i]],map=map[[i]], run = FALSE, pre.clean=FALSE))
            f$sdrep <- list(value = numeric(0))
            f$opt <- list(objective = NA,convergence=0)
            ##f$pl$missing[] <- lo
            class(f) <- "sam"
            f
        }))
        names(newFitS) <- getStockNames(fit)
                
        ## ADD OTHER OPTIONS!
        ee <- attr(fit,"m_envir")
        cc <- as.list(attr(fit,"m_call"))[-1]
        ee$newFitS <- newFitS
        ee$sharedObs <- dat$sharedObs
        cc$x <- newFitS
        cc$shared_data <- sharedObs
        ##cc$map <- mmap ## NOT IMPLEMENTED??
        cc$parlist <- mpl
        cc$run <- FALSE
        obj <- do.call(multisam.fit,cc, envir = ee)
        newFitM <- newFitS
        class(newFitM) <- "msam"
        attr(newFitM,"m_opt") <- list(par = obj$par,
                                      objective = NA,
                                      convergence = 0,
                                      he = attr(fit,"m_opt")$he)
        ## pl
        attr(newFitM,"m_pl") <- mpl
        ## sdrep
        ## plMap <- mpl
        ## map <- obj$env$map
        ## with.map <- intersect(names(plMap), names(map))
        ## applyMap <- function(par.name) {
        ##     tapply(plMap[[par.name]], map[[par.name]], mean)
        ## }
        ## plMap[with.map] <- sapply(with.map, applyMap, simplify = FALSE)         
        ## p <- unlist(plMap)
        ## names(p) <- names(obj$env$last.par)          
        ## obj$env$last.par <- obj$env$last.par.best <- p
        attr(newFitM,"m_obj") <- obj        
         ## rep
        attr(newFitM,"m_rep") <- obj$report(obj$env$last.par.best)
        ## SDREPORT
        obj2 <- TMB::MakeADFun(obj$env$data, obj$env$parameters, type = "ADFun", 
                               ADreport = TRUE, DLL = obj$env$DLL, silent = obj$env$silent)
        sdv <- obj2$fn(obj$env$last.par.best[obj$env$lfixed()])
        sdrep <- list(value = sdv,
                      sd = rep(0,length(sdv)))
        
        idxL <- c(grep("_lastLogN$",names(sdrep$value)), grep("_lastLogF$",names(sdrep$value)),
                 grep("_lastLogSW$",names(sdrep$value)),grep("_lastLogCW$",names(sdrep$value)),
                 grep("_lastLogitMO$",names(sdrep$value)),grep("_lastLogNM$",names(sdrep$value)))
        sdrep$estY <- sdrep$value[idxL]
        sdrep$covY <- attr(fit,"m_sdrep")$covY

        idxBL <- c(grep("_beforeLastLogN$",names(sdrep$value)), grep("_beforeLastLogF$",names(sdrep$value)),
                 grep("_beforeLastLogSW$",names(sdrep$value)),grep("_beforeLastLogCW$",names(sdrep$value)),
                 grep("_beforeLastLogitMO$",names(sdrep$value)),grep("_beforeLastLogNM$",names(sdrep$value)))
        sdrep$estYm1 <- sdrep$value[idxBL]
        sdrep$covYm1 <- attr(fit,"m_sdrep")$covYm1 ##diag(0,length(sdrep$estYm1))

        sdrep$estYYm1 <- sdrep$value[c(idxL,idxBL)]
        sdrep$covYYm1 <- attr(fit,"m_sdrep")$covYYm1##diag(0,length(sdrep$estYYm1))

        ## covRecPars
        ## idx <- grep("_rec_pars$",names(sdrep$value))
        sdrep$covRecPars <- attr(fit,"m_sdrep")$covRecPars #diag(0,length(idx))
        ## colnames(sdrep$covRecPars) <- rownames(sdrep$covRecPars) <- names(sdrep$value)[idx]
        attr(newFitM,"m_sdrep") <- sdrep
        
        ## plsd
        attr(newFitM,"m_plsd") <- lapply(mpl,function(x) 0*x)
        ## data
        attr(newFitM,"m_data") <- obj$env$data
        ##low
        attr(newFitM,"m_low") <- attr(fit,"m_low")
        ##high
        attr(newFitM,"m_high") <- attr(fit,"m_high")
        ##corStructure
        attr(newFitM,"corStructure") <- attr(fit,"corStructure")  
        ##partialCors
        attr(newFitM,"partialCors") <- attr(fit,"partialCors")
        ## lower
        attr(newFitM,"lower") <- attr(fit,"lower")
        ## upper
        attr(newFitM,"upper") <- attr(fit,"upper")
        ## newtonsteps
        attr(newFitM,"newtonsteps") <- attr(fit,"newtonsteps")
        ## nlminb.control
        attr(newFitM,"nlminb.control") <- attr(fit,"nlminb.control")
        ## dotargs
        attr(newFitM,"dotargs") <- attr(fit,"dotargs")
        ## call
        attr(newFitM,"m_call") <- attr(fit,"m_call")
        ## envir
        attr(newFitM,"m_envir") <- attr(fit,"m_envir")
        ## corParameters
        attr(newFitM,"corParameters") <- NA

    }
    newFitM
}


reduce_shared <- function(data, year = NULL, fleet = NULL, age = NULL){
    yOrig <- min(as.numeric(data$aux[, "year"])):max(as.numeric(data$aux[,"year"]))
    aOrig <- min(as.numeric(data$aux[data$aux[, "age"] >= 0,"age"])):max(as.numeric(data$aux[, "age"]))
    nam <- c("year", "fleet", "age")[c(length(year) > 0, length(fleet) > 0, length(age) > 0)]
    if ((length(year) == 0) & (length(fleet) == 0) & (length(age) == 0)) {
        idx <- rep(TRUE, nrow(data$aux))
    }
    else {
        idx <- !do.call(paste, as.data.frame(data$aux[, nam, drop = FALSE])) %in% do.call(paste, as.data.frame(cbind(year = year, fleet = fleet, age = age)))
    }
    data$aux <- data$aux[idx, ]
    data$auxData <- data$auxData[idx, ]
    data$logobs <- data$logobs[idx]
    data$weight <- data$weight[idx]
    suf <- sort(unique(data$aux[, "fleet"]))
    data$noFleets <- length(suf)
    data$fleetTypes <- data$fleetTypes[suf]
    oldYears <- data$years
    data$years <- min(as.numeric(data$aux[, "year"])):max(as.numeric(data$aux[,"year"]))
    ages <- min(as.numeric(data$aux[data$aux[, "age"] >= 0, "age"])):max(as.numeric(data$aux[,"age"]))
    data$noYears <- length(data$years)
    mmfun <- function(f, y, ff) {
        idx <- which(data$aux[, "year"] == y & data$aux[, "fleet"] == 
                     f)
        ifelse(length(idx) == 0, NA, ff(idx) - 1)
    }
    data$idx1 <- outer(suf, data$years, Vectorize(mmfun, c("f","y")), ff = min)
    data$idx2 <- outer(suf, data$years, Vectorize(mmfun, c("f","y")), ff = max)
    data$idxCor <- data$idxCor[suf, match(data$years, oldYears)]
    data$nobs <- length(data$logobs[idx])
    data
}



updateAssessment <- function(OM, EM, knotRange, AdviceLag, intermediateFleets){
    if(is(OM,"sam")){
        if(is(EM,"msam"))
            stop("The code cannot handle a single-stock OM and multi-stock EM at the moment.")
        datNew <- OM$data
        confNew <- EM$conf
        if(AdviceLag > 0){
            if(length(intermediateFleets) == 0){ ## Reduce fully
                datNew <- reduce(datNew,
                                 year = tail(datNew$years, AdviceLag),
                                 conf = confNew)
                confNew <- attr(datNew,"conf")
            }else{ ## Keep some fleets in first year
                ## First remove intermediate years 2+
                if(AdviceLag > 1){
                    datNew <- reduce(datNew,
                                     year = tail(datNew$years, AdviceLag-1),
                                     conf = confNew)
                    confNew <- attr(datNew,"conf")
                }
                ## Remove fleets not in intermediateFleets for first intermediate year
                if(is.character(intermediateFleets)){
                    intermediateFleets <- match(intermediateFleets, attr(datNew,"fleetNames"))
                    if(any(is.na(intermediateFleets)))
                        stop("intermediateFleets names does not match model fleet names")
                }
                datNew <- reduce(datNew,
                                 year = tail(datNew$years, 1),
                                 fleet = setdiff(seq_along(datNew$fleetTypes), intermediateFleets),
                                 conf = confNew)
                confNew <- attr(datNew,"conf")            
            }
        }
        ## Prepare parameters
        dp <- defpar(datNew, confNew)
        pl <- EM$pl
        for (nn in intersect(names(dp), names(pl)))
            if (length(dp[[nn]]) == length(pl[[nn]])){
                dp[[nn]][] <- pl[[nn]][]
            }## else if(nn %in% c("logN","logF")){
        ##     nnIndx <- pmin(ncol(dp[[nn]]),ncol(pl[[nn]]))
        ##     dp[[nn]][,seq_len(nnIndx)] <- pl[[nn]][,seq_len(nnIndx),drop=FALSE]
        ##     if(ncol(pl[[nn]]) < ncol(dp[[nn]]))
        ##         dp[[nn]][,-seq_len(nnIndx)] <- do.call("cbind", replicate(ncol(dp[[nn]]) - ncol(pl[[nn]]),pl[[nn]][,ncol(pl[[nn]]),drop=FALSE],simplify=FALSE))
        ## }

        if(length(knotRange) > 0 &&
           (EM$conf$stockRecruitmentModelCode == 90 ||
            EM$conf$stockRecruitmentModelCode == 91 ||
            EM$conf$stockRecruitmentModelCode == 92)){
            capture.output(EM_New <- getSplineRecBreaks(datNew, confNew, dp, srmc=EM$conf$stockRecruitmentModelCode, knotRange = knotRange, returnFit=TRUE, sim.condRE=c(TRUE,FALSE), silent = TRUE, newtonsteps=0))
        }else{
            capture.output(EM_New <- sam.fit(datNew, confNew, dp, silent = TRUE, newtonsteps=0))
        }
        return(EM_New)
    }else if(is(OM,"msam")){
        datNew <- attr(OM,"m_data")
        ## Handle advice year lag        
        if(AdviceLag > 0){
            if(length(intermediateFleets) == 0){ ## Reduce fully
                ## Reduce stock data
                for(i in seq_along(datNew$sam)){
                    datNew$sam[[i]] <- reduce(datNew$sam[[i]],
                                              year = tail(datNew$sam[[i]]$years, AdviceLag),
                                              conf = confNew[[i]])
                    confNew[[i]] <- attr(datNew[[i]],"conf")
                }
                ## Reduce shared data
                if(datNew$sharedObs$hasSharedObs){
                    
                }
                ## Reduce genetic data
            }else{ ## Keep some fleets in first year
                ## First remove intermediate years 2+
                if(AdviceLag > 1){
                    for(i in seq_along(datNew$sam)){
                        datNew$sam[[i]] <- reduce(datNew$sam[[i]],
                                                  year = tail(datNew$sam[[i]]$years, AdviceLag-1),
                                                  conf = confNew[[i]])
                        confNew[[i]] <- attr(datNew[[i]],"conf")
                    }
                    ## Reduce shared data
                    if(datNew$sharedObs$hasSharedObs){

                    }
                }
                ## Remove fleets not in intermediateFleets for first intermediate year
                if(is.character(intermediateFleets)){
                    intermediateFleets <- match(intermediateFleets, attr(datNew,"fleetNames"))
                    if(any(is.na(intermediateFleets)))
                        stop("intermediateFleets names does not match model fleet names")
                }
                for(i in seq_along(datNew$sam)){
                    datNew$sam[[i]] <- reduce(datNew$sam[[i]],
                                              year = tail(datNew$sam[[i]]$years, 1),
                                              fleet = setdiff(seq_along(datNew$fleetTypes), intermediateFleets),
                                              conf = confNew[[i]])
                    confNew[[i]] <- attr(datNew[[i]],"conf")
                }
                ## Reduce shared data
                if(datNew$sharedObs$hasSharedObs){

                }
            }
        }
        if(is(EM,"sam")){
            ## Collect data

            ## Handle advice year lag

            ## Setup new fit
        }else if(is(EM,"msam")){
            confNew <- lapply(EM,function(x) x$conf)
            ## Prepare data
            ## 1) Prep single stock
            sdat <- lapply(seq_along(EM), function(i){
                nn <- intersect(setdiff(names(datNew$sam[[i]]),names(confNew[[i]])),names(EM[[i]]$data))
                datNew$sam[[i]][nn]
            })
            ## 2) Prep shared data
            ## Setup new fit
            ## 1) Fake single stock
            newFitS <- do.call("c",lapply(seq_along(sdat), function(i){
                dd <- sdat[[i]]
                isMis <- is.na(EM[[i]]$data$logobs)
                ## Set missing obs
                if(all(isMis)){
                    dd$logobs[] <- NA
                }else{
                    ## NOTE: This works because newer years are at the end
                    ## => if #5 is missing with data from 1960-2024, #5 is also missing with data from 1960-1925
                    dd$logobs[which(isMis)] <- NA
                }
                pl <- defpar(dd, confNew[[i]])
                suppressWarnings(f <- sam.fit(dd, confNew[[i]], pl[[i]],map=map[[i]], run = FALSE, pre.clean=FALSE))
                f$sdrep <- list(value = numeric(0))
                f$opt <- list(objective = NA,convergence=0)
                ##f$pl$missing[] <- lo
                class(f) <- "sam"
                f
            }))
            names(newFitS) <- getStockNames(EM)
            ## 2) Run multistock
            ee <- attr(fit,"m_envir")
            cc <- as.list(attr(fit,"m_call"))[-1]
            ee$newFitS <- newFitS
            ee$sharedObs <- datNew$sharedObs
            cc$x <- as.name("newFitS")
            cc$shared_data <- as.name("sharedObs")
            EM_New <- do.call(multisam.fit,cc, envir = ee)
        }
        return(EM_New)
    }
    ## We should not reach this point
    NULL
}

## TODO (from VT):
## - Add benchmark process: Refit a number of scenarios and choose automatically, refit reference points, etc.
## - Ways to simulate changes in selectivity
## - Implement restriction on the change in TAC (e.g., -20% to +25%)
## - Implement escapement clause on TAC change restriction (e.g., only if SSB > Btrigger)

## Management plans for inspiration:
## Blue Whiting: https://d3b1dqw2kzexi.cloudfront.net/media/8742/agreed-record-blue-whiting-2017.pdf
## Celtic herring: https://www.ices.dk/sites/pub/Publication%20Reports/Advice/2018/Special_requests/eu.2018.03.pdf


##' Management strategy evaluation using SAM models
##'
##' @param OM sam.fit that will work as operating model
##' @param EM sam.fit that will work as estimation model
##' @param nYears Number of years to run simulation
##' @param AdviceForecastSettings Settings to do forecast that determines advice
##' @param AdviceYears Number of years advice given at a time. How advice is given is determined by AdviceForecastSettings
##' @param AdviceLag Lag between assessment and advice 
##' @param initialAdvice Advice in the first AdviceLag years
##' @param implementationError Function to add implementation error (i.e, transform advice to target catch)
##' @param knotRange Range of spline knot values to try
##' @param intermediateFleets Fleets that are available in the (first) intermediate year
##' @param OMselectivityFixed Fix selectivity in OM?
##' @param ... arguments passed on to addSimulatedYears
##' @return a list with MSE result
MSE <- function(OM,
                EM,
                nYears,
                AdviceForecastSettings,
                AdviceYears = 1,
                AdviceLag = 0,
                initialAdvice = NA,
                implementationError = function(x) x + 1e-4,
                knotRange = 3,
                intermediateFleets = numeric(0),
                OMselectivityFixed = FALSE,
                singleStockManagementToPopulation = function(x,OM){
                    v <- tail(catchtable(OM),1)[(seq_len(length(OM)*3)-1)%%3 == 0]
                    sum(x) * v / sum(v)
                },
                adviceToManagement = function(x) x,
                adviceMethod = c("Basic","ICES","Empirical"),
                EMReferencePoints = NULL,
                ...){

    adviceMethod <- match.arg(adviceMethod)
    if(adviceMethod == "ICES" && is.null(EMReferencePoints) && !is.list(EMReferencePoints)){
        stop("A list of reference points Ftarget, Btrigger and Blimit must be given for ICES-like advice")
    }else if(adviceMethod == "ICES" && is.list(EMReferencePoints)){
        if(is(EM,"msam")){
            if(length(EM) != length(EMReferencePoints))
                stop("EMReferencePoints should have the same number of stocks as the EM")
            if(any(sapply(EMReferencePoints,function(x) any(is.na(match(c("Ftarget","Btrigger","Blimit"),names(x)))))))
                stop("In EMReferencePoints, all stocks should have Ftarget, Btrigger, and Blimit")
        }else{
            if(any(is.na(match(c("Ftarget","Btrigger","Blimit"),names(x)))))
                stop("For a single stock EM, EMReferencePoints should have Ftarget, Btrigger, and Blimit")
        }
    }
    
    ## Assume that OM is a multi stock model
    if(!is(OM,"msam"))
        stop("This is for multiStock OM")
    nStocksOM <- length(OM)
    #nStocksEM <- ifelse(is(EM,"msam"),length(EM),1)
##### Check input #####
    ## Operating model should have catch in the final year (i.e., no intermediate year)
    catchInFinalYear <- sapply(OM,function(x) max(x$data$aux[x$data$aux[,"fleet"] %in% which(x$data$fleetTypes==0),"year"]) == max(x$data$years))
    if(!all(catchInFinalYear))
        stop("Operating model must have catches in the final year")
    
    ## AdviceLag should be non-negative
    AdviceLag <- pmax(0, AdviceLag)

    ## Selectivity in operating model
    if(OMselectivityFixed){
        OM_pl <- OM$pl
        trueSel <- as.numeric(exp(t(OM_pl$logF[OM$conf$keyLogFsta[1,]+1,ncol(OM_pl$logF)])) / tail(fbartable(OM)[,1],1))
    }else{
        trueSel <- NULL
    }

##### Make sure number of years match AdviceYears #####
    nYOld <- nYears
    nYears <- max(seq(1,nYears + AdviceYears-1, by = AdviceYears)) + (AdviceYears - 1)
    if(nYOld < nYears)
        message(sprintf("nYears changed to %d to fit the AdviceYears increments.",nYears))

##### Prepare for output #####
    ssb <- array(NA,c(nYears+AdviceLag,nStocksOM+1,5))
    fbar <- array(NA,c(nYears+AdviceLag,nStocksOM+1,5))
    rec <- array(NA,c(nYears+AdviceLag,nStocksOM+1,5))
    catch <- array(NA,c(nYears+AdviceLag,nStocksOM+1,5))

    
    dimnames(ssb) <- dimnames(fbar) <- dimnames(rec) <- dimnames(catch) <- list(seq(do.call(max,lapply(OM,function(x)x$data$years)) + 1,len = nYears + AdviceLag),
                                                                                c(getStockNames(OM),"Total"),
                                                                                c("Advice","True","Estimate","Low","High"))

##### Copy OM and EM for safe overwriting #####
    OM_update <- OM
    EM_update <- EM
    msg <- "OK"

##### Insert initialAdvice in results table #####
    if(is(EM,"msam")){
        for(s in seq_along(fit)){
            catch[seq_len(AdviceLag),s,"Advice"] <- rep(initialAdvice[[s]],length.out = AdviceLag)
        }
    }else{
        v <- singleStockAdviceToSubstock(rep(Reduce("+",initialAdvice),length.out = AdviceLag),OM)
        for(s in seq_along(fit)){
            catch[seq_len(AdviceLag),s,"Advice"] <- v[s,]
        }
    }
    ## Total
    catch[seq_len(AdviceLag),"Total","Advice"] <- rep(Reduce("+",initialAdvice),length.out = AdviceLag)

##### Helper function to convert advice (number) to forecast constraint for addSimulatedYears #####
    AdviceToCatchConstraint <- function(x, OM){
        lapply(x, function(y) ifelse(is.na(y),NA_real_,sprintf("C=%f",implementationError(singleStockAdviceToPopulation(y,OM)))))
    }

###### Make sure AdviceForecast is long enough #####
    ## Default: Need to forecast AdviceLag years + AdviceYears years
    yx <- 0
    ## If year.base is set, we need to check if we need another year
    if(!is.null(AdviceForecastSettings$year.base)){
        ## We can't use a specific year as year.base
        if(is.numeric(AdviceForecastSettings$year.base))
            stop("Advice forecast year.base cannot be set to a specific year in the MSE. Use 'lastCatchYear', 'secondLastYear', or 'lastYear'.")
        ## If year.base is last catch year, and we have intermediate year fleets, we need to forecast one more year
        if(AdviceForecastSettings$year.base == "lastCatchYear" &&
           length(AdviceForecastSettings$intermediateFleets) > 0){            
            yx <- 1
            ## Else, if year.base is secondLastYear, we need to forecast one more year
        }else if(AdviceForecastSettings$year.base == "secondLastYear"){
            yx <- 1
        }
        ## Otherwise, we are good.
    }
    ## Check AdviceForecast is long ennough
    if(any(sapply(AdviceForecastSettings$constraints,length) < AdviceLag + AdviceYears + yx)){
        warning("Length of AdviceForecastSettings$constraints should equal AdviceLag + AdviceYears + 1. Modifying to match.")
        ## If it's not set, insert NA to forecast using model
        if(is.null(AdviceForecastSettings$constraints)){
            if(is(EM,"msam")){
                AdviceForecastSettings$constraints <- replicate(nStocksOM,rep(NA, length.out = AdviceLag + AdviceYears + yx), simplify = FALSE)
            }else{
                AdviceForecastSettings$constraints <- rep(NA, length.out = AdviceLag + AdviceYears + yx)
            }
        }else{
            ## Otherwise, repeat the last constraint to make long enough
            if(is(EM,"msam")){
                AdviceForecastSettings$constraints <- lapply(AdviceForecastSettings$constraints,
                                                             function(x) c(x,
                                                                           rep(tail(x,1),
                                                                               length.out = AdviceLag + AdviceYears + yx - length(x))))
            }else{
                AdviceForecastSettings$constraints <- c(AdviceForecastSettings$constraints,
                                                        rep(tail(AdviceForecastSettings$constraints,1),
                                                            length.out = AdviceLag + AdviceYears + yx - length(AdviceForecastSettings$constraints)))
            }
        }
    }

    if(AdviceLag > 0){
        capture.output(OM_update <- try({addSimulatedYears(OM_update, 
                                                           constraints = AdviceToCatchConstraint(initialAdvice,OM_update),
                                                          ...)}))
        iy <- head(rownames(ssb),AdviceLag)
        ssb[iy,,"True"] <- ssbtable(OM_update,addTotal=TRUE)[iy,seq(1,by=3,length.out=nStocksOM+1)]
        fbar[iy,,"True"] <- fbartable(OM_update,addTotal=TRUE)[iy,seq(1,by=3,length.out=nStocksOM+1)]
        rec[iy,,"True"] <- rectable(OM_update,addTotal=TRUE)[iy,seq(1,by=3,length.out=nStocksOM+1)]
        catch[iy,,"True"] <- catchtable(OM_update,addTotal=TRUE)[iy,seq(1,by=3,length.out=nStocksOM+1)]
    }
    ## Run assessment loop
    for(i in seq(1,nYears-(AdviceYears-1), by = AdviceYears)){ # Index over assessment year
        yr <- rownames(ssb)[seq(i,len=AdviceYears)]
        yr_tac <- rownames(ssb)[seq(i+AdviceLag,len=AdviceYears)]
        
        cat("\n\n\nSimulation year",i,"\n")
        cat("\tAdvice year",yr_tac,"\n")
        cat("\tLast true data year",tail(rownames(fbartable(OM_update)),1),"\n")
        cat("\tLast observed data year",tail(rownames(fbartable(EM_update)),1),"\n")
        cat("\tAssessment year",yr[1],"\n")
        cat("\tCurrent SSB",ssb[yr[1],,"True"],"\n")

        fcThisYear <- AdviceForecastSettings
        ## Insert advice in forecast constraints if requested
        if(is(EM,"msam")){
            fcThisYear$constraints <- lapply(seq_along(fcThisYear$constraints), function(s){
                v <- gsub("%ADVICE%",catch[yr,s,"Advice"],fcThisYear$constraints[[i]])
                v[grepl("C=NA",fcThisYear$constraints)] <- NA
                v
            })
        }else{
            fcThisYear$constraints <- gsub("%ADVICE%",catch[yr,"Total","Advice"],fcThisYear$constraints)
            fcThisYear$constraints[grepl("C=NA",fcThisYear$constraints)] <- NA
        }
        if(forecastMethod == "Basic"){
            adviceForecast <- try({do.call(modelforecast, c(list(fit = EM_update, progress=FALSE), fcThisYear))})
        }else if(forecastMethod == "ICES"){
            if(is(EM,"msam")){
                cstr0 <- fcThisYear$constraints
                ## Forecast with F=Ftarget
                for(s in seq_along(EM))
                    fcThisYear$constraints[[s]] <- c(head(fcThisYear$constraints[[s]],-1),
                                                     sprintf("F=%f",EMReferencePoints[[s]]$Ftarget),
                                                     "F=1*")
                adviceForecast <- try({do.call(modelforecast, c(list(fit = EM_update, progress=FALSE), fcThisYear))})
                afFTab <- lapply(adviceForecast,attr, which = "tab")
                tabLab <- lapply(adviceForecast,attr, which = "estimateLabel")
                ## Check if SSB at beginning of advice year is < Btrigger
                redoForecast <- FALSE
                for(s in seq_along(EM)){
                    fcorr <- afFTab[[s]][yr_tac,sprintf("ssb:%s",tabLab[[s]])] / EMReferencePoints[[s]]$Btrigger
                    if(fcorr < 1){
                        fcThisYear$constraints[[s]][length(fcThisYear$constraints[[s]])-1] <- sprintf("F=%f",EMReferencePoints[[s]]$Ftarget * fcorr)
                        redoForecast <- TRUE
                    }
                }
                if(redoForecast){
                    adviceForecast <- try({do.call(modelforecast, c(list(fit = EM_update, progress=FALSE), fcThisYear))})
                    afFTab <- lapply(adviceForecast,attr, which = "tab")
                    tabLab <- lapply(adviceForecast,attr, which = "estimateLabel")
                }
                ## Check if SSB at beginning of year after advice is < Blim
                redoForecast <- FALSE
                for(s in seq_along(EM)){
                    nssb <- afFTab[[s]][as.character(as.numeric(yr_tac)+1),sprintf("ssb:%s",tabLab[[s]])]
                    if(nssb < EMReferencePoints[[s]]$Blimit){
                        fcThisYear$constraints[[s]][length(fcThisYear$constraints[[s]])-1] <- sprintf("SSB=%f",EMReferencePoints[[s]]$Blimit)
                        redoForecast <- TRUE
                    }
                }
                if(redoForecast){
                    adviceForecast <- try({do.call(modelforecast, c(list(fit = EM_update, progress=FALSE), fcThisYear))})
                    afFTab <- lapply(adviceForecast,attr, which = "tab")
                    tabLab <- lapply(adviceForecast,attr, which = "estimateLabel")
                }
                ## Check again
                redoForecast <- FALSE
                for(s in seq_along(EM)){
                    nssb <- afFTab[[s]][as.character(as.numeric(yr_tac)+1),sprintf("ssb:%s",tabLab[[s]])]
                    if(nssb < EMReferencePoints[[s]]$Blimit){
                        fcThisYear$constraints[[s]][length(fcThisYear$constraints[[s]])-1] <- sprintf("F=%f",1e-9)
                        redoForecast <- TRUE
                    }
                }
                if(redoForecast){
                    adviceForecast <- try({do.call(modelforecast, c(list(fit = EM_update, progress=FALSE), fcThisYear))})
                    afFTab <- lapply(adviceForecast,attr, which = "tab")
                    tabLab <- lapply(adviceForecast,attr, which = "estimateLabel")
                }
            }else{

            }
        }
        if(methods::is(adviceForecast,"try-error")){
            msg <- "Advice forecast error"
            break;
        }
        afFTab <- lapply(adviceForecast,attr, which = "tab")
        tabLab <- lapply(adviceForecast,attr, which = "estimateLabel")
        if(is(EM,"msam")){
            for(s in seq_along(EM)){
                cat("\tAdvice",paste(yr_tac,afFTab[[s]][yr_tac,sprintf("catch:%s",tabLab[[s]])],sep=": ",collapse="; "),"\n\n\n")        
                catch[yr_tac,,"Advice"] <- afFTab[[s]][yr_tac,sprintf("catch:%s",tabLab[[s]])]
                fbar[yr_tac,,"Advice"] <- afFTab[[s]][yr_tac,sprintf("fbar:%s",tabLab[[s]])]
                ssb[yr_tac,,"Advice"] <- afFTab[[s]][yr_tac,sprintf("ssb:%s",tabLab[[s]])]
                rec[yr_tac,,"Advice"] <- afFTab[[s]][yr_tac,sprintf("rec:%s",tabLab[[s]])]
            }
            ## Total
## MADE IT TO HERE!!
        }else{

        }
        ## Simulate next year (implement different rules for AdviceYears > 0)
        capture.output(OM_update <- try({addSimulatedYears(OM_update, 
                                                          constraints = AdviceToCatchConstraint(rep(catch[yr_tac,,"Advice"],length.out = AdviceYears)),
                                                          ...)}))
        if(methods::is(OM_update,"try-error")){
            msg <- "Adding simulated year error"
            break;
        }
        ssb[yr_tac,,"True"] <- ssbtable(OM_update)[yr_tac,1]
        fbar[yr_tac,,"True"] <- fbartable(OM_update)[yr_tac,1]
        rec[yr_tac,,"True"] <- rectable(OM_update)[yr_tac,1]
        catch[yr_tac,,"True"] <- catchtable(OM_update)[yr_tac,1]
        
        ## New assessment
        EM_update <- try({updateAssessment(OM_update, EM_update, knotRange, AdviceLag, intermediateFleets)})
        if(!methods::is(EM_update,"sam")){
            msg <- "Assessment error"
            break;
        }
        ssb[yr,,c("Estimate","Low","High")] <- ssbtable(EM_update)[yr,,drop=FALSE]
        fbar[yr,,c("Estimate","Low","High")] <- fbartable(EM_update)[yr,,drop=FALSE]
        rec[yr,,c("Estimate","Low","High")] <- rectable(EM_update)[yr,,drop=FALSE]
        catch[yr,,c("Estimate","Low","High")] <- catchtable(EM_update)[yr,,drop=FALSE]        
    }
    ## Last year assessments
    if(msg == "OK")
        for(i0 in seq_len(AdviceLag)){
            yr <- rownames(ssb)[nYears + i0]       
            cat("Assessment for last intermediate year",i0,"/",AdviceLag,"\n")
            cat("\tLast true data year",tail(rownames(fbartable(OM_update)),1),"\n")
            cat("\tLast observed data year",tail(rownames(fbartable(OM_update)),1),"\n")
            cat("\tAssessment year",yr[1],"\n")
            EM_update <- try({updateAssessment(OM_update, EM_update, knotRange, AdviceLag - i0,intermediateFleets)})
            if(methods::is(EM_update,"sam")){
                msg <- "Assessment error"
                break;
            }
            ssb[yr,,c("Estimate","Low","High")] <- ssbtable(EM_update)[yr,,drop=FALSE]
            fbar[yr,,c("Estimate","Low","High")] <- fbartable(EM_update)[yr,,drop=FALSE]
            rec[yr,,c("Estimate","Low","High")] <- rectable(EM_update)[yr,,drop=FALSE]
            catch[yr,,c("Estimate","Low","High")] <- catchtable(EM_update)[yr,,drop=FALSE]        
        }
    

    ## Return
    list(OM = OM_update,
         EM = EM_update,
         ssb = ssb,
         fbar = fbar,
         rec = rec,
         catch = catch,
         msg = msg)

}
