## TODO:
## - [ ] Where should genetic data be added?
## - [ ] Add split + separate assessments EM
## - [ ] Extend with custom simulation model (e.g. Kat/NS example)
## - [ ] Add function / map to combine fleets
## - [ ] Single fleet model needs to aggregate data
## - [ ] Add spict EM
## - [ ] Add empirical (life history, chr, rfb, rb, ome-over-two rbf) EM
## - [ ] Add options for misspecification in EM

cAdd <- function(x,y){
    as.character(as.numeric(as.character(x))+as.numeric(as.character(y)))
}


addSimulatedYears <- function(fit, constraints, ...){
    cat("Hello from generic addSimulatedYears\n")
    UseMethod("addSimulatedYears")
}

## Function to add simulated years to fit for single-stock sam
##' @method addSimulatedYears sam
##' @export
addSimulatedYears.sam <- function(fit, constraints,resampleFirst=FALSE,trueSel=NULL,refit=FALSE,silent=TRUE, resampleParameters = FALSE, deterministicF=TRUE, maxTrueF = 3.0,maxScaleF=1.5, ...){
    cat("Hello from addSimulatedYears.sam\n")
    require(stockassessment)
    if(resampleParameters){
        p0 <- rmvnorm(1, fit$opt$par, solve(fit$opt$he), TRUE)
        fit$obj$fn(p0)
        pl0 <- fit$obj$env$parList(par = fit$obj$env$last.par)        
    }else{
        pl0 <- NULL
    }
    ## Check that catch is < TSB * maxCatchOFTSB
    ## TSB <- tail(tsbtable(fit),1)[1]
    ## for(i in seq_along(constraints)){
    ##     v <- as.numeric(gsub("(C=)(.+)","\\2",constraints[[i]]))
    ##     constraints[[i]] <- sprintf("C=%f",pmin(v,TSB[[i]][1] * maxCatchOfTSB))
    ## }
    constraints <- paste0(constraints,sprintf("|F=%f & F=%f*",maxTrueF,maxScaleF))
    doSim <- modelforecast(fit, constraints, nosim=1, returnObj=2,addDataYears=TRUE,resampleFirst=resampleFirst, useModelLastN = FALSE,customSel = trueSel, custom_pl = pl0, deterministicF=deterministicF, ...)
    v <- doSim()
    obj <- environment(doSim)$obj
    ## Remove biology if reported without "dat."
    nmRm <- gsub("^dat\\.","",names(v)[grep("^dat\\.",names(v))])
    v[intersect(nmRm,names(v))] <- NULL
    ## Remove "dat." from names of biology
    
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
##' @method addSimulatedYears msam
##' @export
addSimulatedYears.msam <- function(fit, constraints,resampleFirst=FALSE,trueSel=NULL,refit=FALSE,silent=TRUE, resampleParameters = FALSE, deterministicF=TRUE, maxTrueF = 3.0, maxScaleF=1.5, ...){
    cat("Adding simulated year with constraints:\n")
    cat("Model year range:",attr(fit,"m_data")$minYearAll,attr(fit,"m_data")$maxYearAll,"(",attr(fit,"m_data")$maxYearAll-attr(fit,"m_data")$minYearAll+1,")","\n")    
    cat(paste(paste0("\t",getStockNames(fit),": ",constraints),collapse="\n"),"\n")
    if(resampleParameters){
        cat("\tResampling parameters\n")
        p0 <- rmvnorm(1, attr(fit,"m_opt")$par, solve(attr(fit,"m_opt")$he), TRUE)
        obj <- attr(fit,"m_obj")
        obj$fn(p0)
        pl0 <- obj$env$parList(par = obj$env$last.par)        
    }else{
        pl0 <- NULL # attr(fit,"m_pl")
    }
    if(!is.na(maxTrueF) || !is.na(maxScaleF))
        for(i in seq_along(constraints)){
            str <- paste(c(sprintf("F=%f",maxTrueF),sprintf("F=%f*",maxScaleF))[c(!is.na(maxTrueF),!is.na(maxScaleF))],collapse=" & ")
            constraints[[i]] <- paste0(constraints[[i]],str)
        }
    cat("\tUpdated constraints with maxTrueF\n")
    cat(paste(paste0("\t",getStockNames(fit),": ",constraints),collapse="\n"),"\n")
    doSim <- modelforecast(fit, constraints, nosim=1, returnObj=2,addDataYears=TRUE,resampleFirst=resampleFirst, useModelLastN = FALSE,customSel = trueSel, custom_pl = pl0, deterministicF=deterministicF, ...)
    cat("Ready to simulate...\n")
    v <- doSim()
    cat("Get obj...\n")
    obj <- environment(doSim)$obj
    ## Remove biology if reported without "dat."
    nmRm <- gsub("^dat\\.","",names(v)[grep("^dat\\.",names(v))])
    v[intersect(nmRm,names(v))] <- NULL
    tmpSW <- lapply(fit,function(x) x$data$stockMeanWeight)
    tmpPM <- lapply(fit,function(x) x$data$propMat)
    tmpNM <- lapply(fit,function(x) x$data$natMor)
    tmpCW <- lapply(fit,function(x) x$data$catchMeanWeight)
    ## Remove "dat." from names of biology 
    names(v) <- gsub("dat\\.","",names(v))    
    dat <- obj$env$data ## Already updated dimensions by modelforecast
    ## dat <- dat[!duplicated(names(dat))]
    ## dat$sam <- obj$env$data$sam
    cat("\tUpdate everything\n")
    for(i in seq_along(fit)){
        ##Update biology (except biopar) + logobs
        nms <- intersect(names(dat$sam[[i]]), names(v))
        for(nn in nms)
            dat$sam[[i]][[nn]] <- v[[nn]][[i]]
        dat$sam[[i]]$years <- min(as.numeric(dat$sam[[i]]$aux[,"year"])):(max(as.numeric(dat$sam[[i]]$aux[,"year"])))
        dat$sam[[i]]$noYears <- length(dat$sam[[i]]$years)
        ## Update names of biology
        dmnm <- list(dat$sam[[i]]$years, fit[[i]]$conf$minAge:fit[[i]]$conf$maxAge, dimnames(fit[[i]]$catchMeanWeight)[3])        
        ## Save observed biology
        if(!is.null(v[["obs_stockMeanWeight"]])){
            dat$sam[[i]]$stockMeanWeight <- v[["obs_stockMeanWeight"]][[i]][seq_along(dat$sam[[i]]$years),]
        }##else if(!is.null(v[["dat.stockMeanWeight"]])){
        ##     dat$sam[[i]]$stockMeanWeight <- v[["dat.stockMeanWeight"]][[i]]
        ## }
        dat$sam[[i]]$stockMeanWeight[1:nrow(tmpSW[[i]]),] <- tmpSW[[i]]
        if(!is.null(v[["obs_propMat"]])){
            dat$sam[[i]]$propMat <- v[["obs_propMat"]][[i]][seq_along(dat$sam[[i]]$years),]
        }## else if(!is.null(v[["dat.propMat"]])){
        ##     dat$sam[[i]]$propMat <- v[["dat.propMat"]][[i]]
        ## }
        dat$sam[[i]]$propMat[1:nrow(tmpPM[[i]]),] <- tmpPM[[i]]
        if(!is.null(v[["obs_natMor"]])){
            dat$sam[[i]]$natMor <- v[["obs_natMor"]][[i]][seq_along(dat$sam[[i]]$years),]
        }## else if(!is.null(v[["dat.natMor"]])){
        ##     dat$sam[[i]]$natMor <- v[["dat.natMor"]][[i]]
        ## }
        dat$sam[[i]]$natMor[1:nrow(tmpNM[[i]]),] <- tmpNM[[i]]        
        if(!is.null(v[["obs_catchMeanWeight"]])){
            dat$sam[[i]]$catchMeanWeight <- v[["obs_catchMeanWeight"]][[i]][seq_along(dat$sam[[i]]$years),,]
        }## else if(!is.null(v[["dat.catchMeanWeight"]])){
        ##     dat$sam[[i]]$catchMeanWeight <- v[["dat.catchMeanWeight"]][[i]]
        ## }
        dat$sam[[i]]$catchMeanWeight[1:nrow(tmpCW[[i]]),,] <- tmpCW[[i]]
        ## Bio using averages
        dat$sam[[i]]$propM <- v[["bio_propM"]][[i]]
        dat$sam[[i]]$propF <- v[["bio_propF"]][[i]]
        dat$sam[[i]]$landFrac <- v[["bio_landFrac"]][[i]]
        ## Update dimnames
        cat("propMat",dim(dat$sam[[i]]$propMat),length(dmnm[[1]]),length(dmnm[[2]]),"\n")
        dimnames(dat$sam[[i]]$propMat) <- dmnm[1:2] #list(dmnm[[1]][nrow(dat$sam[[i]]$propMat)],dmnm[[2]])
        cat("SW",dim(dat$sam[[i]]$stockMeanWeight),length(dmnm[[1]]),length(dmnm[[2]]),"\n")
        dimnames(dat$sam[[i]]$stockMeanWeight) <- dmnm[1:2] #list(dmnm[[1]][nrow(dat$sam[[i]]$stockMeanWeight)],dmnm[[2]])
        cat("natMor",dim(dat$sam[[i]]$natMor),length(dmnm[[1]]),length(dmnm[[2]]),"\n")
dimnames(dat$sam[[i]]$natMor) <- dmnm[1:2] #list(dmnm[[1]][nrow(dat$sam[[i]]$natMor)],dmnm[[2]])
        cat("propM",dim(dat$sam[[i]]$propM),length(dmnm[[1]]),length(dmnm[[2]]),"\n")
        dimnames(dat$sam[[i]]$propM) <- dmnm[1:2] #list(dmnm[[1]][nrow(dat$sam[[i]]$propM)],dmnm[[2]])
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
        d <- dat$sam[[i]][nn]
        attr(d,"fleetNames") <- attr(fit[[i]]$data,"fleetNames")
        d
    })
    mpl <- obj$env$parList(par=obj$env$last.par.best)
    ## Replace simulated values
    mpl$logN <- combineParameter(v$logN)
    mpl$logF <- combineParameter(v$logFs)
    mpl$logitFseason <- combineParameter(v$logitFseason)
    ## BIOPAR
    map <- lapply(fit,function(x) x$obj$env$map)
    mmap <- obj$env$map
    if(refit){
        cat("\tRefit\n")
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
        cat("\trun multisam.fit\n")
        newFit <- do.call(multisam.fit,cc, envir = ee)
    }else{
        cat("\tDo not refit\n")
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

        ## Make sure parameter lengths fit
        ref_pars <- collect_pars(newFitS)
        for (nm in intersect(names(mpl), names(ref_pars))) {
            if (length(mpl[[nm]]) != length(ref_pars[[nm]])) {
                mpl[[nm]] <- ref_pars[[nm]]
            }
            if (is.integer(mpl[[nm]])) storage.mode(mpl[[nm]]) <- "double"
        }
        
        ## ADD OTHER OPTIONS!
        ee <- attr(fit,"m_envir")
        cc <- as.list(attr(fit,"m_call"))[-1]
        ee$newFitS <- newFitS
        ee$sharedObs <- dat$sharedObs
        cc$x <- newFitS
        cc$shared_data <- as.name("sharedObs")
        ##cc$map <- mmap ## NOT IMPLEMENTED??
        ## Anything related to shared selectivity needs to be removed
        mpl$shared_phbeta <- combineParameter(replicate(length(newFitS),numeric(0),FALSE))
        cc$shared_selectivity <- 0
        cc$shared_seasonality <- 0
        cc$parlist <- mpl
        cc$run <- FALSE        
        cat("\tmake obj\n")
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
        cat("\tsdreport\n")        
        obj2 <- TMB::MakeADFun(obj$env$data, obj$env$parameters, type = "ADFun", 
                               ADreport = TRUE, DLL = obj$env$DLL, silent = obj$env$silent)
        sdv <- obj2$fn()
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
    cat("\tAll done\n")    
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


ICESAdviceForecast <- function(EM, ...){
    UseMethod("ICESAdviceForecast")
}



ICESAdviceForecast.sam <- function(EM_update,OM_update,fcThisYear,EMReferencePoints,adviceRules, yr,yr_tac,backcorrected = FALSE, ...){
    ## When to check SSB against MSYBtrigger and Blim
    dySSBAR <- NA_real_
    dySSBZC <- NA_real_
    ## If given, use that
    if(!is.null(EMReferencePoints$dySSBAR)){
        dySSBAR <- EMReferencePoints$dySSBAR
    }else if(any(EM_update$data$propM > 0) || any(EM_update$data$propF > 0)){ ## Spawning in year
        ##SSB year before advice year for advice rule
        dySSBAR <- -1
    }else{
        dySSBAR <- 0
    }
    ## "ICES advice will be based on bringing the stock above Blim at the end of the projection year"
    if(!is.null(EMReferencePoints$dySSBZC)){
        dySSBZC <- EMReferencePoints$dySSBARZC
    }else if(mean(EM_update$data$propF) > 0.5){ 
        ## If substantial fishing before spawning, use advice year to check for zero catch advice              
        dySSBZC <- 0
    }else{
        ## Otherwise, use advice year + 1
        dySSBZC <- 1
    }
                                        #cstr0 <- fcThisYear$constraints
    ## Forecast with F=Ftarget
    ## F = FMSY when SSB is at or above MSY Btrigger
    cat("\t\tInitial ICES forecast...\n")
    fcThisYear$constraints <- c(head(fcThisYear$constraints,-2),
                                sprintf("F=%f",EMReferencePoints$Ftarget),
                                "F=1*")
    ##fcThisYear$fastFixedF <- TRUE
    cat("nosim:",fcThisYear$nosim,"\n")
    if(backcorrected){
        adviceForecast <- try({do.call(stockassessment:::backcorrected_modelforecast.sam, c(list(fit = EM_update, progress=FALSE), fcThisYear))})
    }else{
        adviceForecast <- try({do.call(modelforecast, c(list(fit = EM_update, progress=FALSE), fcThisYear))})
    }
    afFTab <- attr(adviceForecast,"tab")
    tabLab <- attr(adviceForecast,"estimateLabel")
    ## Check if SSB at beginning of advice year is < Btrigger
    adviceRules[yr_tac,"Total"] <- "ICES MSY"
    redoForecast <- FALSE
    ## Check if SSB < MSY Btrigger at the "beginning of advice year"
    fcorr <- afFTab[cAdd(yr_tac,dySSBAR),sprintf("ssb:%s",tabLab)] / EMReferencePoints$Btrigger
    cat("\n\n",afFTab,"\n\n")
    cat("fcorr:", fcorr, cAdd(yr_tac,dySSBAR), sprintf("ssb:%s",tabLab), afFTab[cAdd(yr_tac,dySSBAR),sprintf("ssb:%s",tabLab)] , EMReferencePoints$Btrigger,"\n")
    if(fcorr < 1){
        ## F = FMSY × SSB/MSY Btrigger when the stock is below MSY Btrigger and above Blim (or below Blim but forecasted SSB > Blim
        fcThisYear$constraints[length(fcThisYear$constraints)-1] <- sprintf("F=%f",EMReferencePoints$Ftarget * fcorr)
        redoForecast <- TRUE
        adviceRules[yr_tac,"Total"] <- "ICES MSY PA (Below MSYBtrigger)"
    }
    if(redoForecast){
        cat("\t\tICES forecast, below MSYBtrigger...\n")
        if(backcorrected){
            adviceForecast <- try({do.call(stockassessment:::backcorrected_modelforecast.sam, c(list(fit = EM_update, progress=FALSE), fcThisYear))})
        }else{
            adviceForecast <- try({do.call(modelforecast, c(list(fit = EM_update, progress=FALSE), fcThisYear))})
        }
        afFTab <- attr(adviceForecast,"tab")
        tabLab <- attr(adviceForecast,"estimateLabel")
    }
    ## Check if SSB at "end of projection year" is < Blim
    redoForecast <- FALSE
    nssb <- afFTab[cAdd(yr_tac,dySSBZC),sprintf("ssb:%s",tabLab)]
    cat("ssb: ",nssb,"\n")
    cat("Blim: ", EMReferencePoints$Blimit,"\n")
    if(nssb < EMReferencePoints$Blimit){
        ##If so, forecast to reach Blim
        fcThisYear$constraints[length(fcThisYear$constraints)-1] <- sprintf("SSB=%f",EMReferencePoints$Blimit)
        redoForecast <- TRUE
        adviceRules[yr_tac,"Total"] <- "ICES MSY PA (Below Blim)"
    }
    if(redoForecast){
        cat("\t\tICES forecast, below Blim...\n")
        if(backcorrected){
            adviceForecast <- try({do.call(stockassessment:::backcorrected_modelforecast.sam, c(list(fit = EM_update, progress=FALSE), fcThisYear))})
        }else{
            adviceForecast <- try({do.call(modelforecast, c(list(fit = EM_update, progress=FALSE), fcThisYear))})
        }        
        afFTab <- attr(adviceForecast,"tab")
        tabLab <- attr(adviceForecast,"estimateLabel")
    }
    ## Check again
    redoForecast <- FALSE
    nssb <- afFTab[cAdd(yr_tac,dySSBZC),sprintf("ssb:%s",tabLab)]
    cat("ssb: ",nssb,"\n")
    cat("Blim: ", EMReferencePoints$Blimit,"\n")
    if(!is.finite(nssb) || nssb < EMReferencePoints$Blimit){
        ## If SSB is below Blim, zero catch advice 
        fcThisYear$constraints[length(fcThisYear$constraints)-1] <- sprintf("F=%f",1e-4)
        redoForecast <- TRUE
        adviceRules[yr_tac,"Total"] <- "ICES MSY PA (Zero catch)"
    }
    if(redoForecast){
        cat("\t\tICES forecast, Zero catch advice...\n")
        if(backcorrected){
            adviceForecast <- try({do.call(stockassessment:::backcorrected_modelforecast.sam, c(list(fit = EM_update, progress=FALSE), fcThisYear))})
        }else{
            adviceForecast <- try({do.call(modelforecast, c(list(fit = EM_update, progress=FALSE), fcThisYear))})
        }
        afFTab <- attr(adviceForecast,"tab")
        tabLab <- attr(adviceForecast,"estimateLabel")
    }    
    return(list(adviceForecast = adviceForecast,
                afFTab = afFTab,
                tabLab = tabLab,
                adviceRules = adviceRules))
}



ICESAdviceForecast.msam <- function(EM_update,OM_update,fcThisYear,EMReferencePoints,adviceRules, yr,yr_tac,backcorrected=FALSE, ...){
    ## When to check SSB against MSYBtrigger and Blim
    dySSBAR <- rep(NA_real_,length(EM_update))
    dySSBZC <- rep(NA_real_,length(EM_update))
    for(s in seq_along(EM_update)){
        ## If given, use that
        if(!is.null(EMReferencePoints[[s]]$dySSBAR)){
            dySSBAR[s] <- EMReferencePoints[[s]]$dySSBAR
        }else if(any(EM_update[[s]]$data$propM > 0) || any(EM_update[[s]]$data$propF > 0)){ ## Spawning in year
            ##SSB year before advice year for advice rule
            dySSBAR[s] <- -1
        }else{
            dySSBAR[s] <- 0
        }
        ## "ICES advice will be based on bringing the stock above Blim at the end of the projection year"
        if(!is.null(EMReferencePoints[[s]]$dySSBZC)){
            dySSBZC[s] <- EMReferencePoints[[s]]$dySSBARZC
        }else if(mean(EM_update[[s]]$data$propF) > 0.5){ 
            ## If substantial fishing before spawning, use advice year to check for zero catch advice              
            dySSBZC[s] <- 0
        }else{
            ## Otherwise, use advice year + 1
            dySSBZC[s] <- 1
        }
    }
                                        #cstr0 <- fcThisYear$constraints
    ## Forecast with F=Ftarget
    ## F = FMSY when SSB is at or above MSY Btrigger
    cat("\t\tInitial ICES forecast...\n")
    for(s in seq_along(EM_update))
        fcThisYear$constraints[[s]] <- c(head(fcThisYear$constraints[[s]],-2),
                                         sprintf("F=%f",EMReferencePoints[[s]]$Ftarget),
                                         "F=1*")
    ##fcThisYear$fastFixedF <- TRUE
    if(backcorrected){
        adviceForecast <- try({do.call(multiStockassessment:::backcorrected_modelforecast.msam, c(list(fit = EM_update, progress=FALSE), fcThisYear))})
    }else{
        adviceForecast <- try({do.call(modelforecast, c(list(fit = EM_update, progress=FALSE), fcThisYear))})
    }
    afFTab <- lapply(adviceForecast,attr, which = "tab")
    tabLab <- lapply(adviceForecast,attr, which = "estimateLabel")
    ## Check if SSB at beginning of advice year is < Btrigger
    adviceRules[cAdd(yr_tac,0),seq_along(EM_update)] <- "ICES MSY"
    redoForecast <- FALSE
    for(s in seq_along(EM_update)){
        ## Check if SSB < MSY Btrigger at the "beginning of advice year"
        nssb <- afFTab[[s]][cAdd(yr_tac,dySSBAR[s]),sprintf("ssb:%s",tabLab[[s]])]
        cat("\t\tStock: ",s,"\n")
        cat("\t\t\tforecast ssb: ",nssb,"\n")
        cat("\t\t\tBlim: ", EMReferencePoints[[s]]$Blimit,"\n")
        fcorr <- nssb / EMReferencePoints[[s]]$Btrigger                    
        if(fcorr < 1){
            ## F = FMSY × SSB/MSY Btrigger when the stock is below MSY Btrigger and above Blim (or below Blim but forecasted SSB > Blim           
            fcThisYear$constraints[[s]][length(fcThisYear$constraints[[s]])-1] <- sprintf("F=%f",EMReferencePoints[[s]]$Ftarget * fcorr)
            redoForecast <- TRUE
            adviceRules[cAdd(yr_tac,0),s] <- "ICES MSY PA (Below MSYBtrigger)"
        }
    }
    if(redoForecast){
        cat("\t\tICES forecast, below MSYBtrigger...\n")
        cat("\t\tNew constraints: ",paste(fcThisYear$constraints,collapse=", "),"\n")
        if(backcorrected){
            adviceForecast <- try({do.call(multiStockassessment:::backcorrected_modelforecast.msam, c(list(fit = EM_update, progress=FALSE), fcThisYear))})
        }else{
            adviceForecast <- try({do.call(modelforecast, c(list(fit = EM_update, progress=FALSE), fcThisYear))})
        }
        afFTab <- lapply(adviceForecast,attr, which = "tab")
        tabLab <- lapply(adviceForecast,attr, which = "estimateLabel")
    }
    ## Check if SSB at "end of projection year" is < Blim
    redoForecast <- FALSE
    for(s in seq_along(EM_update)){
        nssb <- afFTab[[s]][cAdd(yr_tac,dySSBZC[s]),sprintf("ssb:%s",tabLab[[s]])]
        cat("Stock: ",s,"\n")
        cat("\tssb: ",nssb,"\n")
        cat("\tBlim: ", EMReferencePoints[[s]]$Blimit,"\n")
        if(nssb < EMReferencePoints[[s]]$Blimit){
                                        #If so, forecast to reach Blim (SWITCH TO FIND AN F THAT GIVES 50% > Blim? - Should be the same??)
            fcThisYear$constraints[[s]][length(fcThisYear$constraints[[s]])-1] <- sprintf("SSB=%f",EMReferencePoints[[s]]$Blimit)
            redoForecast <- TRUE
            adviceRules[cAdd(yr_tac,0),s] <- "ICES MSY PA (Below Blim)"
        }
    }
    if(redoForecast){
        cat("\t\tICES forecast, below Blim...\n")
        cat("\t\tNew constraints: ",paste(fcThisYear$constraints,collapse=", "),"\n")
        if(backcorrected){
            adviceForecast <- try({do.call(multiStockassessment:::backcorrected_modelforecast.msam, c(list(fit = EM_update, progress=FALSE), fcThisYear))})
        }else{
            adviceForecast <- try({do.call(modelforecast, c(list(fit = EM_update, progress=FALSE), fcThisYear))})
        }
        afFTab <- lapply(adviceForecast,attr, which = "tab")
        tabLab <- lapply(adviceForecast,attr, which = "estimateLabel")
    }
    ## Check again
    redoForecast <- FALSE
    for(s in seq_along(EM_update)){
        nssb <- afFTab[[s]][cAdd(yr_tac,dySSBZC[s]),sprintf("ssb:%s",tabLab[[s]])]
        cat("Stock: ",s,"\n")
        cat("\tssb: ",nssb,"\n")
        cat("\tBlim: ", EMReferencePoints[[s]]$Blimit,"\n")
        if(nssb < EMReferencePoints[[s]]$Blimit){
            ## If SSB is below Blim, zero catch advice 
            fcThisYear$constraints[[s]][length(fcThisYear$constraints[[s]])-1] <- sprintf("F=%f",1e-4)
            redoForecast <- TRUE
            adviceRules[cAdd(yr_tac,0),s] <- "ICES MSY PA (Zero catch)"
        }
    }
    if(redoForecast){
        cat("\t\tICES forecast, Zero catch advice...\n")
        cat("\t\tNew constraints: ",paste(fcThisYear$constraints,collapse=", "),"\n")
        if(backcorrected){
            adviceForecast <- try({do.call(multiStockassessment:::backcorrected_modelforecast.msam, c(list(fit = EM_update, progress=FALSE), fcThisYear))})
        }else{
            adviceForecast <- try({do.call(modelforecast, c(list(fit = EM_update, progress=FALSE), fcThisYear))})
        }
        afFTab <- lapply(adviceForecast,attr, which = "tab")
        tabLab <- lapply(adviceForecast,attr, which = "estimateLabel")
    }
    ## Precautionary reduction
    ## NOTE: Handle the special case when intermediate year has F=0!
    cat(adviceRules[cAdd(yr_tac,0),seq_along(EM_update)],"\n")
    hasZeroCA <- any(adviceRules[cAdd(yr_tac,0),seq_along(EM_update)] == "ICES MSY PA (Zero catch)")    
    hasBelowTrigger <- any(adviceRules[cAdd(yr_tac,0),seq_along(EM_update)] != "ICES MSY")
    cat("hasZeroCA: ",hasZeroCA,"\n")
    cat("hasBelowTrigger: ",hasBelowTrigger,"\n")
    FRedu <- sapply(1:3,function(q){
        lastF <- afFTab[[q]][cAdd(yr_tac,-1),sprintf("fbar:%s",tabLab[[s]])]
        newF <- afFTab[[q]][cAdd(yr_tac,0),sprintf("fbar:%s",tabLab[[s]])]
        Fmsy <- EMReferencePoints[[s]]$Ftarget
        newF / ifelse(lastF==0,Fmsy,lastF)
    })
    FRedu[adviceRules[cAdd(yr_tac,0),seq_along(EM_update)] == "ICES MSY"] <- Inf
    maxRedu <- min(FRedu,na.rm=TRUE) ## NOTE: largest reduction is minimum fraction   
    stockWithRedu <- which.min(FRedu)
    cat("maxRedu: ",maxRedu,"; stockWithRedu: ", stockWithRedu, "\n")
    redoForecast <- FALSE    
    for(s in seq_along(EM_update)){
        if(!is.null(EMReferencePoints[[s]]$PA) && EMReferencePoints[[s]]$PA){
            if(hasZeroCA){ ## If one has zero catch advice, all gets zero catch advice
                cat("\t\tStock ",s," precautionary reduction - zero catch...\n")
                fcThisYear$constraints[[s]][length(fcThisYear$constraints[[s]])-1] <- sprintf("F=%f",1e-4)
                if(adviceRules[cAdd(yr_tac,0),s] != "ICES MSY PA (Zero catch)"){
                    adviceRules[cAdd(yr_tac,0),s] <- "ICES Precautionary reduction (Zero catch)"
                }
                redoForecast <- TRUE
            }else if(hasBelowTrigger){ ## If one is below the trigger all get the same reduction in F                
                if(s != stockWithRedu){
                    cat("\t\tStock ",s," precautionary reduction - below Btrigger...\n")
                    if(is.null(EMReferencePoints[[s]]$PAcompare) || EMReferencePoints[[s]]$PAcompare == "intermediateYear"){
                        fcThisYear$constraints[[s]][length(fcThisYear$constraints[[s]])-1] <- sprintf("F=%f", pmax(afFTab[[s]][cAdd(yr_tac,-1),sprintf("fbar:%s",tabLab[[s]])] * maxRedu,1e-4) )
                    }else if(EMReferencePoints[[s]]$PAcompare == "target"){
                        fcThisYear$constraints[[s]][length(fcThisYear$constraints[[s]])-1] <- sprintf("F=%f", pmax(EMReferencePoints[[s]]$Ftarget * maxRedu,1e-4) )
                    }else{
                        stop("PAcompare should be 'intermediateYear' or 'target'")
                    }
                    adviceRules[cAdd(yr_tac,0),s] <- sprintf("ICES Precautionary reduction (%.2f%%)",(1-maxRedu)*100)
                    redoForecast <- TRUE
                }
            }
        }
    }
    if(redoForecast){
        cat("\t\tICES forecast, Precautionary reduction...\n")
        cat("\t\tNew constraints: ",paste(fcThisYear$constraints,collapse=", "),"\n")
        for(xxx in 1:length(fcThisYear$constraints[[s]]))
            cat(fcThisYear$constraints[[xxx]],"\n")
        if(backcorrected){
            adviceForecast <- try({do.call(multiStockassessment:::backcorrected_modelforecast.msam, c(list(fit = EM_update, progress=FALSE), fcThisYear))})
        }else{
            adviceForecast <- try({do.call(modelforecast, c(list(fit = EM_update, progress=FALSE), fcThisYear))})
        }
        afFTab <- lapply(adviceForecast,attr, which = "tab")
        tabLab <- lapply(adviceForecast,attr, which = "estimateLabel")
    }
    
    return(list(adviceForecast = adviceForecast,
                afFTab = afFTab,
                tabLab = tabLab,
                adviceRules = adviceRules))
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
        ## ConfNew
        if(is(EM,"msam")){
            confNew <- lapply(EM,function(x) x$conf)
            map <- lapply(EM,function(x) x$obj$env$map)
            mmap <- attr(EM,"m_obj")$env$map
        }else if(is(EM,"sam")){
            require(stockassessment)
            confNew <- EM$conf
            map <- EM$obj$env$map
            mmap <- NULL
        }
        
        if(is(EM,"sam")){
            ## Map OM fleets to EM fleets
            OM2EM <- attr(EM,"OM2EM")
            if(is.null(OM2EM))
                stop("The EM needs an OM2EM attribute to map the multi-stock data to single-stock")
            ofn <- attr(OM[[1]]$data,"fleetNames")
            ofn2 <- attr(OM,"m_data")$sharedObs
            efn <- OM2EM$fleetLevels
            cat(ofn,"\n")
            cat(efn,"\n")
            ofac <- factor(ofn,ofn,efn)
            if(datNew$sharedObs$hasSharedObs){
                ## Collapse from shared obs
                dObs <- cbind(data.frame(obs = exp(datNew$sharedObs$logobs),
                                         fleetNames = ofac[datNew$sharedObs$aux[,2]]),
                              as.data.frame(datNew$sharedObs$aux)
                              )
            }else{
                ## Collapse from stock data
                naux <- do.call(rbind,lapply(datNew$sam,function(x) x$aux))
                dObs <- cbind(data.frame(obs = exp(do.call(c,lapply(datNew$sam,function(x) x$logobs))),
                                         fleetNames = ofac[naux[,2]]),
                              as.data.frame(naux))                
            }
            ff <- ((unique(dObs$fleetNames[!is.na(as.character(dObs$fleetNames))])))
            ff <- ff[order(match(ff,attr(EM$data,"fleetNames")))]
            names(ff) <- ff
            SingleObs <- lapply(ff,function(f) stripAttr(xtabs(obs~year + age,data = dObs, fleetNames==f & age >= 0)))            
            ## Get fleet types
            ft <- EM$data$fleetTypes ##sapply(split(datNew$sam[[1]]$fleetTypes,efn),function(x) x[1])
            ## Get times
            fst <- EM$data$sampleTimesStart  ##sapply(split(datNew$sam[[1]]$sampleTimesStart,efn),function(x) mean(x))
            fet <- EM$data$sampleTimesEnd  ##sapply(split(datNew$sam[[1]]$sampleTimesEnd,efn),function(x) mean(x))
            for(i in 1:length(SingleObs))
                attr(SingleObs[[i]],"times") <- c(fst[i],fet[i])
            ## Collapse biology
            
            ## Setup new fit
            addNames <- function(x){                
                ## y <- OM[[1]]$data$years
                ## a <- seq(OM[[1]]$conf$minAge, OM[[1]]$conf$maxAge, 1)
                ## dimnames(x)[1:2] <- list(head(y,nrow(x)),a)
                x
            }
            datSS <- setup.sam.data(surveys = SingleObs[!(ft %in% c(0,7))],
                                    residual.fleets=SingleObs[(ft %in% c(0))],
                                    prop.mature=addNames(Reduce("+",lapply(datNew$sam,function(x)x$propMat))/3), 
                                    stock.mean.weight=addNames(Reduce("+",lapply(datNew$sam,function(x)x$stockMeanWeight))/3), 
                                    catch.mean.weight=addNames(Reduce("+",lapply(datNew$sam,function(x)x$catchMeanWeight))/3), 
                                    dis.mean.weight=addNames(Reduce("+",lapply(datNew$sam,function(x)x$disMeanWeight))/3), 
                                    land.mean.weight=addNames(Reduce("+",lapply(datNew$sam,function(x)x$landMeanWeight))/3),
                                    prop.f=addNames(Reduce("+",lapply(datNew$sam,function(x)x$propF))/3), 
                                    prop.m=addNames(Reduce("+",lapply(datNew$sam,function(x)x$propM))/3), 
                                    natural.mortality=addNames(Reduce("+",lapply(datNew$sam,function(x)x$natMor))/3), 
                                    land.frac=addNames(Reduce("+",lapply(datNew$sam,function(x)x$landFrac))/3))


            ## Handle advice year lag        
            if(AdviceLag > 0){
                if(length(intermediateFleets) == 0){ ## Reduce fully
                    datSS <- reduce(datSS, year = tail(datSS$years,AdviceLag), conf = confNew)
                    confNew <- attr(datSS,"conf")
                }else{ ## Keep some fleets in first year
                    ## First remove intermediate years 2+
                    if(AdviceLag > 1){
                        datSS <- reduce(datSS, year = tail(datSS$years,AdviceLag-1), conf = confNew)
                        confNew <- attr(datSS,"conf")
                    }
                    ## Remove fleets not in intermediateFleets for first intermediate year
                    if(is.character(intermediateFleets)){
                        intermediateFleets <- match(intermediateFleets, attr(datSS,"fleetNames"))
                        if(any(is.na(intermediateFleets)))
                            stop("intermediateFleets names does not match model fleet names")
                    }
                    datSS <- reduce(datSS,
                                    year = tail(datSS$years,AdviceLag),
                                    fleet = setdiff(seq_along(datSS$fleetTypes), intermediateFleets),
                                    conf = confNew)
                    confNew <- attr(datSS,"conf")
                }
            }                
            dp <- defpar(datSS,confNew)
            cc <- as.list(attr(EM,"call"))[-1]
            ee <- new.env()
            parent.env(ee) <- globalenv()
            list2env(as.list(attr(EM,"envir"),ee))
            list2env(as.list(getNamespace("stockassessment")),ee)
            ee$datSS <- datSS
            ee$confNew <- confNew
            ee$dp <- dp
            cc$data <- as.name("datSS")            
            cc$conf <- as.name("confNew")
            cc$parameters <- as.name("dp")
            cc$silent <- TRUE
            EM_New <- do.call(sam.fit,cc, envir = ee)
            attr(EM_New$data,"fleetNames") <- attr(EM$data,"fleetNames")
            attr(EM_New,"envir") <- attr(EM,"envir")
            attr(EM_New,"OM2EM") <- attr(EM,"OM2EM")
            return(EM_New)
        }else if(is(EM,"msam")){
            ## Map OM fleets to EM fleets
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
                        datNew$sharedObs <- reduce_shared(datNew$sharedObs,
                                                          year = tail(datNew$sam[[i]]$years, AdviceLag))
                    }
                    ## Reduce genetic data
                }else{ ## Keep some fleets in first year
                    ## First remove intermediate years 2+
                    if(AdviceLag > 1){
                        for(i in seq_along(datNew$sam)){
                            datNew$sam[[i]] <- reduce(datNew$sam[[i]],
                                                      year = tail(datNew$sam[[i]]$years, AdviceLag-1),
                                                      conf = confNew[[i]])
                            confNew[[i]] <- attr(datNew$sam[[i]],"conf")
                        }
                        ## Reduce shared data
                        if(datNew$sharedObs$hasSharedObs){
                            datNew$sharedObs <- reduce_shared(datNew$sharedObs,
                                                              year = tail(datNew$sam[[1]]$years, AdviceLag-1))
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
                                                  fleet = setdiff(seq_along(datNew$sam[[i]]$fleetTypes), intermediateFleets),
                                                  conf = confNew[[i]])
                        confNew[[i]] <- attr(datNew$sam[[i]],"conf")
                    }
                    ## Reduce shared data
                    if(datNew$sharedObs$hasSharedObs){
                        datNew$sharedObs <- reduce_shared(datNew$sharedObs,
                                                          year = tail(datNew$sam[[1]]$years, 1),
                                                          fleet = setdiff(seq_along(datNew$sam[[1]]$fleetTypes), intermediateFleets))
                    }
                }
            }

            ##confNew <- lapply(EM,function(x) x$conf) # Updated above
            ## Prepare data
            ## Map OM fleets to EM fleets
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
            ee <- attr(EM,"m_envir")
            cc <- as.list(attr(EM,"m_call"))[-1]
            ee$newFitS <- newFitS
            ee$sharedObs <- datNew$sharedObs
            ee$mpl <- attr(EM,"m_pl")
            cc$x <- as.name("newFitS")            
            if(!is.null(cc$shared_data)) ## Check that EM has shared obs
                cc$shared_data <- as.name("sharedObs")
            cc$parlist <- as.name("mpl")
            cc$silent <- TRUE
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


splitCatch <- function(C,fit, Type = c("KeepF","KeepC"), nTail = 1){
    Type <- match.arg(Type)
    if(Type == "KeepF"){rep <- attr(fit,"m_rep")
        pl <- attr(fit,"m_pl")
        NN <- sapply(fit,function(s) s$conf$maxAge-s$conf$minAge+1)
        logN <- lapply(ntable(fit,returnList=TRUE),log)
        logF <- rep$logFs
        if(!is.null(rep$catchMeanWeight)){
            CW <- rep$catchMeanWeight
        }else{
            CW <- lapply(fit, function(x) x$data$catchMeanWeight)
        }
        logHazard_M_breakpoints <- lapply(rep$mort,function(x) x$logHazard_M_breakpoints)
        logHazard_F_breakpoints <- lapply(rep$mort,function(x) x$logHazard_F_breakpoints)
        activeHazard_F <- lapply(rep$mort, function(x) x$activeHazard_F)
        activeHazard_breakpoints <- lapply(rep$mort, function(x) x$activeHazard_breakpoints)
        sampleTimesStart <- lapply(fit, function(x) x$data$sampleTimesStart)
        sampleTimesEnd <- lapply(fit, function(x) x$data$sampleTimesEnd)
        activeHazardMap_risk <- lapply(rep$mort, function(x) x$activeHazardMap_risk)
        logspace_add <- function(logx,logy){
            r <- pmax(logx,logy)
            ii <- is.finite(r)
            r[ii] <- pmax(logx[ii],logy[ii]) + log1p(exp(-abs(logx[ii]-logy[ii])))
            r
        }
        logspace_sum <- function(logx){
            if(length(logx)<=1) return(logx)
            if(length(logx)==2) return(logspace_add(logx[1],logx[2]))
            Mx <- max(logx)
            Mx + log(sum(exp(logx-Mx)))
        }
        logspace_1m <- function(logx) log1p(-exp(logx))
        logspace_1p <- function(logx) log1p(exp(logx))
        getFbar <- function(eta,s,y){
            fr <- fit[[s]]$conf$fbarRange - fit[[s]]$conf$minAge + 1
            rowMeans(exp(logF[[s]][[y]][,fr[1]:fr[2]] + eta))    
        }
        predCatch <- function(eta, s, y){
            lN <- logN[[s]][y,]
            lHM <- logHazard_M_breakpoints[[s]][,y,,,drop=FALSE] ## Age, Fleet, Season/brkpnt
            dim(lHM) <- dim(lHM)[-2]
            lHF <- logHazard_F_breakpoints[[s]][,y,,,drop=FALSE] + eta
            dim(lHF) <- dim(lHF)[-2]
                        
            logHazard_breakpoints <- logspace_add(apply(lHM,c(1,3),logspace_sum), apply(lHF,c(1,3),logspace_sum))
            logSurvBefore <- matrix(0,nrow(logHazard_breakpoints),dim(lHF)[2])
            for(f in 1:ncol(logSurvBefore)){
                t0 <- 0
                t1 <- sampleTimesStart[[s]][f]
                logS0 <- logSurvBefore[,f]
                for(t in head(seq_along(activeHazard_breakpoints[[s]]),-1)){
                    ## Skip interval if it ends before time interval
                    if(activeHazard_breakpoints[[s]][t+1] <= t0)
                        next;
                    ## Break loop if interval starts after time interval
                    if(activeHazard_breakpoints[[s]][t] >= t1)
                        break;
                    ## Otherwise, look at overlap between intervals	   
                    Astart = pmax(activeHazard_breakpoints[[s]][t],t0)
                    Aend = pmin(activeHazard_breakpoints[[s]][t+1],t1)
                    ## Hazard is constant, so cumulative hazard is hazard times interval length
                    logS0 = logS0 - exp(logHazard_breakpoints[,t]) * (Aend - Astart)
                }
                logSurvBefore[,f] <- logS0
            }
            logCIF <- matrix(-Inf,nrow(logHazard_breakpoints),dim(lHF)[2])
            for(f in 1:ncol(logSurvBefore)){
                t0 <- sampleTimesStart[[s]][f]
                t1 <- sampleTimesEnd[[s]][f]
                logS0 <- numeric(nrow(logCIF))
                vlogCIF <- rep(-Inf,nrow(logCIF))
                for(t in head(seq_along(activeHazard_breakpoints[[s]]),-1)){
                    ## Skip interval if it ends before time interval
                    if(activeHazard_breakpoints[[s]][t+1] <= t0)
                        next;
                    ## Break loop if interval starts after time interval
                    if(activeHazard_breakpoints[[s]][t] >= t1)
                        break;
                    ## Otherwise, look at overlap between intervals	   
                    ## Full interval
                    Astart = pmax(activeHazard_breakpoints[[s]][t],t0)
                    Aend = pmin(activeHazard_breakpoints[[s]][t+1],t1)
                    lCIF_F_brk <- lHF[,f,t] - logHazard_breakpoints[,t] + logspace_1m(-exp(logHazard_breakpoints[,t])*(Aend-Astart))
                    tmp <- logS0 + lCIF_F_brk
                    vlogCIF = logspace_add(vlogCIF, tmp)
                    logS0 <- logS0 - exp(logHazard_breakpoints[,t])*(Aend-Astart)
                }
                logCIF[,f] <- vlogCIF
            }
            ## Output catch, sum of all fleets
            ## logN(a,y) + mort.logFleetSurvival_before(a,y,f-1) + mort.fleetLogCumulativeIncidence(a,y,f-1);
            exp(logspace_sum(lN + logSurvBefore[,fit[[s]]$data$fleetTypes==0] + logCIF[,fit[[s]]$data$fleetTypes==0] + log(t(CW[[s]][y,,]))))      
        }
        a <- nlminb(0, function(a){
            Ca <- sapply(seq_along(fit),function(s) predCatch(a,s,nrow(logN[[1]])))
            (sum(C) - sum(Ca))^2
        })
        vv <- sapply(seq_along(fit),function(s) predCatch(a$par,s,nrow(logN[[1]])))
        return(vv)
    }else{
        v <- colMeans(tail(catchtable(fit),nTail))[(seq_len(length(fit)*3)-1)%%3 == 0]
        return(unname(sum(C) * v / sum(v)))
    }
    return(unname(C * rep(1/length(fit),length.out = length(fit))))
}


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
                AdviceLag = 1,
                initialAdvice = NA,
                implementationError = function(x) x + 1e-4,
                knotRange = 3,
                intermediateFleets = numeric(0),
                OMselectivityFixed = FALSE,
                managementToPopulation = function(x,OM){                    
                    splitCatch(x,OM,"KeepF",1)
                },
                adviceToManagement = function(advice,adviceLast) advice,
                adviceMethod = c("Basic","ICES","Empirical"),
                EMReferencePoints = NULL,
                maxTrueF = 3.0,
                maxScaleF = 1.5,
                inherentImplementationError = FALSE,
                backcorrected = FALSE,
                ...){

    cat("Setting up the MSE\n")
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
            if(any(is.na(match(c("Ftarget","Btrigger","Blimit"),names(EMReferencePoints)))))
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
        if(is(OM,"sam")){
            OM_pl <- OM$pl
            trueSel <- as.numeric(exp(t(OM_pl$logF[OM$conf$keyLogFsta[1,]+1,ncol(OM_pl$logF)])) / tail(fbartable(OM)[,1],1))
        }else{
        }
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
    adviceRules <- matrix(NA_character_,nYears+AdviceLag,nStocksOM+1)
    
    dimnames(ssb) <- dimnames(fbar) <- dimnames(rec) <- dimnames(catch) <- list(seq(do.call(max,lapply(OM,function(x)x$data$years)) + 1,len = nYears + AdviceLag),
                                                                                c(getStockNames(OM),"Total"),
                                                                                c("Advice","True","Estimate","Low","High"))
    dimnames(adviceRules) <- list(seq(do.call(max,lapply(OM,function(x)x$data$years)) + 1,len = nYears + AdviceLag),
                                  c(getStockNames(OM),"Total"))

    EMoutput <- vector("list",nYears+AdviceLag)
    AFoutput <- vector("list",nYears+AdviceLag)
    names(EMoutput) <- names(AFoutput) <- seq(do.call(max,lapply(OM,function(x)x$data$years)) + 1,len = nYears + AdviceLag)

##### Copy OM and EM for safe overwriting #####
    OM_update <- OM
    EM_update <- EM
    msg <- "OK"

##### Insert initialAdvice in results table #####
    if(is(EM,"msam")){
        if(!is.matrix(initialAdvice))
            initialAdvice <- matrix(initialAdvice,AdviceLag,length(OM),byrow=TRUE)
        for(s in seq_along(EM)){
            catch[seq_len(AdviceLag),s,"Advice"] <- initialAdvice[,s]
        }
        ## Total
        catch[seq_len(AdviceLag),"Total","Advice"] <- rowSums(initialAdvice)
        adviceRules[seq_len(AdviceLag),seq_along(OM)] <- "Initial advice"
    }else{
        ## v <- singleStockManagementToSubstock(rep(Reduce("+",initialAdvice),length.out = AdviceLag),OM)
        ## for(s in seq_along(EM)){
        ##     catch[seq_len(AdviceLag),s,"Advice"] <- v[s,]
        ## }
        catch[seq_len(AdviceLag),"Total","Advice"] <- rep(initialAdvice,length.out = AdviceLag)
        adviceRules[seq_len(AdviceLag),"Total"] <- "Initial advice"
    }

##### Helper function to convert advice (number) to forecast constraint for addSimulatedYears #####
    splitYears <- function(x) split(x,row(x))
    AdviceToCatchConstraint <- function(x,xlast, OM){
        cat("... in AdviceToCatchConstraint\n")
        cat("x: ",x,"\n")
        cat("xlast: ",xlast,"\n")
        ## By year
        xll <- splitYears(xlast)
        xl <- splitYears(x)
        v <- vector("list",length(xl))
        cat("a2m: ",adviceToManagement(xl[[1]],xll[[1]]),"\n")
        cat("m2p: ", managementToPopulation(adviceToManagement(xl[[1]],xll[[1]]),OM),"\n")        
        v[[1]] <- sprintf("C=%f",implementationError(managementToPopulation(adviceToManagement(xl[[1]],xll[[1]]),OM)))
        cat("v[[1]]: ",v[[1]],"\n")
        if(length(xl) > 1)
            for(i in 2:length(xl))
                v[[i]] <- sprintf("C=%f",implementationError(managementToPopulation(adviceToManagement(xl[[i]],xl[[i-1]]),OM)))
        ## By stock
        v2 <- do.call(rbind,v)
        r <- split(v2,col(v2))
        names(r) <- getStockNames(OM)
        r
    }
    drop3D <- function(x){
        if(is.matrix(x))
            return(x)
        dim(x) <- dim(x)[1:2]
        x
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
           length(intermediateFleets) > 0){            
            yx <- 1
            ## Else, if year.base is secondLastYear, we need to forecast one more year
        }else if(AdviceForecastSettings$year.base == "secondLastYear"){
            yx <- 1
        }
        ## Otherwise, we are good.
    }
    ## Check AdviceForecast is long ennough
    if((is(EM,"msam") && any(sapply(AdviceForecastSettings$constraints,length) < AdviceLag + AdviceYears + yx)) || (is(EM,"sam") && (length(AdviceForecastSettings$constraints) < AdviceLag + AdviceYears + yx) )){
        warning(sprintf("Length of AdviceForecastSettings$constraints should equal AdviceLag + AdviceYears + %d = %d. Modifying to match.",yx,AdviceLag+AdviceYears+yx))
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
    cat("Ready to start\n")
    if(AdviceLag > 0){
        cat("Update OM to handle advice lag\n")
        ## Update OM
        iy <- head(rownames(ssb),AdviceLag)
        if(methods::is(EM_update,"msam")){
            prevC <- matrix(do.call(c,lapply(catchtable(EM_update,FALSE,addTotal=TRUE,returnList=TRUE),function(x){
                x[as.character(min(as.numeric(iy))-1),1]
            })),nrow=1)
        }else{
            prevC <- matrix(do.call(cbind,c(replicate(length(OM_update),NA,FALSE),tail(catchtable(EM_update,FALSE)[as.character(min(as.numeric(iy))-1),1],1))),nrow=1)
        }
        cat(prevC,"\n")
        a2c <- AdviceToCatchConstraint(drop3D(catch[iy,,"Advice",drop=FALSE]),prevC,OM_update)
        cat(unlist(a2c),"\n")
        cat("...run add simulated years \n")
        cat(class(OM_update),"\n")
        OM_update <- try({addSimulatedYears(OM_update, 
                                            constraints = a2c,
                                            deterministicF = !inherentImplementationError,
                                            maxTrueF = maxTrueF,
                                            maxScaleF = maxScaleF,
                                            ...)})
        ssb[iy,,"True"] <- ssbtable(OM_update,addTotal=TRUE)[iy,seq(1,by=3,length.out=nStocksOM+1)]
        fbar[iy,,"True"] <- fbartable(OM_update,addTotal=TRUE)[iy,seq(1,by=3,length.out=nStocksOM+1)]
        rec[iy,,"True"] <- rectable(OM_update,addTotal=TRUE)[iy,seq(1,by=3,length.out=nStocksOM+1)]
        catch[iy,,"True"] <- catchtable(OM_update,addTotal=TRUE)[iy,seq(1,by=3,length.out=nStocksOM+1)]
    }
    cat("Starting assessment loop\n")      
    ## Run assessment loop
    for(i in seq(1,nYears-(AdviceYears-1), by = AdviceYears)){ # Index over assessment year
        yr <- rownames(ssb)[seq(i,len=AdviceYears)]
        yr_tac <- rownames(ssb)[seq(i+AdviceLag,len=AdviceYears)]
        
        cat("\n\n\nSimulation year",i,"\n")
        cat("\tAdvice year",yr_tac,"\n")
        cat("\tLast true data year",tail(rownames(fbartable(OM_update)),1),"\n")
        cat("\tLast observed data year",tail(rownames(fbartable(EM_update)),1),"\n")
        cat("\tAssessment year",yr[1],"\n")
        cat("\tCurrent SSB",ssb[tail(rownames(fbartable(OM_update)),1),,"True"],"\n")

        cat("\t\tUpdating assessment...\n")
        ## Update assessment
        capture.output(EM_update <- try({updateAssessment(OM_update, EM_update, knotRange, AdviceLag, intermediateFleets)}))
        
        if(!methods::is(EM_update,"sam") && !methods::is(EM_update,"msam")){
            msg <- "Assessment error"
            break;
        }
        cat("\t\tSaving output...\n")
        
        if(is(EM,"msam")){
            EMoutput[[yr]] <- list(fbar = fbartable(EM_update,returnList=TRUE,addTotal=TRUE),
                                   ssb = ssbtable(EM_update,returnList=TRUE,addTotal=TRUE),
                                   catch = catchtable(EM_update,returnList=TRUE,addTotal=TRUE),
                                   tsb = tsbtable(EM_update,returnList=TRUE,addTotal=TRUE),
                                   rec = rectable(EM_update,returnList=TRUE,addTotal=TRUE))
            tab <- ssbtable(EM_update,addTotal=TRUE,returnList=TRUE)
            for(s in seq_along(tab))
                ssb[yr,s,c("Estimate","Low","High")] <- tab[[s]][yr,,drop=FALSE]
            tab <- fbartable(EM_update,addTotal=TRUE,returnList=TRUE)
            for(s in seq_along(tab))
                fbar[yr,s,c("Estimate","Low","High")] <- tab[[s]][yr,,drop=FALSE]
            tab <- rectable(EM_update,addTotal=TRUE,returnList=TRUE)
            for(s in seq_along(tab))
                rec[yr,s,c("Estimate","Low","High")] <- tab[[s]][yr,,drop=FALSE]
            tab <- catchtable(EM_update,addTotal=TRUE,returnList=TRUE)
            for(s in seq_along(tab))
                catch[yr,s,c("Estimate","Low","High")] <- tab[[s]][yr,,drop=FALSE] 
        }else if(is(EM,"sam")){
            EMoutput[[yr]] <- list(fbar = list(Total=fbartable(EM_update)),
                                   ssb = list(Total=ssbtable(EM_update)),
                                   catch = list(Total=catchtable(EM_update)),
                                   tsb = list(Total=tsbtable(EM_update)),
                                   rec = list(Total=rectable(EM_update)))
            cat(rownames(ssbtable(OM_update)),"\n")
            cat(rownames(ssbtable(EM_update)),"\n")
            ssb[yr,"Total",c("Estimate","Low","High")] <- ssbtable(EM_update)[yr,,drop=FALSE]
            fbar[yr,"Total",c("Estimate","Low","High")] <- fbartable(EM_update)[yr,,drop=FALSE]
            rec[yr,"Total",c("Estimate","Low","High")] <- rectable(EM_update)[yr,,drop=FALSE]
            catch[yr,"Total",c("Estimate","Low","High")] <- catchtable(EM_update)[yr,,drop=FALSE] 
        }
     
        cat("\t\tPreparing forecast...\n")

        fcThisYear <- AdviceForecastSettings
        ## Insert advice in forecast constraints if requested
        if(is(EM,"msam")){
            fcThisYear$constraints <- lapply(seq_along(fcThisYear$constraints), function(s){
                v <- gsub("%ADVICE%",catch[yr,s,"Advice"],fcThisYear$constraints[[s]])
                if(any(grepl("C=NA",fcThisYear$constraints)))
                    v[grepl("C=NA",fcThisYear$constraints)] <- NA
                v
            })
        }else{
            fcThisYear$constraints <- gsub("%ADVICE%",catch[yr,"Total","Advice"],fcThisYear$constraints)
            fcThisYear$constraints[grepl("C=NA",fcThisYear$constraints)] <- NA
            cat("FCC: ",fcThisYear$constraints,"\n")
        }
        if(adviceMethod == "Basic"){
            if(backcorrected && is(EM,"sam")){
                adviceForecast <- try({do.call(stockassessment:::backcorrected_modelforecast.sam, c(list(fit = EM_update, progress=FALSE), fcThisYear))})
            }else if(backcorrected && is(EM,"msam")){
                adviceForecast <- try({do.call(multiStockassessment:::backcorrected_modelforecast.msam, c(list(fit = EM_update, progress=FALSE), fcThisYear))})
            }else{
                adviceForecast <- try({do.call(modelforecast, c(list(fit = EM_update, progress=FALSE), fcThisYear))})
            }
            if(methods::is(adviceForecast,"try-error")){
                msg <- "Advice forecast error"
                break;
            }
            if(is(EM,"msam")){
                adviceRules[yr_tac,seq_along(OM)] <- "Basic"
            }else{
                adviceRules[yr_tac,"Total"] <- "Basic"
            }
        }else if(adviceMethod == "ICES"){
            cat("\t\tICES type advice...\n")
            tmp <- try({ICESAdviceForecast(EM_update,OM_update,fcThisYear,EMReferencePoints,adviceRules,yr,yr_tac,backcorrected,...)})
            if(is(tmp,"try-error")){
                msg <- "Advice forecast error"
                break;
            }
            adviceForecast <- tmp$adviceForecast
            adviceRules <- tmp$adviceRules
            if(methods::is(adviceForecast,"try-error")){
                msg <- "Advice forecast error"
                break;
            }
            afFTab <- tmp$afFTab
            tabLab <- tmp$tabLab
            AFoutput[[yr]] <- afFTab
        }
        cat("\t\tSaving output...\n")
        ## Save advice
        if(is(EM,"msam")){
            for(s in seq_along(EM)){
                cat("\tAdvice",paste(yr_tac,afFTab[[s]][yr_tac,sprintf("catch:%s",tabLab[[s]])],sep=": ",collapse="; "),"\n\n\n")        
                catch[yr_tac,s,"Advice"] <- afFTab[[s]][yr_tac,sprintf("catch:%s",tabLab[[s]])]
                fbar[yr_tac,s,"Advice"] <- afFTab[[s]][yr_tac,sprintf("fbar:%s",tabLab[[s]])]
                ssb[yr_tac,s,"Advice"] <- afFTab[[s]][yr_tac,sprintf("ssb:%s",tabLab[[s]])]
                rec[yr_tac,s,"Advice"] <- afFTab[[s]][yr_tac,sprintf("rec:%s",tabLab[[s]])]
            }
            ## Total (Does this make sense?)
            catch[yr_tac,"Total","Advice"] <- sum(adviceToManagement(catch[yr_tac,,"Advice"],catch[cAdd(yr_tac,-1),,"Advice"]))
            #fbar[yr_tac,"Total","Advice"] <- afFTab[[s]][yr_tac,sprintf("fbar:%s",tabLab[[s]])]
            #ssb[yr_tac,"Total","Advice"] <- afFTab[[s]][yr_tac,sprintf("ssb:%s",tabLab[[s]])]
            #rec[yr_tac,"Total","Advice"] <- afFTab[[s]][yr_tac,sprintf("rec:%s",tabLab[[s]])]
        }else{
            ## cat(yr_tac,sprintf("catch:%s",tabLab),"\n")
            ## return(afFTab)
            cat("\tAdvice",paste(yr_tac,afFTab[yr_tac,sprintf("catch:%s",tabLab)],sep=": ",collapse="; "),"\n\n\n")        
            catch[yr_tac,"Total","Advice"] <- afFTab[yr_tac,sprintf("catch:%s",tabLab)]
            fbar[yr_tac,"Total","Advice"] <- afFTab[yr_tac,sprintf("fbar:%s",tabLab)]
            ssb[yr_tac,"Total","Advice"] <- afFTab[yr_tac,sprintf("ssb:%s",tabLab)]
            rec[yr_tac,"Total","Advice"] <- afFTab[yr_tac,sprintf("rec:%s",tabLab)]
        }
        ## Simulate next year (implement different rules for AdviceYears > 0)
        cat("\t\tUpdate OM...\n")
        cat("Catch constraints:\n",do.call(c,AdviceToCatchConstraint(drop3D(catch[yr_tac,,"Advice",drop=FALSE]),drop3D(catch[cAdd(yr_tac,-1),,"Advice",drop=FALSE]),OM_update)),"\n\n")
        capture.output(OM_update <- try({addSimulatedYears(OM_update, 
                                                           constraints = AdviceToCatchConstraint(drop3D(catch[yr_tac,,"Advice",drop=FALSE]),drop3D(catch[cAdd(yr_tac,-1),,"Advice",drop=FALSE]),OM_update),
                                                           deterministicF = !inherentImplementationError,
                                                           maxTrueF = maxTrueF,
                                                           maxScaleF = maxScaleF,
                                                          ...)}))
        if(methods::is(OM_update,"try-error")){
            msg <- "Adding simulated year error"
            break;
        }
        cat("\t\tSaving output...\n")
        tab <- ssbtable(OM_update,addTotal=TRUE,returnList=TRUE)
        for(s in seq_along(tab))
            ssb[cAdd(yr,1),s,"True"] <- tab[[s]][cAdd(yr,1),1]
        tab <- fbartable(OM_update,addTotal=TRUE,returnList=TRUE)
        for(s in seq_along(tab))
            fbar[cAdd(yr,1),s,"True"] <- tab[[s]][cAdd(yr,1),1]
        tab <- rectable(OM_update,addTotal=TRUE,returnList=TRUE)
        for(s in seq_along(tab))
            rec[cAdd(yr,1),s,"True"] <- tab[[s]][cAdd(yr,1),1]
        tab <- catchtable(OM_update,addTotal=TRUE,returnList=TRUE)
        for(s in seq_along(tab))
            catch[cAdd(yr,1),s,"True"] <- tab[[s]][cAdd(yr,1),1]     
    }

    ## Last year assessments
    if(msg == "OK")
        for(i0 in seq_len(AdviceLag)){
            yr <- rownames(ssb)[nYears + i0]       
            cat("Assessment for last intermediate year",i0,"/",AdviceLag,"\n")
            cat("\tLast true data year",tail(rownames(fbartable(OM_update)),1),"\n")
            cat("\tLast observed data year",tail(rownames(fbartable(OM_update)),1),"\n")
            cat("\tAssessment year",yr[1],"\n")
            ##capture.output(EM_update <- try({updateAssessment(OM_update, EM_update, knotRange, AdviceLag - i0,intermediateFleets)}))
            capture.output(EM_update <- try({updateAssessment(OM_update, EM_update, knotRange, AdviceLag, intermediateFleets)}))
            if(methods::is(EM_update,"sam")){
                msg <- "Assessment error"
                break;
            }
            if(is(EM,"msam")){
                EMoutput[[yr]] <- list(fbar = fbartable(EM_update,returnList=TRUE,addTotal=TRUE),
                                       ssb = ssbtable(EM_update,returnList=TRUE,addTotal=TRUE),
                                       catch = catchtable(EM_update,returnList=TRUE,addTotal=TRUE),
                                       tsb = tsbtable(EM_update,returnList=TRUE,addTotal=TRUE),
                                       rec = rectable(EM_update,returnList=TRUE,addTotal=TRUE))
                tab <- ssbtable(EM_update,addTotal=TRUE,returnList=TRUE)
                for(s in seq_along(tab))
                    ssb[yr,s,c("Estimate","Low","High")] <- tab[[s]][yr,,drop=FALSE]
                tab <- fbartable(EM_update,addTotal=TRUE,returnList=TRUE)
                for(s in seq_along(tab))
                    fbar[yr,s,c("Estimate","Low","High")] <- tab[[s]][yr,,drop=FALSE]
                tab <- rectable(EM_update,addTotal=TRUE,returnList=TRUE)
                for(s in seq_along(tab))
                    rec[yr,s,c("Estimate","Low","High")] <- tab[[s]][yr,,drop=FALSE]
                tab <- catchtable(EM_update,addTotal=TRUE,returnList=TRUE)
                for(s in seq_along(tab))
                    catch[yr,s,c("Estimate","Low","High")] <- tab[[s]][yr,,drop=FALSE] 
            }else if(is(EM,"sam")){
                EMoutput[[yr]] <- list(fbar = list(Total=fbartable(EM_update)),
                                       ssb = list(Total=ssbtable(EM_update)),
                                       catch = list(Total=catchtable(EM_update)),
                                       tsb = list(Total=tsbtable(EM_update)),
                                       rec = list(Total=rectable(EM_update)))
                ssb[yr,"Total",c("Estimate","Low","High")] <- ssbtable(EM_update)[yr,,drop=FALSE]
                fbar[yr,"Total",c("Estimate","Low","High")] <- fbartable(EM_update)[yr,,drop=FALSE]
                rec[yr,"Total",c("Estimate","Low","High")] <- rectable(EM_update)[yr,,drop=FALSE]
                catch[yr,"Total",c("Estimate","Low","High")] <- catchtable(EM_update)[yr,,drop=FALSE] 
            }
            
        }
    

    ## Return
    list(OM = OM_update,
         EM = EM_update,
         ssb = ssb,
         fbar = fbar,
         rec = rec,
         catch = catch,
         EMoutput = EMoutput,
         AFoutput = AFoutput,
         adviceRules = adviceRules,
         msg = msg)

}



MSE_FixedConstraint <- function(OM,
                                constraint,
                                nYears,
                                AdviceYears = 1, ## These are to make sure the output is similar to F!=0
                                AdviceLag = 1,
                                initialAdvice = 0,
                                inherentImplementationError = FALSE,                  
                                ...){
    attach(getNamespace("multiStockassessment"))
    cat("Setting up the MSE\n")
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
    ## F=0, so selectivity does not matter

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
    adviceRules <- matrix(NA_character_,nYears+AdviceLag,nStocksOM+1)
    
    dimnames(ssb) <- dimnames(fbar) <- dimnames(rec) <- dimnames(catch) <- list(seq(do.call(max,lapply(OM,function(x)x$data$years)) + 1,len = nYears + AdviceLag),
                                                                                c(getStockNames(OM),"Total"),
                                                                                c("Advice","True","Estimate","Low","High"))
    dimnames(adviceRules) <- list(seq(do.call(max,lapply(OM,function(x)x$data$years)) + 1,len = nYears + AdviceLag),
                                  c(getStockNames(OM),"Total"))

    EMoutput <- vector("list",nYears+AdviceLag)
    AFoutput <- vector("list",nYears+AdviceLag)
    names(EMoutput) <- names(AFoutput) <- seq(do.call(max,lapply(OM,function(x)x$data$years)) + 1,len = nYears + AdviceLag)

##### Copy OM and EM for safe overwriting #####
    OM_update <- OM
    msg <- "OK"

##### Insert initialAdvice in results table #####
    if(!is.matrix(initialAdvice))
        initialAdvice <- matrix(initialAdvice,AdviceLag,length(OM),byrow=TRUE)
    for(s in seq_along(OM)){
        catch[seq_len(AdviceLag),s,"Advice"] <- initialAdvice[,s]
    }
    ## Total
    catch[seq_len(AdviceLag),"Total","Advice"] <- rowSums(initialAdvice)
    adviceRules[seq_len(AdviceLag),seq_along(OM)] <- "Initial advice"
    cat("Ready to start\n")
    if(AdviceLag > 0){
        cat("Update OM to handle advice lag\n")
        ## Update OM
        iy <- head(rownames(ssb),AdviceLag)
        OM_update <- try({addSimulatedYears(OM_update, 
                                            constraints = constraint,
                                            deterministicF = !inherentImplementationError,
                                            ## maxTrueF = maxTrueF,
                                            ## maxScaleF = maxScaleF,
                                            ...)})
        ssb[iy,,"True"] <- ssbtable(OM_update,addTotal=TRUE)[iy,seq(1,by=3,length.out=nStocksOM+1)]
        fbar[iy,,"True"] <- fbartable(OM_update,addTotal=TRUE)[iy,seq(1,by=3,length.out=nStocksOM+1)]
        rec[iy,,"True"] <- rectable(OM_update,addTotal=TRUE)[iy,seq(1,by=3,length.out=nStocksOM+1)]
        catch[iy,,"True"] <- catchtable(OM_update,addTotal=TRUE)[iy,seq(1,by=3,length.out=nStocksOM+1)]
    }
    cat("Starting assessment loop\n")      
    ## Run assessment loop
    for(i in seq(1,nYears-(AdviceYears-1), by = AdviceYears)){ # Index over assessment year
        yr <- rownames(ssb)[seq(i,len=AdviceYears)]
        yr_tac <- rownames(ssb)[seq(i+AdviceLag,len=AdviceYears)]
        
        cat("\n\n\nSimulation year",i,"\n")
        cat("\tAdvice year",yr_tac,"\n")
        cat("\tLast true data year",tail(rownames(fbartable(OM_update)),1),"\n")
        cat("\tAssessment year",yr[1],"\n")
        cat("\tCurrent SSB",ssb[tail(rownames(fbartable(OM_update)),1),,"True"],"\n")

        cat("\t\tNot doing an assessment for F=0...\n")
        cat("\t\tNo output to save...\n")
        
        cat("\t\tNo need for a forecast...\n")
        adviceRules[yr_tac,] <- "F=0"
        cat("\t\tSaving advice...\n")
        ## Save advice
        catch[yr_tac,,"Advice"] <- 0
        fbar[yr_tac,,"Advice"] <- 0
        ## ssb[yr_tac,,"Advice"] <- afFTab[[s]][yr_tac,sprintf("ssb:%s",tabLab[[s]])]
        ## rec[yr_tac,,"Advice"] <- afFTab[[s]][yr_tac,sprintf("rec:%s",tabLab[[s]])]
        ## Simulate next year (implement different rules for AdviceYears > 0)
        cat("\t\tUpdate OM...\n")
        capture.output(OM_update <- try({addSimulatedYears(OM_update, 
                                                           constraints =  constraint,
                                                           deterministicF = !inherentImplementationError,
                                                           ## maxTrueF = maxTrueF,
                                                           ## maxScaleF = maxScaleF,
                                                          ...)}))
        if(methods::is(OM_update,"try-error")){
            msg <- "Adding simulated year error"
            break;
        }
        cat("\t\tSaving output...\n")
        tab <- ssbtable(OM_update,addTotal=TRUE,returnList=TRUE)
        for(s in seq_along(tab))
            ssb[cAdd(yr,1),s,"True"] <- tab[[s]][cAdd(yr,1),1]
        tab <- fbartable(OM_update,addTotal=TRUE,returnList=TRUE)
        for(s in seq_along(tab))
            fbar[cAdd(yr,1),s,"True"] <- tab[[s]][cAdd(yr,1),1]
        tab <- rectable(OM_update,addTotal=TRUE,returnList=TRUE)
        for(s in seq_along(tab))
            rec[cAdd(yr,1),s,"True"] <- tab[[s]][cAdd(yr,1),1]
        tab <- catchtable(OM_update,addTotal=TRUE,returnList=TRUE)
        for(s in seq_along(tab))
            catch[cAdd(yr,1),s,"True"] <- tab[[s]][cAdd(yr,1),1]     
    }    

    ## Return
    list(OM = OM_update,
         EM = NA,
         ssb = ssb,
         fbar = fbar,
         rec = rec,
         catch = catch,
         EMoutput = EMoutput,
         AFoutput = AFoutput,
         adviceRules = adviceRules,
         msg = msg)

}
