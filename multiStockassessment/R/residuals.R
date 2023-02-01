
##' @export
residuals.msam <- function(object, discrete = FALSE, trace = TRUE, ...) {
    oldObj <- attr(object,"m_obj")
    dat <- oldObj$env$data
    logobsList <- lapply(dat$sam,function(x)x$logobs)
    auxList <- lapply(dat$sam,function(x){
        cbind(x$aux,fleetType = x$fleetTypes[x$aux[,"fleet"]])
    })
    auxDList <- lapply(dat$sam,function(x){
        x$auxData
    })
    if(dat$sharedObs$hasSharedObs){
        logobsList <- c(logobsList, SharedObservations=list(dat$sharedObs$logobs))
        tmp <- cbind(dat$sharedObs$aux, fleetType=dat$sharedObs$fleetTypes[dat$sharedObs$aux[,"fleet"]])
        auxList <- c(auxList, SharedObservations=list(tmp))
        auxDList <- c(auxDList, SharedObservations=list(dat$sharedObs$auxData))
    }
    logobs <- do.call("c",logobsList)
    stock <- rep(seq_along(logobsList)-1, times = sapply(logobsList,length))    
    indx <- do.call("c",lapply(logobsList,function(x)seq_along(x)-1))
    year <- do.call("c",lapply(auxList,function(x)x[,"year"]))
    dpart <- do.call("c",lapply(auxDList,function(x){
        if(ncol(x) < 5)
            return(rep(0,nrow(x)))
        x[,5]
        }))
    fleetType <- do.call("c",lapply(auxList,function(x)x[,"fleetType"]))
    fleet <- do.call("c",lapply(as.list(seq_along(auxList)),
                                function(i){
                                    prefix <- gsub(" ","0",
                                                   formatC(i,width=nchar(as.character(length(auxList)))))
                                    paste0(prefix,"_",auxList[[i]][,"fleet"])
                                }))
    age <- do.call("c",lapply(auxList,function(x)x[,"age"]))
    xtr <- ifelse(fleetType >= 80, dpart,0)
    kp <- which(!is.na(logobs) & fleetType < 80)
    sortVals <- order(year[kp],stock[kp],fleet[kp],age[kp], xtr[kp])
    dat$fake_obs <- logobs[kp][sortVals]
    dat$fake_stock <- stock[kp][sortVals]
    dat$fake_indx <- indx[kp][sortVals]
    dat$fake_includeOthers <- 0
    dat$doResiduals <- 1
    ##rSubset <- which(!is.na(dat$fake_obs))

    cat(sprintf("One-observation-ahead residuals for %d stock%s%s. Total number of observations: %d\n",    length(object),ifelse(length(object)>1,"s",""),ifelse(!is.null(dat$sharedObs$logobs)," with shared observations",""),nobs(object)))
   
    
    ran <- c("logN", "logF", "missing","logSW", "logCW", "logitMO", "logNM", "logP","logitFseason", "shared_logFscale","shared_missingObs", "logGst", "logGtrip")
    residObj1 <- TMB::MakeADFun(dat,attr(object,"m_pl"),oldObj$env$map, random=ran, DLL="multiStockassessment")
    residObj1$fn()
    resid1 <- TMB::oneStepPredict(residObj1,"fake_obs","fake_keep",
                                 discrete = discrete,
                                 trace = trace,
                                 reverse = TRUE,
                                 method = "oneStepGaussianOffMode")## ,
                                 ## ...)

    ## Part 2
    kp2 <- which(!is.na(logobs) & fleetType >= 80)
    if(length(kp2) > 0){
        sortVals2 <- order(year[kp2],stock[kp2],fleet[kp2],age[kp2], xtr[kp2])
        dat$fake_obs <- logobs[kp2][sortVals2]
        dat$fake_stock <- stock[kp2][sortVals2]
        dat$fake_indx <- indx[kp2][sortVals2]
        dat$fake_includeOthers <- 1
        dat$doResiduals <- 1
        ##rSubset <- which(!is.na(dat$fake_obs))

        residObj2 <- TMB::MakeADFun(dat,attr(object,"m_pl"),oldObj$env$map, random=ran, DLL="multiStockassessment")
        residObj2$fn()
        resid2 <- TMB::oneStepPredict(residObj2,"fake_obs","fake_keep",
                                      discrete = FALSE,
                                      trace = trace,
                                      reverse = FALSE,
                                      method = "oneStepGaussianOffMode")
    }


    ## Should check for subset!.
    residCollect <- as.data.frame(lapply(colnames(resid1), function(nm){
        rep(NA_real_, length(logobs))
    }))
    colnames(residCollect) <- colnames(resid1)
    residCollect[kp[sortVals],] <- resid1
    if(length(kp2) > 0)
        residCollect[kp2[sortVals2],] <- resid2
    
    ##resid$residual[is.na(resid$observation)] <- NA
    cat("One-observation-ahead residuals. Done\n")
    ret <- cbind(do.call("rbind",lapply(auxList,as.data.frame)), residCollect, stock = stock, do.call("rbind",lapply(auxDList,as.data.frame)))
    stockNames <- getStockNames(object)
    fleetNamesList <- lapply(as.list(1:length(object)),
                             function(i){
                                 paste(stockNames[i],attributes(object[[i]]$data)$fleetNames,sep=" - ")
                             })
    if(dat$sharedObs$hasSharedObs){
        fleetNamesList <- c(fleetNamesList, list(paste("Shared observation",attr(attr(object,"m_data")$sharedObs,"fleetNames"),sep=" - ")))
    }
    attr(ret,"fleetNames") <- unlist(fleetNamesList)
    class(ret)<- c("msamres","samres","data.frame")
    ret
}


