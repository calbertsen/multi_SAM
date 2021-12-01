
##' @export
residuals.msam <- function(object, discrete = FALSE, trace = TRUE, ...) {
    oldObj <- attr(object,"m_obj")
    dat <- oldObj$env$data
    logobsList <- lapply(dat$sam,function(x)x$logobs)
    auxList <- lapply(dat$sam,function(x)x$aux)
    if(dat$sharedObs$hasSharedObs){
        logobsList <- c(logobsList, SharedObservations=list(dat$sharedObs$logobs))
        auxList <- c(auxList, SharedObservations=list(dat$sharedObs$aux))
    }
    logobs <- do.call("c",logobsList)
    stock <- rep(seq_along(logobsList)-1, times = sapply(logobsList,length))
    indx <- do.call("c",lapply(logobsList,function(x)seq_along(x)-1))
    year <- do.call("c",lapply(auxList,function(x)x[,"year"]))
    fleetNum <- do.call("c",lapply(auxList,function(x)x[,"fleet"]))
    fleet <- do.call("c",lapply(as.list(seq_along(auxList)),
                                function(i){
                                    prefix <- gsub(" ","0",
                                                   formatC(i,width=nchar(as.character(length(auxList)))))
                                    paste0(prefix,"_",auxList[[i]][,"fleet"])
                                }))
    age <- do.call("c",lapply(auxList,function(x)x[,"age"]))
    sortVals <- order(year,fleet,age)
    dat$fake_obs <- logobs[sortVals]
    dat$fake_stock <- stock[sortVals]
    dat$fake_indx <- indx[sortVals]
    dat$doResiduals <- 1

    cat(sprintf("One-observation-ahead residuals for %d stock%s%s. Total number of observations: %d\n",    length(object),ifelse(length(object)>1,"s",""),ifelse(!is.null(dat$sharedObs$logobs)," with shared observations",""),nobs(object)))
   
    
    ran <- c("logN", "logF", "missing", "shared_logFscale")
    residObj <- TMB::MakeADFun(dat,attr(object,"m_pl"),oldObj$env$map, random=ran, DLL="multiStockassessment")
    
    resid <- TMB::oneStepPredict(residObj,"fake_obs","fake_keep",
                                 discrete = discrete,
                                 trace = trace, ...)
    resid$residual[is.na(resid$observation)] <- NA
    cat("One-observation-ahead residuals. Done\n")
    ret <- cbind(year = year[sortVals], fleet = as.numeric(factor(fleet))[sortVals], age = age[sortVals], resid, stock = stock[sortVals], stockFleet = fleetNum[sortVals])
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


