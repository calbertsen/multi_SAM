
##' @export
residuals.msam <- function(object, discrete = FALSE, trace = TRUE, ...) {
    cat(sprintf("One-observation-ahead residuals for %d stock%s. Total number of observations: %d\n",    length(object),ifelse(length(object)>1,"s",""),nobs(object)))
    oldObj <- attr(object,"m_obj")
    dat <- oldObj$env$data
    
    logobs <- do.call("c",lapply(dat$sam,function(x)x$logobs))
    stock <- do.call("c",lapply(as.list(1:length(dat$sam)),function(i)rep(i-1,length(dat$sam[[i]]$logobs))))
    indx <- do.call("c",lapply(dat$sam,function(x)seq_along(x$logobs)-1))
    year <- do.call("c",lapply(dat$sam,function(x)x$aux[,"year"]))
    fleetNum <- do.call("c",lapply(dat$sam,function(x)x$aux[,"fleet"]))
    fleet <- do.call("c",lapply(as.list(1:length(dat$sam)),
                                function(i){
                                    prefix <- gsub(" ","0",
                                                   formatC(i,width=nchar(as.character(length(dat$sam)))))
                                    paste0(prefix,"_",dat$sam[[i]]$aux[,"fleet"])
                                }))
    age <- do.call("c",lapply(dat$sam,function(x)x$aux[,"age"]))
    sortVals <- order(year,fleet,age)
    dat$fake_obs <- logobs[sortVals]
    dat$fake_stock <- stock[sortVals]
    dat$fake_indx <- indx[sortVals]
    dat$doResiduals <- 1

    ran <- c("logN", "logF", "missing")
    residObj <- TMB::MakeADFun(dat,attr(object,"m_pl"),oldObj$env$map, random=ran, DLL="multiStockassessment")
    
    resid <- TMB::oneStepPredict(residObj,"fake_obs","fake_keep",
                                 discrete = discrete,
                                 trace = trace, ...)
    resid$residual[is.na(resid$observation)] <- NA
    cat("One-observation-ahead residuals. Done\n")
    ret <- cbind(year = year[sortVals], fleet = as.numeric(factor(fleet))[sortVals], age = age[sortVals], resid, stock = stock[sortVals], stockFleet = fleetNum[sortVals])
    stockNames <- multiStockassessment:::getStockNames(object)
    fleetNamesList <- lapply(as.list(1:length(object)),
                             function(i){
                                 paste(stockNames[i],attributes(object[[i]]$data)$fleetNames,sep=" - ")
                                 })       
    attr(ret,"fleetNames") <- unlist(fleetNamesList)
    class(ret)<- c("msamres","samres","data.frame")
    ret
}


