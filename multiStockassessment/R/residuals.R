
##' @export
residuals.msam <- function(object, ... ) {
    oldObj <- attr(object,"m_obj")
    dat <- oldObj$env$data

    logobs <- do.call("c",lapply(dat$sam,function(x)x$logobs))
    stock <- do.call("c",lapply(as.list(1:length(dat$sam)),function(i)rep(i-1,length(dat$sam[[i]]$logobs))))
    indx <- do.call("c",lapply(dat$sam,function(x)seq_along(x$logobs)-1))
    year <- do.call("c",lapply(dat$sam,function(x)x$aux[,"year"]))
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
    residObj <- TMB::MakeADFun(dat,attr(object,"m_pl"),oldObj$env$map, random=ran, DLL="multiStockassessment", ...)
    resid <- TMB::oneStepPredict(residObj,"fake_obs","fake_keep",discrete = FALSE, trace = 1)
    resid
}
