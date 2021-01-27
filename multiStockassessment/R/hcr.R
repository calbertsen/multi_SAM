
##' @rdname hcr
##' @export
##' 
hcr <- function(fit,
                ...){
    UseMethod("hcr")
}

##' @rdname hcr
##' @method hcr sam
##' @importFrom stockassessment hcr
##' @export
hcr.sam <- function(fit,...){
    stockassessment::hcr(fit, ...)
}
                

##' @rdname hcr
##' @method hcr msam
##' @export
hcr.msam <- function(fit,                     
                     nYears = 20,
                     Ftarget,
                     Btrigger,
                     Flim = 0,
                     Blim = 0,
                     Flow = 0,
                     Blow = 0,
                     nosim = 10000,
                     ave.years = lapply(fit,function(x)max(x$data$years)+(-4:0)),
                     rec.years = lapply(fit,function(x)numeric(0)), #max(fit$data$years)+(-9:0),
                     preForecast = list(),
                     ...){

    nStocks <- length(fit)

    nYears <- rep(nYears, len = nStocks)
    Ftarget <- rep(Ftarget, len = nStocks)
    Btrigger <- rep(Btrigger, len = nStocks)
    Flim <- rep(Flim, len = nStocks)
    Blim <- rep(Blim, len = nStocks)
    Flow <- rep(Flow, len = nStocks)
    Blow <- rep(Blow, len = nStocks)
    hcrConf <- split(cbind(Ftarget,Flim,Flow,Blim,Blow,Btrigger), seq_along(fit))
    
    hcrArg <- vector("list",nStocks)
    if(length(preForecast) > 0){
        preForecast <- lapply(preForecast, targetToList, nStocks = nStocks)
        preYears <- sapply(preForecast[[1]], length)
        preForecast <- lapply(preForecast, function(xx)
            lapply(seq_along(xx),
                   function(i) c(xx[[i]], rep(NA, nYears[i])))
            )
    }else{
        preYears <- numeric(nStocks)
    }
    hcr <- lapply(seq_along(nYears), function(i) rep(c(NA,1),c(preYears[i],nYears[i])))
    
   args <- c(list(fit = fit,
                 hcr = hcr,
                 hcrConf = hcrConf,
                 nosim = nosim,
                 ave.years = ave.years,
                 rec.years = rec.years),
              preForecast,
              list(...)
              )
    f <- do.call("modelforecast", args)

    hcrFuns <- sapply(seq_len(nStocks), function(i){ function(ssb) .Call("hcrR", ssb, c(Ftarget[i], Flim[i], Flow[i], Blim[i], Blow[i], Btrigger[i]), PACKAGE = "stockassessment") })
    
    r <- list(Ftarget = Ftarget,
              Flim = Flim,
              Flow = Flow,
              Blim = Blim,
              Blow = Blow,
              Btrigger = Btrigger,
              hcr = hcrFuns,
              forecast = f)
    attr(r,"fit") <- fit
    class(r) <- "msam_hcr"
    r
}

