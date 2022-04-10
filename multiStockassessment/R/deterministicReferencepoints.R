##' Deterministic reference points
##'
##' @param fit 
##' @param referencepoints 
##' @param catchType 
##' @param nYears 
##' @param Fsequence 
##' @param aveYears 
##' @param selYears 
##' @param biasCorrect 
##' @param newton.control 
##' @param ... 
##' @return 
##' @author Christoffer Moesgaard Albertsen
##' @importFrom stockassessment deterministicReferencepoints
##' @method deterministicReferencepoints msam
##' @export
deterministicReferencepoints.msam <- function(fit,
                                              referencepoints,
                                              catchType = "landing",
                                              nYears = 300,
                                              Fsequence = seq(0,2, len = 50),
                                              aveYears = lapply(fit,function(x)max(x$data$years)+(-9:0)),
                                              selYears = lapply(fit,function(x)max(x$data$years)),
                                              biasCorrect = FALSE,
                                              newton.control = list(),
                                              ...){

    .checkInput <- function(x, checkFun = identity){
        if(!is.list(x))
            x <- replicate(nStocks, list(x))
        if(!length(x) == nStocks)
            stop("Wrong input")
        lapply(x, checkFun)
    }
    nStocks <- length(fit)

    referencepoints <- .checkInput(referencepoints)
    catchType <- .checkInput(catchType, function(y){
        v <- pmatch(y, c("catch","landing","discard"))
        if(length(v) != 1 || any(is.na(v))) stop("Invalid catchType")
        v})
    nYears <- .checkInput(nYears)
    Fsequence <- .checkInput(Fsequence, function(fs){
        if(!all(diff(fs) > 0) || !all(fs >= 0))
            stop("Values of Fsequence must be positive and increasing.")
        if(!isTRUE(all.equal(fs[1],0, check.attributes = FALSE, use.names = FALSE)))
            warning("The first value of Fsequence should be 0.")
        fs
    })
    aveYears <- .checkInput(aveYears)
    aveYearsIn <- aveYears
    for(i in seq_along(aveYears))
        aveYears[[i]] <- match(aveYears[[i]], fit[[i]]$data$years)-1
    selYears <- .checkInput(selYears)
    selYearsIn <- selYears
    for(i in seq_along(selYears))
        selYears[[i]] <- match(selYears[[i]], fit[[i]]$data$years)-1

    .refpointEnum <- getFromNamespace(".refpointEnum","stockassessment")
    .refpointNames <- getFromNamespace(".refpointNames","stockassessment")
    .refpointParserS <- getFromNamespace(".refpointParser","stockassessment")
    .refpointMergerS <- getFromNamespace(".refpointMerger","stockassessment")
    .refpointStartingValueS <- getFromNamespace(".refpointStartingValue","stockassessment")
    .refpointOutput <- getFromNamespace(".refpointOutput","stockassessment")

    fayTabs <- faytable(fit,returnList=TRUE)
    fbarTabs <- fbartable(fit,returnList=TRUE)

    .makeOneRpArgs <- function(i){
        ## Parse input reference points
        rpArgs <- Reduce(.refpointMergerS,
                         lapply(referencepoints[[i]], .refpointParserS, nYears = nYears[[i]], aveYears = aveYears[[i]], selYears = selYears[[i]], logCustomSel = numeric(0), catchType = catchType[[i]] - 1),
                         list())
        ## Add starting values    
        rpArgs <- lapply(seq_along(rpArgs), function(j).refpointStartingValueS(rpArgs[[j]], fit = fit[[i]], Fsequence = Fsequence[[i]], fay = fayTabs[[i]], fbar = fbarTabs[[i]][,1]))
        ## Add Fsequence
        rp0 <- list(rpType = -1,
                    xVal = log(Fsequence[[i]]),
                    nYears = nYears[[i]],
                    aveYears = aveYears[[i]],
                    selYears = selYears[[i]],
                    logCustomSel = numeric(0),
                    catchType = catchType[[i]] - 1,
                    logF0 = log(Fsequence[[i]]))
        r <- c(list(rp0), rpArgs)
        attr(r,"newton_config") <- newton.control
        r
    }
    rpArgs <- lapply(seq_len(nStocks), .makeOneRpArgs)
    
    obj0 <- attr(fit,"m_obj")
    argsIn <- as.list(obj0$env)[methods::formalArgs(TMB::MakeADFun)[methods::formalArgs(TMB::MakeADFun) != "..."]]
    argsIn$silent <- obj0$env$silent
    argsIn$parameters <- attr(fit,"m_pl")
    for(i in seq_along(fit))
        argsIn$data$sam[[i]]$referencepoints <- rpArgs[[i]]
    
    args <- argsIn

    objSDR <- do.call(TMB::MakeADFun, args)
    objSDR$fn(attr(fit,"m_opt")$par)
    sdr <- TMB::sdreport(objSDR, objSDR$par, attr(fit,"m_opt")$he,
                         bias.correct= biasCorrect,
                         skip.delta.method = biasCorrect,
                         bias.correct.control = list(sd = TRUE,
                                                     split = objSDR$env$ADreportIndex()[grepl("referencepoint_[[:digit:]]+_.+",names(objSDR$env$ADreportIndex()))]
                                                     ))
    ssdr <- summary(sdr)
    ssdr
    res <- vector("list",nStocks)
    ii <- grep("referencepoint_",rownames(ssdr))    
    ssdrReduced <- ssdr[ii,]
    ssdrRL <- split(ssdrReduced, gsub("(SAM_[[:digit:]]+)(_referencepoint.+)","\\1",rownames(ssdrReduced)))
    for(i in seq_along(res)){
        ssdrR <- ssdr[grepl(sprintf("SAM_%d_referencepoint",i-1),rownames(ssdr)),]
        if(nrow(ssdrR) > 0){            
            res[[i]] <- .refpointOutput(ssdrR,rpArgs[[i]][-1], fit[[i]], biasCorrect, aveYearsIn[[i]], selYearsIn[[i]], Fsequence[[i]], referencepoints[[1]])
        }
    }
    names(res) <- getStockNames(fit)
    class(res) <- "msam_referencepoints"
    res
}
