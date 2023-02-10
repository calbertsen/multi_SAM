##' @method stochasticReferencepoints msam
##' @importFrom stockassessment stochasticReferencepoints
##' @export
stochasticReferencepoints.msam <- function(fit,
                                          referencepoints,
                                          method = "Q0.5",
                                          catchType = "catch",
                                          nYears = 300,
                                          Frange = c(0,2),
                                          nosim = 1000,
                                          aveYears = c(), #max(fit$data$years)+(-9:0),
                                          selYears = max(fit$data$years),
                                          newton.control = list(),
                                          seed = .timeToSeed(),
                                          formula = ~ibc(F,5),
                                          nosim_ci = 200,
                                          derivedSummarizer = NA,
                                          nTail = 1,
                                          constraint = "F=%f",
                                          deterministicF = TRUE,
                                          ...){
    warning("The function is not fully implemented yet")
    oldSeed <- NULL
    if(exists(".Random.seed"))
        oldSeed <- .Random.seed
    on.exit(set.seed(oldSeed))
    set.seed(seed)

    .refpointSMethodParser <- getFromNamespace(".refpointSMethodParser","stockassessment")
    .refpointMerger <- getFromNamespace(".refpointMerger","stockassessment")
    .refpointParser <- getFromNamespace(".refpointParser","stockassessment")
    .refpointCheckRecruitment <- getFromNamespace(".refpointCheckRecruitment","stockassessment")

    MT <- .refpointSMethodParser(method, formula = formula, derivedSummarizer=derivedSummarizer)

    
    catchType <- pmatch(catchType,c("catch","landing","discard"))-1
    if(is.na(catchType))
        stop("Invalid catch type")

    aveYearsIn <- aveYears
    ## aveYears <- match(aveYears, fit$data$years) - 1
    ## if(any(is.na(aveYears)))
    ##     stop("aveYears has years without data.")

    selYearsIn <- selYears
    ## selYears <- match(selYears, fit$data$years) - 1
    ## if(any(is.na(selYears)))
    ##     stop("selYears has years without data.")

    if(!all(Frange >= 0) && Frange[1] < Frange[2] && length(Frange) ==2)
        stop("Wrong Frange")
    if(!nosim > 0)
        stop("nosim must be a positive integer")

    ## Get RPs for best fit
    rpArgs <- Reduce(.refpointMerger,
                     lapply(referencepoints, .refpointParser, cutoff = 0.1),
                     list())
    invisible(lapply(fit, function(ff) lapply(rpArgs,.refpointCheckRecruitment,fit=ff)))

    ## !!! SHOULD BE MSAM SPECIFIC
    doSim <- .getDoSim(logf1=tail(log(fbartable(fit)[,1]),1),
                       fit=fit, nYears = nYears, aveYears = aveYears, selYears = selYears, pl = fit$pl, constraint=constraint,deterministicF=deterministicF,...)

}
