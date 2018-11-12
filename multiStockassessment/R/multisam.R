##' Fit multiple SAM models with correlated survival processes.
##'
##' @title Fit multiple SAM models with correlations
##' @param x samset from the stockassessment package
##' @param corStructure symmetric boolean matrix. True if a (partial) correlation in survival should be estimated between the corresponding age/stock combination
##' @param usePartialCors if TRUE corStructure describes the partial correlations. If FALSE corStructure describes correlations.
##' @param newtonsteps As for stockassessment::sam.fit
##' @param lower As for stockassessment::sam.fit
##' @param upper As for stockassessment::sam.fit
##' @param ... Additional arguments passed to TMB::MakeADFun
##' @return A list of class msam and samset
##' @author Christoffer Moesgaard Albertsen
##' @examples
##' \donttest{
##' if(require(stockassessment)){
##'   data(nscodData)
##'   data(nscodConf)
##'   data(nscodParameters)
##'   fit <- sam.fit(nscodData, nscodConf, nscodParameters)
##'   fits <- c(fit)
##'   obj <- multisam.fit(fits)
##' }
##' }
##' @importFrom stats nlminb optimHess
##' @importFrom methods is
##' @importFrom TMB MakeADFun sdreport
##' @export
multisam.fit <- function(x,corStructure,usePartialCors=TRUE,newtonsteps=3,lower=NULL,upper=NULL,...){
    ## Check input
    requireNamespace("TMB")
    if(!methods::is(x,"samset"))
        stop("x must be a samset.")

    ## Prepare data for TMB
    dat0 <- collect_data(x)
    dat <- list(sam = dat0,
                usePartialCor = as.integer(usePartialCors),
                maxYearAll = max(unlist(lapply(dat0,function(dd)dd$years))),
                minYearAll = min(unlist(lapply(dat0,function(dd)dd$years))),
                maxAgeAll = max(unlist(lapply(dat0,function(dd)dd$maxAge))),
                minAgeAll = min(unlist(lapply(dat0,function(dd)dd$minAge))),
                cons = boolMat2ConstraintList(corStructure),
                fake_obs = numeric(0),
                fake_stock = numeric(0),
                fake_indx = numeric(0),
                doResiduals = as.integer(FALSE)
                )

    ## Prepare parameters for TMB
    pars <- collect_pars(x)
    pars$RE <- rep(0,sum(lower.tri(corStructure)))

    ## Create initial TMB Object
    ran <- c("logN", "logF", "missing")
    obj <- TMB::MakeADFun(dat, pars, random=ran, DLL="multiStockassessment", ...)

    ## Check for unused correlation parameters
    indx <- which(obj$gr()[names(obj$par) %in% "RE"] == 0)
    mvals <- 1:sum(lower.tri(corStructure))
    mvals[indx] <- NA
    map <- list(RE = factor(mvals))

    ## Create TMB Object
    obj <- TMB::MakeADFun(dat, pars, map, random=ran, DLL="multiStockassessment", ...)

    pars2 <- pars
    pars2$RE <- obj$par[names(obj$par)=="RE"]
    ## Fit object
   if(is.null(lower))
        lower <- getLowerBounds(pars2)
    if(is.null(upper))
        upper <- getUpperBounds(pars2)

    lower2 <- rep(-Inf,length(obj$par))
    upper2 <- rep(Inf,length(obj$par))
    for(nn in names(lower)) lower2[names(obj$par)==nn]=lower[[nn]]
    for(nn in names(upper)) upper2[names(obj$par)==nn]=upper[[nn]]
    
    opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
                         control=list(trace=1,
                                      eval.max=2000,
                                      iter.max=1000),
                         lower=lower2,
                         upper=upper2)

    for(i in seq_len(newtonsteps)) { # Take a few extra newton steps
        atLBound <- (opt$par < (lower2 + sqrt(.Machine$double.eps)))
        atUBound <- (upper2 < (opt$par + sqrt(.Machine$double.eps)))
        atBound <- atLBound | atUBound
        g <- as.numeric( obj$gr(opt$par) )
        h <- stats::optimHess(opt$par, obj$fn, obj$gr)
        opt$par[!atBound] <- opt$par[!atBound]- solve(h[!atBound,!atBound], g[!atBound])
        opt$par[atBound] <- (atLBound * lower2 + atUBound * upper2)[atBound]
        opt$objective <- obj$fn(opt$par)
    }

    ## Get report and sdreport
    rep <- obj$report(obj$env$last.par.best)
    sdrep <- TMB::sdreport(obj,opt$par)
    ssdrep <- summary(sdrep)

    ## Do as in stockassessment package
    # Last two states
    idx <- c(which(names(sdrep$value)=="lastLogN"),which(names(sdrep$value)=="lastLogF"))
    sdrep$estY <- sdrep$value[idx]
    sdrep$covY <- sdrep$cov[idx,idx]

    idx <- c(which(names(sdrep$value)=="beforeLastLogN"),which(names(sdrep$value)=="beforeLastLogF"))
    sdrep$estYm1 <- sdrep$value[idx]
    sdrep$covYm1 <- sdrep$cov[idx,idx]

    pl <- as.list(sdrep,"Est")
    plsd <- as.list(sdrep,"Std")
    
    sdrep$cov<-NULL # save memory
    

    ## Prepare result
    res <- x

    ##sdrep

    ## pl

    ## plsd

    ## rep

    ## Add multi SAM info
    attr(res,"m_obj") <- obj
    attr(res,"m_opt") <- opt
    attr(res,"m_rep") <- rep
    attr(res,"m_sdrep") <- sdrep
    attr(res,"m_pl") <- pl
    attr(res,"m_plsd") <- plsd
    attr(res,"m_data") <- dat
    attr(res,"m_low") <- lower2
    attr(res,"m_high") <- upper2
    attr(res,"corStructure") <- corStructure
    attr(res,"partialCors") <- usePartialCors
    attr(res,"correlationParameters") <- ssdrep[rownames(ssdrep)=="RE",]
    class(res) <- c("msam","samset")
    return(res)
}
