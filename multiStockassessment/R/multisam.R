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
##' @return A list of class samset and samsetfit
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
##' @importFrom stats nlminb
##' @importFrom TMB MakeADFun
##' @export
multisam.fit <- function(x,corStructure,usePartialCors=TRUE,newtonsteps=3,lower=NULL,upper=NULL,...){
    requireNamespace("TMB")
    if(class(x) != "samset")
        stop("x must be a samset.")
    
    dat0 <- collect_data(x)
    dat <- list(sam = dat0,
                usePartialCor = as.integer(usePartialCors),
                maxYearAll = max(unlist(lapply(dat0,function(dd)dd$years))),
                minYearAll = min(unlist(lapply(dat0,function(dd)dd$years))),
                maxAgeAll = max(unlist(lapply(dat0,function(dd)dd$maxAge))),
                minAgeAll = min(unlist(lapply(dat0,function(dd)dd$minAge))),
                cons = boolMat2ConstraintList(corStructure)
                )
    pars <- collect_pars(x)
    pars$RE <- rep(0,sum(lower.tri(corStructure)))
    if(is.null(lower))
        lower <- getLowerBounds(pars)
    if(is.null(upper))
        upper <- getUpperBounds(pars)
    ran <- c("logN", "logF", "missing")
    obj <- TMB::MakeADFun(dat, pars, random=ran, DLL="multiStockassessment", ...)
    indx <- which(obj$gr()[names(obj$par) %in% "RE"] == 0)
    mvals <- 1:sum(lower.tri(corStructure))
    mvals[indx] <- NA
    map <- list(RE = factor(mvals))
    obj <- TMB::MakeADFun(dat, pars, map, random=ran, DLL="multiStockassessment", ...)
    
    lower2<-rep(-Inf,length(obj$par))
    upper2<-rep(Inf,length(obj$par))
    for(nn in names(lower)) lower2[names(obj$par)==nn]=lower[[nn]]
    for(nn in names(upper)) upper2[names(obj$par)==nn]=upper[[nn]]
    
    opt <- stats::nlminb(obj$par, obj$fn,obj$gr ,control=list(trace=1, eval.max=2000, iter.max=1000),lower=lower2,upper=upper2)

    ## for(i in seq_len(newtonsteps)) { # Take a few extra newton steps 
    ##     g <- as.numeric( obj$gr(opt$par) )
    ##     h <- optimHess(opt$par, obj$fn, obj$gr)
    ##     opt$par <- opt$par - solve(h, g)
    ##     opt$objective <- obj$fn(opt$par)
    ## }
    ## rep <- obj$report()
    ## sdrep <- sdreport(obj,opt$par)

    ##                                     # Last two states
    ## idx <- c(which(names(sdrep$value)=="lastLogN"),which(names(sdrep$value)=="lastLogF"))
    ## sdrep$estY <- sdrep$value[idx]
    ## sdrep$covY <- sdrep$cov[idx,idx]

    ## idx <- c(which(names(sdrep$value)=="beforeLastLogN"),which(names(sdrep$value)=="beforeLastLogF"))
    ## sdrep$estYm1 <- sdrep$value[idx]
    ## sdrep$covYm1 <- sdrep$cov[idx,idx]

    ## pl <- as.list(sdrep,"Est")
    ## plsd <- as.list(sdrep,"Std")

    ## sdrep$cov<-NULL # save memory

    ## ## Split everything to sam classes

    res <- list(obj=obj,opt=opt)
    class(res) <- c("msamfit")
    return(res)
}
