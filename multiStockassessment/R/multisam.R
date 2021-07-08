##' Fit multiple SAM models with correlated survival processes.
##'
##' @title Fit multiple SAM models with correlations
##' @param x samset from the stockassessment package
##' @param formula formula for covariance matrix covariates (See Details)
##' @param corStructure symmetric boolean matrix. True if a (partial) correlation in survival should be fixed to zero between the corresponding age/stock combination
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
##' @details
##' Function to fit a multi-stock SAM model (Albertsen et al., 2018).
##' @section Covariance covariates
##' The covariance matrix for the abundance (N) process can be described through a formula, similar to linear models. As covariates, the formula can include 'Age', 'Stock', 'Index', and the name of any matrix in the data set of the original SAM fits with the same number of columns as the number of age groups (e.g., stockMeanWeight). In the latter case, column means are used.
##'
##' To fit the models from Albertsen et al. (2018), the formula '~ factor(Index)' should be used in combination with 'suggestCorStructure'. The default values will fit independent models.
##' 
##' @importFrom stats nlminb optimHess
##' @importFrom methods is
##' @importFrom TMB MakeADFun sdreport
##' @references
##' Albertsen, C. M., Nielsen, A. and Thygesen, U. H. (2018) Connecting single-stock assessment models through correlated survival. ICES Journal of Marine Science, 75(1), 235-244. doi: 10.1093/icesjms/fsx114
##' @export
multisam.fit <- function(x,
                         formula = ~-1,
                         corStructure = suggestCorStructure(x,nAgeClose = Inf,noCorInArea=FALSE),
                         usePartialCors=TRUE,
                         newtonsteps=3,
                         rm.unidentified=FALSE,
                         nlminb.control = list(trace=1, eval.max=20000, iter.max=20000),
                         lower=NULL,
                         upper=NULL,
                         starting=NULL,
                         community_formula = ~-1,
                         community_type = 1,
                         shared_data = NULL,
                         shared_conf = NULL,
                         shared_keys = character(0),
                         shared_selectivity = FALSE,
                         shared_stockrecruitment = FALSE,
                         shared_oneFScaleSd = FALSE,
                         shared_initN = FALSE,
                         initN = 0,
                         initF = FALSE,
                         ...){
    ## Check input
    requireNamespace("TMB")
    if(!methods::is(x,"samset"))
        stop("x must be a samset.")

    if(is.null(shared_data))
        shared_data <- list(hasSharedObs = 0L)
    if(is.null(shared_conf))
        shared_conf <- list()
    
    ## Prepare data for TMB
    dat0 <- collect_data(x)
    dat <- list(sam = dat0,
                sharedObs = shared_data,
                usePartialCor = as.integer(usePartialCors),
                maxYearAll = as.integer(max(unlist(lapply(dat0,function(dd)dd$years)))),
                minYearAll = as.integer(min(unlist(lapply(dat0,function(dd)dd$years)))),
                maxAgeAll = as.integer(max(unlist(lapply(dat0,function(dd)dd$maxAge)))),
                minAgeAll = as.integer(min(unlist(lapply(dat0,function(dd)dd$minAge)))),
                cons = boolMat2ConstraintList(corStructure),
                X = build_cov_design(formula, x),
                XCon = build_cov_design(community_formula,  x, type = community_type),
                ConConf = as.integer(community_type),
                fake_obs = numeric(0),
                fake_stock = integer(0),
                fake_indx = integer(0),
                doResiduals = as.integer(FALSE)
                )

    ## Prepare parameters for TMB
    pars <- collect_pars(x)
    pars$RE <- numeric(ncol(dat$X)) ##rep(0,sum(lower.tri(corStructure)))
    pars$betaCon <- numeric(ncol(dat$XCon))

    pars$shared_logSdObs <- numeric(length(dat$sharedObs$fleetTypes))
    if(shared_selectivity){
        lfsDim <- local({xx <- attr(pars$logF,"cdim");  xx[1] <- 0; xx})
        ## newFdimR <- attr(pars$logF,"rdim")
        ## newFdimR[-1] <- 0
        ## newFdimC <- attr(pars$logF,"cdim")
        ## newFdimC[-1] <- 0
        ## pars$logF <- combineParameter(lapply(seq_along(newFdimR), function(i) matrix(0, newFdimR[i], newFdimC[i])))
        initFdim <- local({xx <- attr(pars$logF,"rdim") * (initF>0);  xx[-1] <- 1 * (initF>0); xx})
    }else{
        lfsDim <- attr(pars$logF,"cdim") * 0
        initFdim <- attr(pars$logF,"rdim") * 0
    }
    pars$shared_logFscale <- combineVectors(lapply(lfsDim,numeric))
    ## pars$shared_lfsMean <- numeric(ifelse(shared_selectivity,length(dat$sam)-1,0))
    pars$shared_lfsSd <- numeric(ifelse(shared_selectivity,length(dat$sam)-1,0))
    ## pars$shared_lfsRho <- 2 + numeric(ifelse(shared_selectivity,length(dat$sam)-1,0))
    pars$shared_logitMissingVulnerability <- numeric(sum(is.na(dat$sharedObs$keyFleetStock)))
    
    if(initN == 0){
        pars$initLogN <- combineVectors(lapply(attr(pars$logN,"rdim") * 0,numeric))
    }else if(initN == 1){
        pars$initLogN <- combineVectors(lapply(attr(pars$logN,"rdim") * 1,numeric))
    }else if(initN == 2){
        pars$initLogN <- combineVectors(lapply(rep(1,length(dat$sam)),numeric)) #+ c(11,11,11)
    }
    ##pars$initLogN <- combineVectors(lapply(attr(pars$logN,"rdim") * ifelse(initN,1,0),numeric))
    pars$initLogF <- combineVectors(lapply(initFdim,numeric))
    
    ## Prepare map for TMB
    map0 <- collect_maps(x)

    if(shared_oneFScaleSd){
        map0$shared_lfsSd <- factor(rep(1,length(pars$shared_lfsSd)))
        ## map0$shared_lfsRho <- factor(rep(1,length(pars$shared_lfsRho)))
    }
    if(shared_initN && initN == 2){
        map0$initLogN <- factor(rep(1,length(pars$initLogN)))
    }else if(initN == 2){
        ## map0$initLogN <- factor(c(1,NA,NA))#rep(NA,3))
    }
    
    
    ## Share keys
    keyName <- c("keyLogFpar","keyQpow","keyVarF","keyVarLogN","keyVarObs","keyCorObs")
    parName <- c("logFpar","logQpow","logSdLogFsta","logSdLogN","logSdLogObs","transfIRARdist")
    for(nm in shared_keys)
        map0[[parName[match(nm,keyName)]]] <- factor(unlist(lapply(dat$sam,function(x){
            xx <- sort(unique(as.numeric(x[[nm]])))
            xx[xx >= 0]
        })))
    if(shared_stockrecruitment){
        sr <- sapply(dat$sam,function(x)x$stockRecruitmentModelCode)
        srvd <- attr(pars$rec_pars,"vdim")
        if(sum(srvd) > 0)
            map0$rec_pars <- factor(unlist(lapply(seq_along(srvd), function(i){ if(srvd[i]==0) return(character(0)); paste0(sr[i],"_",seq_len(srvd[i]))})))
    }
    
    ## Create initial TMB Object
    ran <- c("logN", "logF", "missing", "shared_logFscale")
    obj <- TMB::MakeADFun(dat, pars, map0, random=ran, DLL="multiStockassessment", ...)

    ## Check for unused correlation parameters
    indx <- which(obj$gr()[names(obj$par) %in% "RE"] == 0)
    mvals <- seq_along(pars$RE)
    mvals[indx] <- NA
    map <- c(map0, list(RE = factor(mvals)))

    
    if(rm.unidentified){
        gr <- obj$gr()
        safemap <- obj$env$parList(gr)
        safemap <- safemap[!names(safemap)%in%ran]
        safemap <- lapply(safemap, function(x)factor(ifelse(abs(x)>1.0e-15,1:length(x),NA)))
        map <- c(map,safemap)        
    }
    
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

    p0 <- obj$par
    for(nn in names(starting)) p0[names(obj$par)==nn]=starting[[nn]]
    
    opt <- stats::nlminb(p0, obj$fn, obj$gr,
                         control=nlminb.control,
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
    opt$he <- stats::optimHess(opt$par, obj$fn, obj$gr)
    ## Get report and sdreport
    rep <- obj$report(obj$env$last.par.best)
    sdrep <- TMB::sdreport(obj,opt$par, opt$he)
    ssdrep <- summary(sdrep)

    ## Do as in stockassessment package
    # Last two states
    idx <- c(grep("_lastLogN$",names(sdrep$value)), grep("_lastLogF$",names(sdrep$value)))
    sdrep$estY <- sdrep$value[idx]
    sdrep$covY <- sdrep$cov[idx,idx]

    idx <- c(grep("_beforeLastLogN$",names(sdrep$value)), grep("_beforeLastLogF$",names(sdrep$value)))
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
    attr(res,"lower") <- lower2
    attr(res,"upper") <- upper2
    attr(res,"newtonsteps") <- newtonsteps
    attr(res,"nlminb.control") <- nlminb.control
    attr(res,"dotargs") <- list(...)
    corpars <- ssdrep[rownames(ssdrep)=="RE",]
    if(!is.null(obj$env$map$RE)){
        corpars <- corpars[obj$env$map$RE,]
        ##corpars[is.na(corpars)] <- 0
    }
    rownames(corpars) <- colnames(dat$X)
    attr(res,"corParameters") <- corpars
    class(res) <- c("msam","samset")
    return(res)
}
