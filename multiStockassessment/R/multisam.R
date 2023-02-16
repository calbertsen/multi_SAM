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
##' @importFrom TMB MakeADFun sdreport runSymbolicAnalysis
##' @references
##' Albertsen, C. M., Nielsen, A. and Thygesen, U. H. (2018) Connecting single-stock assessment models through correlated survival. ICES Journal of Marine Science, 75(1), 235-244. doi: 10.1093/icesjms/fsx114
##' @export
multisam.fit <- function(x,
                         formula = ~-1,
                         corStructure = suggestCorStructure(x,nAgeClose = 0),
                         usePartialCors=TRUE,
                         newtonsteps=0,
                         rm.unidentified=FALSE,
                         nlminb.control = list(trace=1, eval.max=20000, iter.max=20000),
                         lower=NULL,
                         upper=NULL,
                         starting=NULL,
                         community_formula = ~-1,
                         community_type = 1,
                         shared_data = NULL,
                         ## shared_conf = NULL,
                         shared_keys = character(0),
                         shared_selectivity = 0,
                         shared_seasonality = 0,
                         shared_stockrecruitment = FALSE,
                         shared_oneFScalePars = FALSE,
                         shared_initN = FALSE,
                         shared_fleetParameters = c(),
                         shared_proportionalHazard = NULL,
                         shared_phmap = NULL,
                         skip_stock_observations = FALSE,
                         stockAreas = matrix(1,1,length(x)),
                         genetics_data = prepareGenetics(),
                         genetics_dirichlet = FALSE,
                         genetics_spatialAge = TRUE,
                         genetics_independentStocks = TRUE,
                         initN = 0,
                         initF = FALSE,
                         parlist = NULL,
                         run = TRUE,
                         symbolicAnalysis = FALSE,
                         fullDerived = FALSE,
                         mohn = FALSE,
                         ...){
    mc <- match.call(expand.dots = TRUE)
    envir <- parent.frame(1L)
    ## Check input
    requireNamespace("TMB")
    if(!methods::is(x,"samset"))
        stop("x must be a samset.")

    if(is.null(shared_data)){
        shared_data <- list(hasSharedObs = 0L)
    }else{
        if(is.null(shared_data$covCombine)){
            shared_data$covCombine <- rep(0, length.out = nrow(shared_data$keyFleetStock))
        }else{
            shared_data$covCombine <- rep(shared_data$covCombine, length.out = nrow(shared_data$keyFleetStock))
        }
    }
    ## if(is.null(shared_conf))
    ##     shared_conf <- list()

    ## Prepare data for TMB
    dat0 <- collect_data(x)
    if(!is.list(shared_proportionalHazard)){
            shared_proportionalHazard <- replicate(length(x), shared_proportionalHazard, simplify = FALSE)
            shared_proportionalHazard[1] <- list(NULL)
    }else{
        if(length(shared_proportionalHazard) < length(x))
            stop("When a list, the length of shared_proportionalHazard must match the length of x")
    }
    Xph <- lapply(shared_proportionalHazard, function(f){
        if(is.null(f) || !(shared_selectivity %in% c(3,4,6)))
            return(matrix(0,0,0))
        maxYearAll = as.integer(max(unlist(lapply(dat0,function(dd)dd$years))))
        minYearAll = as.integer(min(unlist(lapply(dat0,function(dd)dd$years))))
        maxAgeAll = as.integer(max(unlist(lapply(dat0,function(dd)dd$maxAge))))
        minAgeAll = as.integer(min(unlist(lapply(dat0,function(dd)dd$minAge))))
        AYgrd <- expand.grid(Age = as.numeric(minAgeAll:maxAgeAll),
                             Year = as.numeric(minYearAll:maxYearAll))
        if(inherits(f,"formula")){
            return(model.matrix(f, AYgrd))
        }
        if(is.matrix(f)){            
            if(ncol(f) > 0 && nrow(AYgrd) == nrow(f))
                return(f)
            stop("shared_proportionalHazard matrix has incorrect dimensions")
        }
        warning("Unsupported shared_proportionalHazard")
        return(matrix(0,0,0))
    })
    if(!(is.matrix(Xph[[1]]) && ncol(Xph[[1]]) == 0 && nrow(Xph[[1]]) == 0))
        stop("First proportional hazard covariate matrix must be 0x0")

    if(any(!corStructure)){
        tms <- terms(formula)
        if(length(attr(tms,"variables")) + attr(tms,"intercept")){
            warning("corStructure used with a formula without parameters. Using ~factor(Index) instead.")
            formula <- ~factor(Index)
        }
    }
    
    dat <- list(sam = dat0,
                sharedObs = shared_data,
                geneticsData = genetics_data,
                reportingLevel = as.integer(fullDerived),
                mohn = as.integer(mohn),
                usePartialCor = as.integer(usePartialCors),
                maxYearAll = as.integer(max(unlist(lapply(dat0,function(dd)dd$years)))),
                minYearAll = as.integer(min(unlist(lapply(dat0,function(dd)dd$years)))),
                maxAgeAll = as.integer(max(unlist(lapply(dat0,function(dd)dd$maxAge)))),
                minAgeAll = as.integer(min(unlist(lapply(dat0,function(dd)dd$minAge)))),
                cons = boolMat2ConstraintList(corStructure),
                X = build_cov_design(formula, x),
                XCon = build_cov_design(community_formula,  x, type = community_type),
                ConConf = as.integer(community_type),
                Xph = Xph,
                shared_F_type = as.integer(shared_selectivity),
                shared_Fseason_type = as.integer(shared_seasonality),
                skip_stock_observations = as.integer(skip_stock_observations),
                stockAreas = stockAreas,
                fake_obs = numeric(0),
                fake_stock = integer(0),
                fake_indx = integer(0),
                fake_includeOthers = 1L,
                doResiduals = as.integer(FALSE)
                )
    ## If there are 0s for commercial fleets in shared_data$key, turn off F in conf$keyLogFsta
    ## if(any(shared_data$keyFleetStock[shared_data$fleetTypes == 0,]==0)){
    ##     ## First column of indx is fleet index, second is stock index
    ##     indx <- which(shared_data$fleetTypes[row(shared_data$keyFleetStock)] == 0 & shared_data$keyFleetStock == 0, arr.ind=TRUE)
    ##     for(ii in seq_len(nrow(indx))){
    ##         dat$sam[[indx[ii,2]]]$keyLogFsta[indx[ii,1],] <- -1
    ##         tmp <- dat$sam[[indx[ii,2]]]$keyLogFsta
    ##         tmp[tmp == -1] <- NA
    ##         tmp[] <- as.integer(factor(tmp))-1
    ##         tmp[is.na(tmp)] <- -1
    ##         dat$sam[[indx[ii,2]]]$keyLogFsta[] <- tmp
    ##         ## season
    ##         dat$sam[[indx[ii,2]]]$keyLogFseason[indx[ii,1],] <- -1
    ##         tmp <- dat$sam[[indx[ii,2]]]$keyLogFseason
    ##         tmp[tmp == -1] <- NA
    ##         tmp[] <- as.integer(factor(tmp))-1
    ##         tmp[is.na(tmp)] <- -1
    ##         dat$sam[[indx[ii,2]]]$keyLogFseason <- tmp
    ##         ## var
    ##         dat$sam[[indx[ii,2]]]$keyVarF[indx[ii,1],] <- -1
    ##         tmp <- dat$sam[[indx[ii,2]]]$keyVarF
    ##         tmp[tmp == -1] <- NA
    ##         tmp[] <- as.integer(factor(tmp))-1
    ##         tmp[is.na(tmp)] <- -1
    ##         dat$sam[[indx[ii,2]]]$keyVarF <- tmp
    ##     }
    ## }
    ## Prepare parameters for TMB
    pars <- collect_pars(x)
    pars$RE <- numeric(ncol(dat$X)) ##rep(0,sum(lower.tri(corStructure)))
    pars$betaCon <- numeric(ncol(dat$XCon))

    ## pars$shared_logSdObs <- numeric(length(dat$sharedObs$fleetTypes))
    
    if(shared_selectivity > 0){
        if(!(length(unique(attr(pars$logF,"cdim"))) == 1 &&
             length(unique(attr(pars$logF,"rdim"))) == 1))
            stop("logF must have same dimension for all stocks when sharing selectivity")
        if(shared_selectivity == 1){    #Scale by vector AR1
            lfsRDim <- attr(pars$logF,"rdim") 
            lfsCDim <- attr(pars$logF,"cdim")
            lfsCDim[1] <- 0
           initFdim <- attr(pars$logF,"rdim") * 0
        }else if(shared_selectivity == 3){                     #Scale by pure parametric
            lfsRDim <- attr(pars$logF,"cdim") * 0
            lfsCDim <- attr(pars$logF,"cdim") * 0
            initFdim <- attr(pars$logF,"rdim") * 0
        }else{                          #Scale by scalar AR/RW (+ parametric)
            lfsRDim <- attr(pars$logF,"rdim") * 0 + 1
            lfsCDim <- attr(pars$logF,"cdim")            
            lfsCDim[1] <- 0
            initFdim <- local({xx <- attr(pars$logF,"rdim") * (initF>0);  xx[-1] <- 1 * (initF>0); xx})
        }
    }else{
        lfsRDim <- lfsCDim <- attr(pars$logF,"cdim") * 0
        initFdim <- attr(pars$logF,"rdim") * 0
    }
    pars$shared_logFscale <- combineMatrices(lapply(seq_along(lfsRDim),function(i) matrix(0,lfsRDim[i],lfsCDim[i])))
    pars$shared_lfsMean <- matrix(0, nrow = lfsRDim,
                                  ncol = length(dat$sam)-1)
    pars$shared_lfsSd <- numeric(ifelse(shared_selectivity %in% c(1,2,4,5,6),length(dat$sam)-1,0))
    pars$shared_lfsRho <- numeric(ifelse(shared_selectivity %in% c(1,2,4),length(dat$sam)-1,0))    
    pars$shared_logitMissingVulnerability <- numeric(sum(is.na(dat$sharedObs$keyFleetStock)))
    pars$shared_missingObs <- numeric(sum(is.na(dat$sharedObs$logobs)))
    pars$shared_phbeta <- combineVectors(lapply(dat$Xph, function(x) rnorm(ncol(x),0,0.01)))
    
    if(initN == 0){
        pars$initLogN <- combineVectors(lapply(attr(pars$logN,"rdim") * 0,numeric))
    }else if(initN == 1){
        pars$initLogN <- combineVectors(lapply(attr(pars$logN,"rdim") * 1,numeric))
    }else if(initN == 2){
        pars$initLogN <- combineVectors(lapply(rep(1,length(dat$sam)),numeric)) #+ c(11,11,11)
    }
    ##pars$initLogN <- combineVectors(lapply(attr(pars$logN,"rdim") * ifelse(initN,1,0),numeric))
    pars$initLogF <- combineVectors(lapply(initFdim,numeric))

    if(!is.null(parlist)){
        ## Message about wrong names in parlist
        for(nm in intersect(names(pars),names(parlist)))
            if(length(parlist[[nm]]) == length(pars[[nm]])){
                pars[[nm]][] <- parlist[[nm]][]
            }else{
                ## Warning/message, does not match
            }
    }

    ## Genetics parameters
    nStockM <- length(x)
    na2 <- function(x,a) ifelse(is.na(x) | is.null(x), a, x)
    nStockG <- na2(ncol(dat$geneticsData$stock2gen),0)
    pars$gen_logKappaSpace <- 0
    pars$gen_logKappaTime <- 0
    n2cor <- function(n) numeric(n*(n-1)/2)
    pars$gen_corparST <- n2cor(pmax(0,nStockG-1))
    pars$gen_logSdST <- numeric(pmax(0,nStockG-1))
    pars$gen_corparAge <- 0
    pars$gen_corparTrip <- n2cor(pmax(0,nStockG-1))
    pars$gen_logSdTrip <- numeric(pmax(0,nStockG-1))
    ad <- attr(dat$geneticsData,"alleleDim")
    pars$gen_alleleFreq <- array(0, dim = c(pmax(ad[1]-1,0),ad[2],nStockG))
    if(length(dat$geneticsData) > 0){   # Get better starting values
        gen_samples <- dat$geneticsData$samples
        gStock <- sapply(gen_samples,function(x)x$keyStock)
        isBaseline <- !is.na(gStock)
        if(any(isBaseline)){
            pars$gen_alleleFreq <- simplify2array(lapply(lapply(split(lapply(gen_samples[isBaseline],function(x)x$alleleCount),gStock[isBaseline]), Reduce, f="+"),function(x){
                x <- x+1
                logP <- log(x / colSums(x)[col(x)])
                logP[-nrow(x),,drop=FALSE] - matrix(logP[nrow(x),],nrow(x)-1,ncol(x), byrow=TRUE)
            }))
        }
    }
    pars$gen_dmScale <- matrix(0, nrow = ad[2], ncol = nStockG)
    pars$gen_muLogP <- matrix(0,dat$maxAgeAll-dat$minAgeAll+1,nStockG)

    maxAge <- max(sapply(x, function(xx) xx$conf$maxAge))
    minAge <- max(sapply(x, function(xx) xx$conf$minAge))
    pars$logGst <- array(0, dim = c(nrow(dat$geneticsData$Qspace),
                                    nrow(dat$geneticsData$Qtime),
                                    ifelse(genetics_spatialAge,maxAge-minAge+1,1),
                                    pmax(0,nStockG-1)))
    
    pars$logGtrip <- matrix(0, pmax(0,nStockG-1), attr(dat$geneticsData,"nTrips"))
    pars$logitArea <- combineParameter(lapply(colSums(stockAreas), function(n) numeric(n-1)))
                                 
    ## Prepare map for TMB
    map0 <- collect_maps(x)

    if(shared_oneFScalePars){
        map0$shared_lfsSd <- factor(rep(1,length(pars$shared_lfsSd)))
        map0$shared_lfsRho <- factor(rep(1,length(pars$shared_lfsRho)))
    }
    if(shared_initN && initN == 2){
        map0$initLogN <- factor(rep(1,length(pars$initLogN)))
    }else if(initN == 2){
        ## map0$initLogN <- factor(c(1,NA,NA))#rep(NA,3))
    }

    if(length(dat$geneticsData$samples) == 0){
        map0$gen_logKappaSpace <- factor(NA * pars$gen_logKappaSpace)
        map0$gen_logKappaTime <- factor(NA * pars$gen_logKappaTime)
        map0$gen_corparST <- factor(NA * pars$gen_corparST)
        map0$gen_logSdST <- factor(NA * pars$gen_logSdST)
        map0$gen_corparAge <- factor(NA * pars$gen_corparAge)
        map0$gen_corparTrip <- factor(NA * pars$gen_corparTrip)
        map0$gen_logSdTrip <- factor(NA * pars$gen_logSdTrip)
        map0$gen_alleleFreq <- factor(NA * pars$gen_alleleFreq)
        map0$gen_dmScale <- factor(NA * pars$gen_dmScale)
        map0$gen_muLogP <- factor(NA * pars$gen_muLogP)
    }else{
        gStock <- sapply(gen_samples,function(x)x$keyStock)
        isBaseline <- !is.na(gStock)
        if(all(isBaseline)){
            map0$gen_logKappaSpace <- factor(NA * pars$gen_logKappaSpace)
            map0$gen_logKappaTime <- factor(NA * pars$gen_logKappaTime)
            map0$gen_corparST <- factor(NA * pars$gen_corparST)
            map0$gen_logSdST <- factor(NA * pars$gen_logSdST)
            map0$gen_corparAge <- factor(NA * pars$gen_corparAge)
            map0$gen_corparTrip <- factor(NA * pars$gen_corparTrip)
            map0$gen_logSdTrip <- factor(NA * pars$gen_logSdTrip)
            map0$gen_muLogP <- factor(NA * pars$gen_muLogP)
        }else{
            gAge <- sapply(gen_samples[isBaseline],function(x)x$keyStock)
            if(all(is.na(gAge)) || !genetics_spatialAge){
                map0$gen_corparAge <- factor(NA * pars$gen_corparAge)                
            }
            ## Is this one parameter too much?
            gmlpMap <- matrix(seq_along(pars$gen_muLogP), nrow(pars$gen_muLogP), ncol(pars$gen_muLogP))
            gmlpMap[,colSums(dat$geneticsData$stock2gen) > 0] <- NA
            map0$gen_muLogP <- factor(gmlpMap)
        }
        if(all(is.na(sapply(gen_samples,function(x)x$keyTrip)))){
            map0$gen_corparTrip <- factor(NA * pars$gen_corparTrip)
            map0$gen_logSdTrip <- factor(NA * pars$gen_logSdTrip)
        }
    }
    if(genetics_dirichlet){
        map0$gen_dmScale <- factor(col(pars$gen_dmScale))
    }else{
        map0$gen_dmScale <- factor(NA * pars$gen_dmScale)
        pars$gen_dmScale[] <- Inf
    }
    if(genetics_independentStocks)
        map0$gen_corparST <- factor(NA * pars$gen_corparST)            

    if(!is.null(shared_phmap)){
        map0$shared_phbeta <- shared_phmap
    }else if(shared_selectivity == 4){        #Both parametric scaling and AR1 level
        ## Need to fix one parameter from parametric scaling since the level/intercept is modelled elsewhere
        vdphb <- attr(pars$shared_phbeta,"vdim")
        mtmp <- split(seq_along(pars$shared_phbeta), factor(rep(seq_along(vdphb), times = vdphb),seq_along(vdphb)))
        mtmp[vdphb > 0] <- lapply(mtmp[vdphb > 0], function(xx){ xx[1] <- NA; xx})            
        map0$shared_phbeta <- factor(multiStockassessment:::combineVectors(mtmp))
    }else if(!(shared_selectivity %in% c(3,4))){
        map0$shared_phbeta <- factor(NA * pars$shared_phbeta)
    }
    
    ## Share parameters
    ## keyName <- c("keyLogFpar","keyQpow","keyVarF","keyVarLogN","keyVarObs","keyCorObs")
    ## parName <- c("logFpar","logQpow","logSdLogFsta","logSdLogN","logSdLogObs","transfIRARdist")
    sn <- getStockNames(x)
    ## Share fleet parameters
    ## if(length(shared_fleetParameters) > 0){
    ##     for(f in shared_fleetParameters){
    ##         ft <- unique(sapply(x, function(y) y$data$fleetType[f]))
    ##         if(length(ft) > 1)
    ##             stop("Fleet type must be the same for all stocks to share parameters")
    ##     }
    ##     ## keyLogFsta (F standard deviation; should not be shared)
    ##     for(nm in c("keyLogFpar","keyQpow","keyVarObs","keyCorObs")){
    ##         lapply(seq_along(dat$sam), function(ii){
    ##             v <- dat$sam[[ii]][[nm]]
    ##             v <- paste(sn[ii],v[-shared_fleetParameters,][],sep="_")
                
    ##         }
    ##         map0[[parName[match(nm,keyName)]]] <- factor(unlist(lapply(dat$sam,function(x){
    ##             xx <- sort(unique(as.numeric(x[[nm]])))
    ##             xx[xx >= 0]
    ##         })))
    ##         }
    ##     }
    ## }    
    ## Share keys
    keyName <- c("keyLogFpar","keyQpow","keyVarF","keyVarLogN","keyVarObs","keyCorObs",
                 "keyStockWeightMean","keyStockWeightObsVar","phiSW","procSdSW",
                 "keyCatchWeightMean","keyCatchWeightObsVar","phiCW","procSdCW",
                 "keyMatureMean","sdMO","phiMO","procSdMO",
                 "keyMortalityMean","keyMortalityObsVar","phiNM","procSdNM",
                 "keyLogFmu","keyLogFrho","seasonMu","seasonRho","seasonSd",
                 "sigmaObsParUS",
                 "rec_transphi")
    parName <- c("logFpar","logQpow","logSdLogFsta","logSdLogN","logSdLogObs","transfIRARdist",
                 "meanLogSW","logSdLogSW","logPhiSW","logSdProcLogSW",
                 "meanLogCW","logSdLogCW","logPhiCW","logSdProcLogCW",
                 "meanLogitMO","logSdMO","logPhiMO","logSdProcLogitMO",
                 "meanLogNM","logSdLogNM","logPhiNM","logSdProcLogNM",
                 "muF","trans_rho_F","seasonMu","seasonLogitRho","seasonLogSd",
                 "sigmaObsParUS",
                 "rec_transphi")
    if(any(!(shared_keys %in% keyName)))
        stop(sprintf("shared keys not valid: %s",paste(shared_keys[!(shared_keys %in% keyName)],collapse=", ")))
    for(nm in shared_keys)
        if(!is.na(match(nm,names(dat$sam[[1]])))){
            map0[[parName[match(nm,keyName)]]] <- factor(unlist(lapply(dat$sam,function(x){
                xx <- sort(unique(as.numeric(x[[nm]])))
                xx[xx >= 0]
            })))
        }else{            
            pnm <- parName[match(nm,keyName)]
            mTmp <- factor(unlist(lapply(splitParameter(pars[[pnm]]),seq_along)))
            if(parName[match(nm,keyName)] %in% names(map0)){
                mTmp[is.na(map0[[parName[match(nm,keyName)]]])] <- NA
            }
            map0[[parName[match(nm,keyName)]]] <- mTmp
            
        }
    if(shared_stockrecruitment){
        sr <- sapply(dat$sam,function(x)x$stockRecruitmentModelCode)
        srvd <- attr(pars$rec_pars,"vdim")
        if(sum(srvd) > 0)
            map0$rec_pars <- factor(unlist(lapply(seq_along(srvd), function(i){ if(srvd[i]==0) return(character(0)); paste0(sr[i],"_",seq_len(srvd[i]))})))
    }
    if(shared_selectivity != 0){
        lfm0 <- lapply(splitParameter(pars$logF),seq_along)
        lfm0[-1] <- lapply(lfm0[-1],function(x) x*NA)
        map0$logF <- factor(unlist(lfm0))
        ## trans_rho should all be the same
        itr0 <- lapply(splitParameter(pars$itrans_rho),seq_along)
        map0$itrans_rho <- factor(unlist(itr0))
    }
    if(shared_seasonality != 0){
        lfm0 <- lapply(splitParameter(pars$logitFseason),seq_along)
        lfm0[-1] <- lapply(lfm0[-1],function(x) x*NA)
        map0$logitFseason <- factor(unlist(lfm0))
    }
    if(skip_stock_observations){
        lfm0 <- lapply(splitParameter(pars$missing),seq_along)
        lfm0 <- lapply(lfm0,function(x) x*NA)
        map0$missing <- factor(unlist(lfm0))
    }

    ## Map F if using shared_selectivity or proportional hazard
    ## if(any(sapply(dat$Xph,nrow) > 0) || shared_selectivity){
    ##     hasPH <- sapply(dat$Xph,nrow) > 0
    ##     hasSS <- sapply(seq_along(x), function(i) i > 1 && shared_selectivity)
    ##     Fmap <- seq_along(pars$logF)
    ##     attr(Fmap,"rdim") <- attr(pars$logF,"rdim")
    ##     attr(Fmap,"cdim") <- attr(pars$logF,"cdim")
    ##     FmapL <- splitParameter(Fmap)
    ##     FmapL[hasPH | hasSS] <- lapply(FmapL[hasPH | hasSS], function(x) rep(NA,length(x)))
    ##     map0$logF <- factor(unlist(FmapL))
    ## }
    
    ## Create initial TMB Object
    ran <- c("logN", "logF", "missing","logSW", "logCW", "logitMO", "logNM", "logP","logitFseason", "shared_logFscale","shared_missingObs", "logGst", "logGtrip")
    ## prf <- c("gen_alleleFreq", "gen_muLogP")
    ## if(all(sapply(pars[prf], length) == 0)){
    prf <- NULL
    ## }else if(all(sapply(pars[prf], length) == 0)){
    ##     prf <- setdiff(prf, prf[sapply(pars[prf], length) == 0])
    ## }    
    ## return(list(dat=dat,pars=pars,map=map0,random=ran))

    obj <- TMB::MakeADFun(dat, pars, map0,
                          random=ran,
                          ##profile = prf,
                          ##regexp = FALSE,
                          DLL="multiStockassessment", ...)
    if(!run)
        return(obj)
 
    ## Check for unused correlation parameters
    if(any(names(obj$par) %in% "RE")){
        indx <- which(obj$gr()[names(obj$par) %in% "RE"] == 0)
        mvals <- seq_along(pars$RE)
        mvals[indx] <- NA
        map <- c(map0, list(RE = factor(mvals)))
    }else{
        map <- map0
    }
    
    if(rm.unidentified){
        obj$fn(obj$par + rnorm(length(obj$par),0,0.01))
        gr <- obj$gr()
        plForMap <- obj$env$parList(gr)
        safemap <- plForMap
        safemap <- safemap[!names(safemap)%in%ran]
        safemap <- lapply(names(safemap), function(nm){
            x <- safemap[[nm]]            
            vv <- 1:length(x)
            if(nm %in% names(map))
                vv <- map[[nm]]            
            factor(ifelse(abs(x)>1.0e-15,vv,NA))
        })
        names(safemap) <- names(plForMap)[!names(plForMap)%in%ran]
        map <- safemap #c(map[setdiff(names(map),names(safemap))],safemap)        
    }
   
    ## Create TMB Object
    obj <- TMB::MakeADFun(dat, pars, map,
                          random=ran,
                          ##profile = prf,
                          DLL="multiStockassessment", ...)
    if(symbolicAnalysis)
        TMB::runSymbolicAnalysis(obj)
    pars2 <- pars
    pars2$RE <- obj$par[names(obj$par)=="RE"]
    ## Fit object
   if(is.null(lower))
        lower <- getLowerBounds(pars2)
    if(is.null(upper))
        upper <- getUpperBounds(pars2)

    lower2 <- rep(-Inf,length(obj$par))
    upper2 <- rep(Inf,length(obj$par))
    for(nn in intersect(names(lower),unique(names(obj$par))))
        lower2[names(obj$par)==nn]=lower[[nn]]
    for(nn in intersect(names(upper),unique(names(obj$par))))
        upper2[names(obj$par)==nn]=upper[[nn]]

    p0 <- obj$par
    for(nn in intersect(names(starting),unique(names(obj$par))))
        p0[names(obj$par)==nn]=starting[[nn]]
    
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
    ##opt$he <- numDeriv::hessian(opt$par, obj$fn, obj$gr)
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

    ## covRecPars
    idx <- grep("_rec_pars$",names(sdrep$value))
    sdrep$covRecPars <- sdrep$cov[idx,idx, drop = FALSE]
    colnames(sdrep$covRecPars) <- rownames(sdrep$covRecPars) <- names(sdrep$value)[idx]

    ## parList    
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
    attr(res,"m_call") <- mc
    attr(res,"m_envir") <- envir
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
